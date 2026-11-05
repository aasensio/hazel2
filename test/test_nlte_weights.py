import hashlib
import os

import pytest

from hazel import nlte_weights


@pytest.fixture
def remote(tmp_path, monkeypatch):
    """Serve a fake checkpoint from a file:// URL with an empty cache."""
    source = tmp_path / 'remote.pth'
    source.write_bytes(b'weights')
    monkeypatch.setattr(nlte_weights, 'URL', source.as_uri())
    monkeypatch.setattr(nlte_weights, 'SHA256', hashlib.sha256(b'weights').hexdigest())
    monkeypatch.setattr(nlte_weights, 'FILENAME', 'test_weights.pth')
    monkeypatch.setenv('XDG_CACHE_HOME', str(tmp_path / 'cache'))
    monkeypatch.delenv('HAZEL_DATA_DIR', raising=False)
    return source


def test_downloads_to_cache_once(remote, tmp_path):
    path = nlte_weights.get_checkpoint()
    assert path == str(tmp_path / 'cache' / 'hazel' / 'test_weights.pth')
    assert open(path, 'rb').read() == b'weights'

    remote.unlink()
    assert nlte_weights.get_checkpoint() == path, "the cached file should be reused"


def test_bad_checksum_leaves_nothing(remote, tmp_path):
    remote.write_bytes(b'tampered')
    with pytest.raises(OSError, match='Checksum'):
        nlte_weights.get_checkpoint()
    assert os.listdir(tmp_path / 'cache' / 'hazel') == []


def test_data_dir_takes_precedence(remote, tmp_path, monkeypatch):
    local = tmp_path / 'local'
    local.mkdir()
    (local / 'test_weights.pth').write_bytes(b'local')
    monkeypatch.setenv('HAZEL_DATA_DIR', str(local))
    assert nlte_weights.get_checkpoint() == str(local / 'test_weights.pth')
