"""
Access to the neural network used for the NLTE treatment of Ca II 8542.

The weights are too large to go in the wheels, so they are downloaded the first
time they are needed and kept in a cache directory. On machines whose compute
nodes have no network access, call ``hazel.download_nlte_model()`` once from a
login node, or set HAZEL_DATA_DIR to a directory that holds the file.
"""
import hashlib
import os
import tempfile
import urllib.request
from pathlib import Path

__all__ = ['download_nlte_model']

FILENAME = '2024-09-13-11_17_34.best.pth'
URL = 'https://cloud.iac.es/index.php/s/FiSGwCMGyATJfKB/download'
SHA256 = '9a8b0aff4b0f174fecb94401476621033e3b976426630e376ba18a65f11c6649'


def _cache_dir():
    root = os.environ.get('XDG_CACHE_HOME') or Path.home() / '.cache'
    return Path(root) / 'hazel'


def _download(destination):
    """
    Download the weights to `destination`. The file is written under a temporary
    name and renamed only once its checksum is verified, so that a failed download
    leaves nothing behind and processes running concurrently (e.g., MPI ranks)
    never read a partial file. The checksum is essential: the file is a pickle, so
    loading it can execute code.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=destination.parent, prefix=FILENAME + '.')
    try:
        sha256 = hashlib.sha256()
        with os.fdopen(fd, 'wb') as out, urllib.request.urlopen(URL, timeout=60) as response:
            while chunk := response.read(1 << 20):
                sha256.update(chunk)
                out.write(chunk)
        if sha256.hexdigest() != SHA256:
            raise OSError(f'Checksum mismatch for the NLTE neural network downloaded from {URL}')
        os.chmod(tmp, 0o644)
        os.replace(tmp, destination)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def get_checkpoint():
    """
    Return the path to the NLTE neural network weights, downloading them if needed.

    Returns
    -------
    str
        Path to the checkpoint. The locations searched are HAZEL_DATA_DIR, the
        ``data`` directory of the package (present in a source checkout) and the
        cache directory (``~/.cache/hazel`` by default).
    """
    candidates = [Path(__file__).parent / 'data' / FILENAME, _cache_dir() / FILENAME]
    if 'HAZEL_DATA_DIR' in os.environ:
        candidates.insert(0, Path(os.environ['HAZEL_DATA_DIR']) / FILENAME)

    for path in candidates:
        if path.exists():
            return str(path)

    destination = candidates[-1]
    print(f'Downloading the NLTE neural network (76 MB) to {destination}')
    _download(destination)
    return str(destination)


def download_nlte_model():
    """
    Download the NLTE neural network now rather than on first use.

    Returns
    -------
    str
        Path to the checkpoint.
    """
    return get_checkpoint()
