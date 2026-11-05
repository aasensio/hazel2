I ran Part 1 steps 1–6 in a scratch environment and they work. Step 7 is optional and I didn't run it exactly as written, though I ran the same kind of build and clean-container test earlier. Parts 2–4 are on GitHub and PyPI, so I couldn't run them, and button names may differ slightly. Your default branch is master, not main.
Part 1: Test the branch on your machine
Nothing here touches your existing environment.
1. Go to the branch (gfortran comes from your conda environment):
conda activate py313
cd ~/Dropbox/GIT/hazel2
git branch --show-current     # must print: pip-wheels

2. Install into a new, separate environment:
python -m venv ~/hazel-test
source ~/hazel-test/bin/activate
pip install -U pip
pip install . pytest

3. Check you're using the installed copy. You must run this outside the repo, otherwise Python imports the local source folder and the test means 
nothing.
cd ~                                                                                                                                              
python -c "import hazel; print(hazel.__version__, hazel.__file__)"
This should print 1.1.0 and a path inside ~/hazel-test/lib/....

4. Run the tests. Use plain pytest, not python -m pytest, which would test the source folder instead.
cd ~/Dropbox/GIT/hazel2
pytest test -q                # expect: 10 passed

5. Test the weights download:
pip install torch --index-url https://download.pytorch.org/whl/cpu
cd ~
python -c "import hazel; print(hazel.download_nlte_model())"   # downloads 76 MB
python -c "import hazel; print(hazel.download_nlte_model())"   # no download this time
Both should print a path in ~/.cache/hazel/.

6. Run a real NLTE synthesis:
cp -r ~/Dropbox/GIT/hazel2/examples ~/hazel-ex
cd ~/hazel-ex/nonmpi/syn
cat > nlte_check.py <<'EOF'
import hazel
tmp = hazel.tools.File_photosphere(mode='single')
tmp.set_default(n_pixel=1, default='falc82')
tmp.save('../../photospheres/falc')
mod = hazel.Model('../../configurations/conf_caii.ini', working_mode='synthesis', verbose=0, root='../../')
mod.use_nlte = False
mod.synthesize()
lte = mod.spectrum['spec1'].stokes[0].copy()
mod.use_nlte = True
mod.atmospheres['ph1'].load_nlte_model(verbose=0)
mod.synthesize()
print('Largest NLTE-LTE difference:', abs(mod.spectrum['spec1'].stokes[0] - lte).max())
EOF
python nlte_check.py          # expect about 0.926

7. Optional: test a real wheel with no compilers. This needs Docker and about 2 GB on your root disk, which has 4.7 GB free.
deactivate
cd ~/Dropbox/GIT/hazel2
./build_wheels.sh
docker run --rm -v "$PWD/wheelhouse:/wh:ro" python:3.13-slim bash -c \
  'pip install /wh/hazel2-*cp313*.whl && cd /tmp && python -c "import hazel; print(hazel.__version__)"'
docker rmi python:3.13-slim

8. Clean up:
deactivate 2>/dev/null
rm -rf ~/hazel-test ~/hazel-ex

Part 2: Merge into master

Merging publishes nothing, and the pull request (PR) lets CI check everything first.

1. Commit and push the branch:
cd ~/Dropbox/GIT/hazel2
git status                    # no stray files such as lower, upper or test outputs
git add -A
git status                    # read the list once more
git commit -m "Build with meson-python and distribute binary wheels as hazel2"
git push -u origin pip-wheels

2. Open the PR at https://github.com/aasensio/hazel2/pull/new/pip-wheels, with master as the base.

3. Wait for CI. The PR starts the Wheels workflow: one sdist job and four wheel jobs. This is the first real test of macOS and Linux aarch64, so a
failure here is expected and harmless. Send me the failing log, and each fix pushed to the branch reruns CI. Build and test also runs on the push, and
it only has to pass for its test step; a Coveralls failure doesn't block anything.

4. Merge once the Wheels jobs are green. Use Create a merge commit, so one Revert click on the PR undoes everything.

5. Update your local copy and environment:
git switch master
git pull
git branch -d pip-wheels
pip uninstall hazel
rm hazel/codes/*.so
pip install meson-python "meson<1.11" ninja
pip install --no-build-isolation -e .

Part 3: Set up PyPI

You only do this once, and no passwords or tokens go into the repo.

1. Create accounts on https://pypi.org and https://test.pypi.org. They are separate sites, and both require two-factor authentication.

2. Add a pending publisher on each site. Go to Your account → Publishing and fill in:

┌──────────────────┬────────────┬───────────────┐
│      Field       │  pypi.org  │ test.pypi.org │
├──────────────────┼────────────┼───────────────┤
│ Project name     │ hazel2     │ hazel2        │
├──────────────────┼────────────┼───────────────┤
│ Owner            │ aasensio   │ aasensio      │
├──────────────────┼────────────┼───────────────┤
│ Repository       │ hazel2     │ hazel2        │
├──────────────────┼────────────┼───────────────┤
│ Workflow name    │ wheels.yml │ wheels.yml    │
├──────────────────┼────────────┼───────────────┤
│ Environment name │ pypi       │ testpypi      │
└──────────────────┴────────────┴───────────────┘

These must match exactly, or the upload is refused. Doing this on pypi.org also reserves the name.

3. Create the GitHub environments. In the repo, go to Settings → Environments and create testpypi and pypi. On pypi, enable Required reviewers and add yourself, leaving "Prevent self-review" off. Every real release will then wait for your approval.

Part 4: Publish

1. Do a TestPyPI dry run. Go to Actions → Wheels → Run workflow and choose master. The button only appears after the merge.

2. Install from TestPyPI, ideally on a machine without compilers:
python -m venv ~/hazel-tp && source ~/hazel-tp/bin/activate
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ hazel2
cd ~ && python -c "import hazel; print(hazel.__version__, hazel.__file__)"
deactivate && rm -rf ~/hazel-tp

3. Release to PyPI. This step is permanent, so only do it once step 2 works:
git switch master && git pull
git tag -a v1.1.0 -m "Hazel 1.1.0"
git push origin v1.1.0
The workflow rebuilds everything, then waits. Open the run, click Review deployments, and approve it. If the tag doesn't match the version in meson.build, the job stops before uploading.

4. Confirm it worked with pip install hazel2 in a fresh environment.

Later releases: change version: in meson.build, merge that, then tag the same number.

Three things you can't easily undo

- Version numbers can't be reused: this applies on both PyPI and TestPyPI, even after deleting a release. If TestPyPI needs a second upload, set the version to 1.1.0rc1 for that test, then go back to 1.1.0 for the real tag. If a bad version reaches PyPI, yank it and release 1.1.1.
- Keep the cloud link alive: every released version has that URL and checksum built in. If you delete the share, NLTE downloads break for all those versions, and replacing the file makes them fail the checksum. For new weights, upload a new file and release a new version with its URL and SHA256 in hazel/nlte_weights.py.
- set_nlte() is still broken: examples/nonmpi/syn/caii_syn.py crashes on it, which is why step 6 sets use_nlte directly. Users will hit this, so it's worth fixing before the release, and I can do that on this branch if you like.