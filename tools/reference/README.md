# Reference (oracle) data

Python scripts here generate JSON fixtures in `test/fixtures/` using Cantera 3.2.

    uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python cantera numpy scipy
    .venv/bin/python tools/reference/<script>.py

Each script writes `test/fixtures/<module>_<name>.json` and is deterministic.
