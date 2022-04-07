"""
Test Jupyter notebooks in the docs.
"""
# stdlib
import logging
from pathlib import Path
import time

# deps
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
import pytest

# package
from watertap.edb import ElectrolyteDB

# Logging setup

_log = logging.getLogger(__name__)

# Constants

JUPYTER_NB_VERSION = 4
EXECUTE_TIMEOUT_SEC = 600

# Exceptions


class NotebookError(Exception):
    pass


# Utility
# =======


def run_notebook(path: Path):
    # parse
    _log.debug(f"begin: parsing '{path}'")
    with path.open("r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=JUPYTER_NB_VERSION)
    _log.debug(f"end: parsing '{path}'")

    # execute
    pp = ExecutePreprocessor(timeout=EXECUTE_TIMEOUT_SEC)
    _log.debug(f"begin: executing '{path}'")
    t0 = time.time()
    try:
        metadata = {"metadata": {"path": path.parent}}
        pp.preprocess(nb, metadata)
    except (CellExecutionError, NameError) as err:
        _log.warning(f"Failed: {path} :: {err}")
        raise NotebookError(f"execution error for '{path}': {err}")
    except TimeoutError as err:
        _log.warning(f"Timeout: {path} :: {err}")
        dur, timeout = time.time() - t0, EXECUTE_TIMEOUT_SEC
        raise NotebookError(f"timeout for '{path}': {dur}s > {timeout}s")
    dur = time.time() - t0
    _log.debug(f"end: executing '{path}' duration={dur}s")


def find_notebooks(path: Path):
    return path.glob("**/*.ipynb")


@pytest.fixture(scope="session")
def docs_root(dirname="docs"):
    p = Path(__file__).parent
    while not (p / dirname).exists():
        pp = p.parent
        if pp == p:
            raise RuntimeError(f"Could not find '{dirname}' dir")
        p = pp
    yield p / dirname


def mongodb():
    try:
        edb = ElectrolyteDB()  # will fail if no DB
        edb.get_base()  # will fail if DB is not loaded
    except Exception as err:
        _log.warning(f"Could not connect to MongoDB: {err}")
    return False


# Tests
# =====


@pytest.mark.skipif(not mongodb(), reason="MongoDB is required")
@pytest.mark.component
def test_edb_notebooks(docs_root):
    print("\n")
    for nb_path in find_notebooks(docs_root / "examples" / "edb"):
        if ".ipynb_checkpoints" in nb_path.parts:
            continue
        print(f"run notebook: {nb_path}")
        run_notebook(nb_path)
