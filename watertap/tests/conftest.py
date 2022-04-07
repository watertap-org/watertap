import logging
from pathlib import Path
from watertap.edb import ElectrolyteDB
import pytest

_log = logging.getLogger(__name__)

DOCS_DIR = "docs"  # expected documentation directory name


def get_docs_root():
    start = Path(__file__)
    p, pp = start, start.parent
    # loop while there is a parent and it's not named 'site_packages'
    while pp != p and pp.name != "site_packages":
        target = p / DOCS_DIR
        if target.exists() and target.is_dir():
            return start, target
        p, pp = pp, pp.parent
    # not found
    return start, None


@pytest.fixture(scope="session")
def docs_root():
    """Find docs root, or call pytest.skip"""
    start, result = get_docs_root()
    if result is None:
        pytest.skip(f"No directory '{DOCS_DIR}' found from '{start}'")
    yield result


def check_for_mongodb() -> bool:
    try:
        edb = ElectrolyteDB()  # will fail if no DB
        edb.get_base()  # will fail if DB is not loaded
    except Exception as err:
        _log.warning(f"Could not connect to MongoDB: {err}")
    return False


@pytest.fixture(scope="module")
def electrolytedb():
    """See if EDB can be instantiated, or call pytest.skip"""
    if not check_for_mongodb():
        pytest.skip("MongoDB is required")
