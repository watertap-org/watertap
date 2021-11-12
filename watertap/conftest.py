from _pytest.config import Config
from _pytest.config.argparsing import Parser


_MARKERS = {
    "unit": "quick tests that do not require a solver, must run in < 2 s",
    "component": "quick tests that may require a solver",
    "integration": "long duration tests",
    "build": "FIXME for building stuff?",
}


def pytest_configure(config: Config):

    for name, descr in _MARKERS.items():
        config.addinivalue_line("markers", f"{name}: {descr}")


def pytest_addoption(parser: Parser):
    parser.addoption(
        "--edb-no-mock",
        help="Force the `edb` fixture to connect to a running MongoDB instance "
             "instead of falling back to mongomock",
        action="store_true",
        default=False,
        dest="edb_no_mock",
    )
