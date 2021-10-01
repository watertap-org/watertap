
_MARKERS = {
    'unit': 'quick tests that do not require a solver, must run in < 2 s',
    'component': 'quick tests that may require a solver',
    'integration': 'long duration tests',
    'build': 'FIXME for building stuff?'
}


def pytest_configure(config):

    for name, descr in _MARKERS.items():
        config.addinivalue_line(
            'markers', f'{name}: {descr}'
        )
