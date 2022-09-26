from ..fsapi import FlowsheetInterface


def pytest_generate_tests(metafunc):
    if "fs_interface" in metafunc.fixturenames:
        by_name = dict(FlowsheetInterface.from_installed_packages())
        metafunc.parametrize(
            "fs_interface",
            list(by_name.values()),
            ids=list(by_name.keys()),
            scope="class",
        )


class TestFlowsheetInterface:
    def test_basic_type(self, fs_interface):
        assert isinstance(fs_interface, FlowsheetInterface)

    def test_build(self, fs_interface):
        fs_interface.build()
        assert fs_interface.dict()
