import pytest

@pytest.fixture
def mocker(pytestconfig):
    """Фикстура для pytest-mock"""
    return pytestconfig.pluginmanager.getplugin("mocker")