import os
from pathlib import Path

from pymf6_tools.make_model import make_input
from pymf6_tools.base_model import make_model_data

from pytest_utils import get_full_model_path, rmtree

def do_test(specific_model_data, model_path):
    try:
        rmtree(model_path)
    except FileNotFoundError:
        pass
    model_data = make_model_data(specific_model_data)
    make_input(model_data)
    found_files = set(path.name for path in model_path.glob('*'))
    found_files.remove('.internal')
    model_file_names = set((model_path / '.internal' / 'model_files'
                        ).read_text().split('\n'))
    assert model_file_names == found_files


def test_base_flow():
    model_path = get_full_model_path('flow_base')
    print(model_path)
    specific_model_data = {
        'model_path': model_path,
        'name': 'flowbase',
        'transport': False,
        'river_active': False,
        'wells_active': False,
        }
    do_test(specific_model_data, model_path)

def test_base_transport():
    model_path = get_full_model_path('transport_base')
    specific_model_data = {
        'model_path': model_path,
        'name': 'transbase',
        'transport': True,
        'river_active': False,
        'wells_active': True,
        }
    do_test(specific_model_data, model_path)

def test_base_river():
    model_path = get_full_model_path('riverbase')
    specific_model_data = {
        'model_path': model_path,
        'name': 'riverbase',
        'transport': True,
        'river_active': True,
        'wells_active': False,
        }
    do_test(specific_model_data, model_path)

def test_base_transport_river():
    model_path = get_full_model_path('transport_river_base')
    specific_model_data = {
        'model_path': model_path,
        'name': 'transriver',
        'transport': True,
        'river_active': True,
        'wells_active': True,
        }
    do_test(specific_model_data, model_path)
