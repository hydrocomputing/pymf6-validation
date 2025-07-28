"""Base model data.

Data for a basic flow model build with `flopy`.
Changing the data allows to quickly create a modified model.
Note: Not all packages are supported, yet.

Assumptions:

1. The first stress period is steady state.
2. Currently, only CHD boundary conditions are supported.
3. For all not supplied values `flopy` default values will be used.
"""


import sys
import numpy as np
from copy import deepcopy


BASE_MODEL_DATA = {
    #  flopy.mf6.ModflowTdis
    'times': (
        10.0,  # perlen (double) is the length of a stress period.
        120,   # nstp (integer) is the number of time steps in a stress period.
        1.0,   # tsmult (double) is the multiplier for the length of successive
               # time steps.
    ),
    'time_units': 'DAYS',
    'length_units': 'meters',
    'repeat_times': 3,  # nper = repeat_times + 1
    #  flopy.mf6.ModflowGwfdis
    'nrow': 15,
    'ncol': 10,
    'nlay': 3,
    'delr': 100.0,
    'delc': 100.0,
    'top': 15.0,
    'botm': [-5.0, -10.0, -15.0],
    #  flopy.mf6.ModflowGwfnpf
    'k': [0.5, 0.000006, 0.5],  # initial value of k
    'k33': [0.1, 0.002, 0.3],  # vertical anisotropy
    #  flopy.mf6.ModflowGwfsto
    'sy': 0.2,
    'ss': 0.000001,
    'initial_head': 10.0,
    # flopy.mf6.ModflowGwfchd(
    'chd': [
        [[(lay, 0, col), 10.] for lay in range(3) for col in range(10)],
        [[(lay, 14, col), 10.] for lay in range(3) for col in range(10)]
    ],
    'transport': False,
    'river': False,
}         

BASE_TRANSPORT_MODEL_DATA = {
    'wells':{},
    'initial_concentration': 1,
    'cnc': [
        [(0, 5, 1), 10.],
        [(0, 6, 1), 10.] # cell_id, conc (const)
    ],
    'scheme': 'UPSTREAM', #'TVD',  # or 'UPSTREAM'
    'longitudinal_dispersivity': 1.0,
    # Ratio of transverse to longitudinal dispersitivity
    'dispersivity_ratio': 1.0,
    'porosity': 0.35,
    'obs': None,
    'chd': [
        [[(lay, 0, col), 10.] for lay in range(3) for col in range(10)],
        [[(lay, 14, col), 10.] for lay in range(3) for col in range(10)]
    ],
}


NRIV = 14

BASE_RIVER_MODEL_DATA = {
    'river_spd': { 
        'rivlay': [0] * NRIV,
        'rivrow': [1, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8, 9],
        'rivcol': [0, 2, 1, 2, 3, 3, 4, 5, 6, 7, 7, 8, 8, 9],
        'rivstg': np.linspace(13, 14, NRIV), 
        'rivbot': np.linspace(7, 10, NRIV), 
        'rivcnd': [0.05] * NRIV  
        } , 
                     
    'river_boundnames': None, 
    'obs_dict': None, # dict, 
    'tsdict': None, # dict,
    'cond': None, 
}         

BASE_WELL_MODEL_DATA = {
    'wells': {
        'wel_out': {'q': (-0.05, -0.5, -0.05), 'coords': (0, 4, 4)},
        'wel_out1': {'q': (-0.05, -0.5, -0.05), 'coords': (0, 6, 4)},
        'wel_out2': {'q': (-0.05, -0.5, -0.05), 'coords': (0, 8, 4)},
              },

}

def make_model_data(
        specific_model_data,
        base_model_data=BASE_MODEL_DATA,
        base_transport_model_data=BASE_TRANSPORT_MODEL_DATA,
        base_river_model_data=BASE_RIVER_MODEL_DATA, 
        base_well_model_data=BASE_WELL_MODEL_DATA  ):
    """Make model data.

    specific_model_data - dictionary with data specific for the current model
                          will merged with `base_model_data`
                          existing keys in `base_model_data` will be overridden
    base_model_data - dictionary with basic model data defaults to
                      `BASE_MODEL_DATA`
    """
    base_model_data=deepcopy(base_model_data)
    base_river_model_data=deepcopy(base_river_model_data)
    base_transport_model_data=deepcopy(base_transport_model_data)
    base_well_model_data=deepcopy(base_well_model_data)
    
    if specific_model_data['transport']:
        base_model_data.update(base_transport_model_data)
    if specific_model_data['river_active']:
        base_model_data.update(base_river_model_data)
    if specific_model_data['wells_active']:
        base_model_data.update(base_well_model_data)
    # old way up to Python 3.8
    if sys.version_info[:2] < (3, 9):
        return {**base_model_data, **specific_model_data}
    # new way starting from Python 3.9
    return base_model_data | specific_model_data
