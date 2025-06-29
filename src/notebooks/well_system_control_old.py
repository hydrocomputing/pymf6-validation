from pymf6.mf6 import MF6
from pymf6.api import States
import os
from pathlib import Path

# Change directory to the folder with the simulation files - nam
dir = os.getcwd()
model_name = "transport_river_base"
model_path = os.path.abspath(os.path.join(dir, 'models', model_name))
print(model_path)
nam_file = os.path.join(model_path, 'mfsim.nam')
nam_file

def run_model(nam_file):
    """ Control script to regulate the pumping wells dynamically.
    The regulation of the system is regulated by groundwater mimimum thresholdr abstraction and
    volume to be treated. """

    # Initialize the MF6 model using the provided nam file
    mf6 = MF6(nam_file=nam_file)
    # mf6 = MF6(nam_file=nam_file, use_modflow_api=True)
    # Set tolerance and head limit values for control - groundwater and river threshold
    tolerance = 0.01  # for the oscillation of the head for both thresholds
    gw_head_limit = 0.5  # fixed threshold according to regulations
    cnc_limit = 1.5  # fixed threshold according to regulations

    # limits for the groundwater
    lower_limit_gw = gw_head_limit - tolerance
    upper_limit_gw = gw_head_limit + tolerance

    # limits for the concentration
    lower_limit_cnc = cnc_limit - tolerance
    upper_limit_cnc = cnc_limit + tolerance

    # Variable to track if the water level has been below the groundwater limits
    been_below_gw = False
    been_below_cnc = False

    # List of wells with coordinates (layer, row, col)
    well_coords_list = [(0, 10, 4), (0, 12, 5), (0, 14, 6)]

    # list of wells to observe the behaviour of the plume
    obs_well_coords_list = [(1, 10, 20), (1, 11, 20), (1, 12, 20)]

    # Dictionary to store information about each well
    wells_info = {coords: {'been_below': False} for coords in well_coords_list}

    # Dictionary to store the info about the observations wells
    wells_info_obs = {coords: {'been_below': False} for coords in obs_well_coords_list}

    # Main simulation loop
    for sim, state in mf6.loop:
        # Check if the start of a new time step
        if state == States.timestep_start:
            # Get the model object
            ml = sim.get_model()
            # Check if stress period (kper == 2)
            if ml.kper == 2:
                # Iteration over each well
                for wel_coords, well_info in wells_info.items():
                    # Retrieve pumping rate and well head information
                    pumping = ml.wel.stress_period_data["flux"]
                    wel_head = ml.X.__getitem__(wel_coords)
                    wel_bc = ml.wel.stress_period_data
                    # Retrieve concentration at observation wells
                    cnc_obs_wells = ml.gwt.output.concentration()
                    cnc_obs = cnc_obs_wells.get_data(kstpkper=None, mflay=None)  # all layers will be included
                    # select the part of the model that needs to be protected from the contamination plume
                    cnc_obs[10:19:]
                    cnc_max_obs = max(cnc_obs)
                    cnc_min_obs = min(cnc_obs)

                    # Adjust pumping rate if the well head is below the lower limit
                    if wel_head <= lower_limit_gw:
                        wel_bc["flux"] = pumping * 0.9
                        Add
                        commentMore
                        actions
                        been_below_gw = True

                    # Adjust pumping rate if the well head is above the limit
                    elif been_below_gw and wel_head >= upper_limit_gw:
                        wel_bc["flux"] = pumping * 1.1

                    # Adjust pumping rate if the well head is below the river stage lower limit
                    if cnc_min_obs <= lower_limit_cnc:
                        wel_bc["flux"] = pumping * 0.9
                        been_below_cnc = True

                    # Adjust pumping rate if the well head is above the limit
                    elif been_below_cnc and cnc_max_obs >= upper_limit_cnc:
                        wel_bc["flux"] = pumping * 1.1


if __name__ == '__main__':
    run_model(nam_file)
    print('done')