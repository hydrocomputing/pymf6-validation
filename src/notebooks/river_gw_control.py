from pymf6.mf6 import MF6
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from flopy.plot import PlotMapView



def run_model(model_path):
    """Control script to regulate the pumping wells dynamically. The regulation of the system is regulated by
     river vital water level and the groundwater mimimum threshold. """
    print('started')
    # Initialize the MF6 model using the provided nam file
    mf6 = MF6(sim_path=model_path)
    print('mf6 initialized')

    # Get the flow models
    flow_models = mf6.models['gwf6']
    gwf = flow_models['gwf_riverbase']

    # Set tolerance and head limit values for control - groundwater and river threshold
    tolerance = 0.01  # for the oscillation of the head for both thresholds
    gw_head_limit = 0.5  # fixed threshold according to regulations

    # limits for the groundwater
    lower_limit_gw = gw_head_limit - tolerance
    upper_limit_gw = gw_head_limit + tolerance

    # Variable to track if the water level has been below the groundwater limits
    been_below_gw = False

    # Dictionary to store information about each well
    mywell_q = {'step':[], 'q':[]}

    for _ in mf6.model_loop():
        if gwf.kper > 1:
            break

    mywell = gwf.packages.wel_0.as_mutable_bc()
    mywell_coords = mywell.nodelist[0]
    # Main simulation loop
    for model in mf6.model_loop():
        print('inside of the loop')
        # Check if the start of a new time step
        if gwf.kper == 2:
            mywell_head = gwf.X[mywell_coords]
            if mywell_head <= lower_limit_gw:
                been_below_gw = True
                mywell.q *= 0.9
            elif been_below_gw and mywell_head >= upper_limit_gw:
                mywell.q *= 1.1
            mywell_q['step'].append(gwf.kstp)
            mywell_q['q'].append(mywell.q[0])

        #print(mywell_q)

    # Save the well flow history
    df = pd.DataFrame(mywell_q)
    df.to_csv("well_q_log.csv", index=False)

    # Plot the well flow
    plt.figure(figsize=(8, 5))
    plt.plot(df['step'], df['q'], marker='x')
    plt.xlabel('Timestep')
    plt.ylabel('Well Flow Rate (q)')
    plt.title('Well Pumping Rate Over Time')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("well_flow_plot.png")
    plt.show()

if __name__ == '__main__':
    dir = os.getcwd()
    model_name = "riverbase"
    model_path = os.path.abspath(os.path.join(dir, 'models', model_name))
    run_model(model_path)
    print('done')
