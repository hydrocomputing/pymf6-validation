from pymf6.mf6 import MF6
import os
import pandas as pd

from matplotlib import pyplot as plt

def run_model(model_path, verbose=False):
    """Control script to regulate the pumping wells dynamically. The regulation of the system is regulated by
     river vital water level and the groundwater mimimum threshold. """
    print('started')
    # Initialize the MF6 model using the provided nam file
    mf6 = MF6(sim_path=model_path)
    print('mf6 initialized')

    # Get the flow models
    flow_models = mf6.models['gwf6']
    gwf = flow_models['gwf_transbase'] # Flow model name 
    transport_models = mf6.models['gwt6']
    gwt = transport_models['gwt_transbase'] # Transport model name

    # Head control parameters
    tolerance_gw = 0.01
    gw_head_limit = 0.5
    lower_limit_gw = gw_head_limit - tolerance_gw
    upper_limit_gw = gw_head_limit + tolerance_gw

    # Concentration control parameters
    tolerance_conc = 0.05
    conc_limit = 0.5
    lower_limit_conc = conc_limit - tolerance_conc
    upper_limit_conc = conc_limit + tolerance_conc

    # State tracking variables
    been_below_gw = False
    been_above_conc = False
    well_node = None
    prev_conc = 0.0  # Stores last known concentration

    # Get well package
    for _ in mf6.model_loop():
        if gwf.kper > 0: # break after 
            break

    wel = gwf.packages.get_package('wel-1').as_mutable_bc()
    well_coords = wel.nodelist[:]
    well_node_obs_coords = wel.nodelist[0]
    well_regulated_1_coords =  wel.nodelist[1]
    well_regulated_2_coords =  wel.nodelist[2]

   # print("Regulated wells:", well_regulated)
    mywell_q = {
        'step': [],
        'head': [],
        'conc': [],
        'q_well1': [],
        'q_well2': []
    }  # Dict of lists per well

    # Run the model loop
    for model in mf6.model_loop():
        if gwf.kper == 2:  # Only operate during stress period 2
            current_head = gwf.X[well_node_obs_coords]
            current_conc = gwt.X[well_node_obs_coords]
            print (current_conc)
            
            # Record system state
            mywell_q['step'].append(gwf.kstp)
            mywell_q['conc'].append(current_conc)
            mywell_q['head'].append(current_head)
            mywell_q['q_well1'].append(wel.q[1])
            mywell_q['q_well2'].append(wel.q[2])

            # Head regulation well 1 
            if current_head <= lower_limit_gw:
                been_below_gw = True
                wel.q[1] *= 0.7
                wel.q[2] *= 0.3
                #wel.q[2] *= 0.9
            elif been_below_gw and current_head >= upper_limit_gw:
                wel.q[1] *= 1.1
                wel.q[2] *= 1.1
                #wel.q[2] *= 1.1

            # Concentration regulation
            if current_conc >= upper_limit_conc:
                been_above_conc = True
                print(wel.q)
                print(wel.q[1])
                print(wel.q[2])
                #wel.q *= 0.9
                wel.q[1] *= 0.7
                wel.q[2] *= 0.3
            elif been_above_conc and current_conc <= lower_limit_conc:
                been_above_conc = False
                wel.q[1] *= 1.1
                wel.q[2] *= 1.1
                
    # Save results
    df = pd.DataFrame(mywell_q)
    df.to_csv("well_control_results.csv", index=False)
    print("Simulation completed")
    return mywell_q


def plot(mywell_q):
    """Simple visualization of results"""
    plt.figure(figsize=(10, 6))
    
    # Plot pumping rates
    plt.plot(mywell_q['step'], mywell_q['q_well1'], 'b-o', label='Well 1 Pumping')
    plt.plot(mywell_q['step'], mywell_q['q_well2'], 'r-o', label='Well 2 Pumping')
    plt.ylabel("Pumping Rate")
    plt.xlabel("Timestep")
    plt.grid(True)
    plt.legend()
    
    # Plot head and concentration
    plt.twinx()
    plt.plot(mywell_q['step'], mywell_q['head'], 'g--', label='Head')
    plt.plot(mywell_q['step'], mywell_q['conc'], 'm--', label='Concentration')
    plt.ylabel("Head/Concentration")
    plt.legend()
    
    plt.title("Well Control Performance")
    plt.tight_layout()
    plt.savefig("well_control_plot.png")
    plt.show()

if __name__ == '__main__':
    model_path = os.path.join(os.getcwd(), 'models', 'transbase')
    results = run_model(model_path=model_path, verbose=True)
    plot(results)