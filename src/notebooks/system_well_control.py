from pymf6.mf6 import MF6
import os

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
    gwf = flow_models['gwf_transbase']
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
    mywell_q = {'step': [], 'q': []}  # Pumping rate history

    for _ in mf6.model_loop():
        if gwf.kper > 1:
            break
    wel = gwf.packages.get_package('wel-1').as_mutable_bc()
    well_node = wel.nodelist[0]

    for model in mf6.model_loop():
        # Flow model operations
        if gwf.kper == 2:  # Target stress period
            current_head = gwf.X[well_node]

            # Groundwater head control
            if current_head <= lower_limit_gw:
                if not been_below_gw:
                    wel.q *= 0.9  # Reduce pumping
                    been_below_gw = True
            elif been_below_gw and current_head >= upper_limit_gw:
                wel.q *= 1.1  # Increase pumping
                been_below_gw = False

            # Concentration control (using previous timestep value)
            if prev_conc >= upper_limit_conc:
                wel.q *= 0.9  # Reduce pumping
                been_above_conc = True
            elif prev_conc <= lower_limit_conc and been_above_conc:
                wel.q *= 1.1  # Increase pumping
                been_above_conc = False

            # Update well properties
            mywell_q['step'].append(gwf.kstp)
            mywell_q['q'].append(wel.q)
            if verbose:
                 print(f"kper={gwf.kper}, kstp={gwf.kstp}: "
                    f"Q={wel.q:.4f}, Head={current_head:.4f}, "
                    f"Prev_Conc={prev_conc:.4f}")
        prev_conc = gwt.X[well_node]  # Update concentration
    return mywell_q

def plot(results):
    # Plot well pumping
    plt.figure(figsize=(8, 5))
    plt.plot(results['step'], results['q'], marker='x')
    plt.title("Pumping Rate Regulation Over Time")
    plt.xlabel("Timestep")
    plt.ylabel("Pumping Rate [L³/T]")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    model_path = os.path.join(os.getcwd(), 'models', 'transbase')
    results = run_model(model_path=model_path)
    plot(results=results)
    # print("Simulation completed. Well rates:")
    # for step, q in zip(results['step'], results['q']):
    #     print(f"  Step {step}: Q = {q:.6f}")