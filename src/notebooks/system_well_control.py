from pymf6.mf6 import MF6
import os

dir = os.getcwd()
model_name = "gwf_transbase"
model_path = os.path.abspath(os.path.join(dir, 'models', model_name))
nam_file = os.path.join(model_path, 'mfsim.nam')

def run_model(nam_file):
    """Control script to regulate the pumping wells dynamically. The regulation of the system is regulated by
     river vital water level and the groundwater mimimum threshold. """
    print('started')
    # Initialize the MF6 model using the provided nam file
    mf6 = MF6(sim_path=model_path)
    print('mf6 initialized')

    # Get the flow models
    flow_models = mf6.models['gwf6']
    gwf = flow_models['transbase']
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

    for model in mf6.model_loop():
        # Flow model operations
        if model.model_type == 'gwf6' and model.name == 'gwf_transbase':
            wel_pkg = model.packages.wel_0.as_mutable_bc()
            well_node = wel_pkg.nodelist[0]  # Store well node
            
            if model.kper == 2:  # Target stress period
                current_head = model.X[well_node]
                current_q = wel_pkg.q[0]  # Current pumping rate
                
                # Groundwater head control
                if current_head <= lower_limit_gw:
                    if not been_below_gw:
                        current_q *= 0.9  # Reduce pumping
                        been_below_gw = True
                elif been_below_gw and current_head >= upper_limit_gw:
                    current_q *= 1.1  # Increase pumping
                    been_below_gw = False
                
                # Concentration control (using previous timestep value)
                if prev_conc >= upper_limit_conc:
                    if not been_above_conc:
                        current_q *= 0.9  # Reduce pumping
                        been_above_conc = True
                elif been_above_conc and prev_conc <= lower_limit_conc:
                    current_q *= 1.1  # Increase pumping
                    been_above_conc = False
                
                # Update well properties
                wel_pkg.q = [current_q]
                mywell_q['step'].append(model.kstp)
                mywell_q['q'].append(current_q)
                print(f"kper={model.kper}, kstp={model.kstp}: "
                      f"Q={current_q:.4f}, Head={current_head:.4f}, "
                      f"Prev_Conc={prev_conc:.4f}")
        
        # Transport model operations
        elif model.model_type == 'gwt6' and model.name == 'gwt_transbase':
            if well_node is not None:
                prev_conc = model.X[well_node]  # Update concentration

    return mywell_q
    # Plot well pumpin 
    plt.figure(figsize=(8, 5))
    plt.plot(results['step'], results['q'], marker='x')
    plt.title("Pumping Rate Regulation Over Time")
    plt.xlabel("Timestep")
    plt.ylabel("Pumping Rate [L³/T]")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    model_path = os.path.join(os.getcwd(), 'models', 'gwf_riverbase')
    nam_file = os.path.join(model_path, 'mfsim.nam')
    results = run_model(nam_file)
    print("Simulation completed. Well rates:")
    for step, q in zip(results['step'], results['q']):
        print(f"  Step {step}: Q = {q:.6f}")