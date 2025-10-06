from pymf6.mf6 import MF6
import os
import pandas as pd

from matplotlib import pyplot as plt


def run_model(model_path, verbose=False):
    """Control script to regulate the pumping wells dynamically. The regulation of the system is regulated by
     concentration threshold at the obesrvation well, the groundwater mimimum threshold and volume of water to treat. """
    print('started')

    # Initialize the MF6 model using the provided nam file
    mf6 = MF6(model_path)
    print('mf6 initialized')

    # Get the flow models
    flow_models = mf6.models['gwf6']
    gwf = flow_models['gwf_pumptreat']  # Flow model name
    transport_models = mf6.models['gwt6']
    gwt = transport_models['gwt_pumptreat']  # Transport model name

    # Get well package
    for _ in mf6.model_loop():
        if gwf.kper > 0:  # break after
            break

    wel = gwf.packages.get_package('wel-1').as_mutable_bc()
    well_coords = wel.nodelist[:]
    well_node_obs_coords = wel.nodelist[0]
    well_regulated_1_coords = wel.nodelist[1]
    well_regulated_2_coords = wel.nodelist[2]
    initial_head = gwf.X[well_node_obs_coords]
    print(initial_head)

    # Concentration control parameters
    tolerance_conc = 0.05
    conc_limit = 0.1
    lower_limit_conc = conc_limit - tolerance_conc
    upper_limit_conc = conc_limit + tolerance_conc

    # Set limits for pumping rate
    min_rate = -0.05  # Minimum extraction rate (m³/day)
    max_rate = -50.0  # Maximum extraction rate (m³/day)

    # State tracking variables
    below_gw = False
    above_conc = False

    # print("Regulated wells:", well_regulated)
    mywell_q = {
        'step': [],
        'head': [],
        'conc': [],
        'source_conc': [],
        'q_well1': [],
        'q_well2': [],
        'q_well3': [],
        'head_state': [],
        'conc_state': [],
        'vol_water': []
    }  # Dict of lists per well

    # Models parameters
    N = 101

    # Run the model loop
    for model in mf6.model_loop():
        if gwf.kper == 1:  # Only operate during stress period 2
            # concentration at the source
            current_head = gwf.X[(0, 41, 26)]
            current_conc = gwt.X[(0, 41, 26)]
            # concentration at the observation well
            current_head = gwf.X[(0, 56, 26)]
            obs_conc = gwt.X[(0, 56, 26)]
            daily_volume = - 12 * (wel.q[0] + wel.q[1] + wel.q[2])
            # print (current_conc)

            # Record system state
            mywell_q['step'].append(gwf.kstp)
            mywell_q['conc'].append(obs_conc)
            mywell_q['source_conc'].append(current_conc)
            mywell_q['head'].append(current_head)
            current_q = wel.q.copy()  # Get current rates
            print(current_q)
            mywell_q['q_well1'].append(wel.q[0])
            mywell_q['q_well2'].append(wel.q[1])
            mywell_q['q_well3'].append(wel.q[2])
            # mywell_q['t_vol'].append(gwf.kstp * wel.q[2]*)
            mywell_q['head_state'].append('below' if below_gw else 'normal')
            mywell_q['conc_state'].append('above' if above_conc else 'normal')
            mywell_q['vol_water'].append(daily_volume)  # since lenght of each period is 12 days
            print('CONCENTRATION AT SOURCE IS', current_conc)
            print('CONCENTRATION AT OBSERVATION WELL IS', obs_conc)
            # else:

            # Concentration regulation
            if obs_conc >= upper_limit_conc:
                above_conc = True
                print(wel.q)
                q = wel.q
                # Control for well 1
                q[0] = q[0] * 1.1
                q[1] = q[1] * 1.2
                q[2] = q[2] * 1.3
                # well control
                for i in range(3):
                    if q[i] < max_rate:
                        q[i] = max_rate
                    elif q[i] > min_rate:
                        q[i] = min_rate
                wel.q = q
                print(f"Step {gwf.kstp}: Conc above limit! Increase pumping")

            elif obs_conc <= lower_limit_conc:
                above_conc = False  # reset state
                print(wel.q)
                q = wel.q
                # Control for well 1
                q[0] = q[0] * 0.9
                q[1] = q[1] * 0.7
                q[2] = q[2] * 0.8
                for i in range(3):
                    if q[i] < max_rate:
                        q[i] = max_rate
                    elif q[i] > min_rate:
                        q[i] = min_rate
                wel.q = q
                print(f"Step {gwf.kstp}: Conc recovered! Reduce pumping")

    # total volume
    total_volume = sum(mywell_q['vol_water'])  # in m³
    print(f"Total volume extracted: {total_volume:.2f} m³")

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
    plt.plot(mywell_q['step'], mywell_q['q_well2'], 'g-o', label='Well 2 Pumping')
    plt.plot(mywell_q['step'], mywell_q['q_well3'], 'r-o', label='Well 3 Pumping')
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
    plt.savefig("well_control_plot.png", dpi=300)
    plt.show()


def plot_state(mywell_q):
    """Enhanced visualization with state tracking"""
    import matplotlib as plt
    import flopy
    import numpy as np
    dir_name = os.getcwd()
    workspace = os.path.join(dir_name + '/models/pymf6', "pumptreat")
    # Set up workspace
    model_name = "pumptreat"
    gwfmodel_name = "gwf_pumptreat"
    # create gwt name
    gwtname = "gwt_pumptreat"

    # create flopy objects
    h1 = 25
    h2 = 24.5
    Nlay = 1
    N = 101
    L = 100.0
    H = 30.0
    k = 1.0
    k33 = 0.3
    q = -0.75
    times = (3000.0, 250, 1.0)
    con_max = 1000.0
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    # Plot 1: Pumping Rates
    ax1.plot(mywell_q['step'], mywell_q['q_well1'], 'b-', label='Well 1')
    ax1.plot(mywell_q['step'], mywell_q['q_well2'], 'r-', label='Well 2')
    ax1.plot(mywell_q['step'], mywell_q['q_well3'], 'g-', label='Well 3')
    ax1.set_ylabel("Pumping Rate [L³/T]")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_title("Pumping Rates")

    # Plot 2: Head with state indicators
    ax2.plot(mywell_q['step'], mywell_q['head'], 'go', label='Head')
    # ax2.axhline(y=24.75, color='gray', linestyle='--', label='Target')
    # ax2.axhline(y=25, color='red', linestyle=':', alpha=0.5, label='Lower Limit')
    # ax2.axhline(y=24, color='blue', linestyle=':', alpha=0.5, label='Upper Limit')

    # Mark head below state
    for i, state in enumerate(mywell_q['head_state']):
        if state == 'below':
            ax2.axvline(x=mywell_q['step'][i], color='orange', alpha=0.2)

    ax2.set_ylabel("Head [m]")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_title("Head at Observation Well")

    # Plot 3: Concentration with state indicators
    ax3.plot(mywell_q['step'], mywell_q['conc'], 'm-', label='Concentration')
    ax3.axhline(y=0.1, color='purple', linestyle='--', label='Target')
    ax3.axhline(y=0.05, color='pink', linestyle=':', alpha=0.7, label='Lower Limit')
    ax3.axhline(y=0.15, color='pink', linestyle=':', alpha=0.7, label='Upper Limit')

    # Mark conc above state
    for i, state in enumerate(mywell_q['conc_state']):
        if state == 'above':
            ax3.axvline(x=results['step'][i], color='red', alpha=0.2)

    ax3.set_xlabel("Timestep")
    ax3.set_ylabel("Concentration")
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    ax3.set_title("Concentration at Observation Well")

    plt.suptitle("Well Control System Performance", fontsize=16)
    plt.tight_layout()
    plt.savefig("well_control_analysis.png", dpi=300)
    plt.show()

    # plot all mess:
    # =======================================================================
    # 1. Plot Boundary Conditions
    # =======================================================================
    print("Plotting boundary conditions...")

    # Create colormap for boundary conditions
    # -1: CHD (blue), 2: WEL (red), 3: CNC (purple), 1: Active (white)
    bc_colors = ['blue', 'white', 'red', 'purple']
    bc_cmap = ListedColormap(bc_colors)

    fig, ax = plt.subplots(figsize=(10, 8))
    modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)

    # Plot ibound with custom colors
    quadmesh = modelmap.plot_ibound(ibound=ibd)

    # Add legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor='blue', label='Constant Head (CHD)'),
        Patch(facecolor='red', label='Well (WEL)'),
        Patch(facecolor='purple', label='Contamination Source (CNC)'),
        Patch(facecolor='white', label='Active Cell')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    ax.set_title("Boundary Conditions")
    plt.savefig(os.path.join(workspace, "boundary_conditions.png"))
    # plt.show()

    # =======================================================================
    # 2. Plot Head Distribution
    # =======================================================================
    print("Plotting head distribution...")

    # Load head data
    head_file = os.path.join(workspace, f"{gwfmodel_name}.hds")
    hds = flopy.utils.HeadFile(head_file)
    head = hds.get_data(kstpkper=(0, 0))  # First stress period

    fig, ax = plt.subplots(figsize=(10, 8))
    modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)

    # Plot head contours
    levels = np.linspace(head.min(), head.max(), 20)
    contours = modelmap.contour_array(head, levels=levels, colors='black')
    plt.clabel(contours, fmt="%.1f", fontsize=8)

    # Plot filled head contours
    head_plot = modelmap.plot_array(head, cmap='viridis', alpha=0.7)
    plt.colorbar(head_plot, shrink=0.5, label='Head (m)')

    # Plot boundary conditions
    modelmap.plot_ibound(ibound=ibd, alpha=0.5)

    ax.set_title("Head Distribution")
    plt.savefig(os.path.join(workspace, "head_distribution.png"))


    print("All plots saved to:", workspace)

    # =======================================================================
    # POST-PROCESSING: Concentration Time Series at (0, 56, 26)
    # =======================================================================
    print("\nGenerating concentration time series plot...")

    # Load concentration data
    conc_file = os.path.join(workspace, f"{gwtname}.ucn")
    conc_obj = flopy.utils.HeadFile(conc_file, text="CONCENTRATION")

    # Get all available times
    times = conc_obj.get_times()

    # Extract concentration at specific cell (layer, row, col) = (0, 56, 26)
    concentrations = []
    for time in times:
        conc_data = conc_obj.get_data(totim=time)
        concentrations.append(conc_data[0, 56, 26])  # Layer 0, Row 56, Column 26

    # Create plot
    plt.figure(figsize=(12, 6))
    plt.plot(times, concentrations, 'b-', linewidth=2, marker='o', markersize=5)

    # Add threshold line (example: 50 μg/L)
    threshold = 0.1
    plt.axhline(y=threshold, color='r', linestyle='--', linewidth=1.5)
    plt.text(times[-1] * 1.02, threshold, f'Threshold: {threshold}',
             color='r', va='top')

    # Format plot
    plt.title(f"Concentration at Cell (Layer 0, Row 56, Column 26)")
    plt.xlabel("Time (days)")
    plt.ylabel("Concentration (μg/L)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save and show
    plot_path = os.path.join(workspace, "concentration_time_series.png")
    plt.savefig(plot_path)
    print(f"Plot saved to: {plot_path}")
    plt.show()


if __name__ == '__main__':
    model_path = os.path.join(os.getcwd(), 'models', 'pumptreat')
    results = run_model(model_path=model_path, verbose=False)
    plot(results)
    plot_state(results)