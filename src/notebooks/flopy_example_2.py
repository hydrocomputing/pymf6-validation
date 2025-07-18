

def example_2_wells():
    # imports
    import os
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import numpy as np
    import flopy

    # For this example, we will set up a temporary workspace.
    # Model input files and output files will reside here.
    dir_name = os.getcwd()
    workspace = os.path.join(dir_name + '/models', "pumptreat")

    # Set up workspace
    model_name = "pumptreat"
    gwfmodel_name = "gwf_pumptreat"

    # create flopy objects
    h1 = 25
    h2 = 24.5
    Nlay = 1
    N = 101
    L = 100.0
    H = 30.0
    k = 1.0
    k33 = 0.3
    q = -50.0
    times = (3000.0, 150, 1.0)
    con_max = 1000.0

    # Create the Flopy simulation object
    sim = flopy.mf6.MFSimulation(
        sim_name=gwfmodel_name, exe_name="C:/Users/lucialabarca/mf6.6.2_win64/bin/mf6.exe", version="mf6", sim_ws=workspace
    )

    # Create the Flopy temporal discretization object
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        sim, pname="tdis", time_units="DAYS", nper=2, perioddata=[(1.0, 1, 1.0), times]
    )

    # Create the Flopy groundwater flow (gwf) model object
    model_nam_file = f"{gwfmodel_name}.nam"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfmodel_name, model_nam_file=model_nam_file, save_flows=True)

    # Solver for GWF (register this one FIRST)
    ims_gwf = flopy.mf6.ModflowIms(
        sim,
        pname="ims_gwf",
        filename=f"{gwfmodel_name}" + ".ims",
        complexity="SIMPLE",
        print_option="SUMMARY"
    )
    sim.register_ims_package(ims_gwf, [gwf.name])  # Register only to GWF

    # Create the discretization package
    bot = np.linspace(-H / Nlay, -H, Nlay)
    delrow = delcol = L / (N - 1)
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=Nlay,
        nrow=N,
        ncol=N,
        delr=delrow,
        delc=delcol,
        top=0.0,
        botm=bot,
    )

    # Create the initial conditions package
    start = h1 * np.ones((Nlay, N, N))
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic", strt=start)

    # Node property flow - MUST HAVE pname
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        pname="npf",  # Critical parameter
        icelltype=1,
        k=k,
        k33=k33,
        save_flows=True
    )

    # Instantiating storage package

    sy = flopy.mf6.ModflowGwfsto.sy.empty(  # pylint: disable-msg=invalid-name
            gwf,
            default_value=0.2
        )
    ss = flopy.mf6.ModflowGwfsto.ss.empty(  # pylint: disable-msg=invalid-name
            gwf, default_value=0.00001
        )


    sto = flopy.mf6.ModflowGwfsto(
            gwf,
            pname='sto',
            save_flows=True,
            iconvert=1,
            ss=ss,
            sy=sy,
            steady_state={0: True},
            transient={index: True for index in range(1, len(times))},
            )

    # Create the constant head package.
    # List information is created a bit differently for
    # MODFLOW 6 than for other MODFLOW versions.  The
    # cellid (layer, row, column, for a regular grid)
    # must be entered as a tuple as the first entry.
    # Remember that these must be zero-based indices!
    chd_kwargs = {
        'auxiliary': 'CONCENTRATION',
        'pname': 'CHD-1'
        }
    chd_rec = []
    #chd_rec.append(((0, int(N), int(N)), h2))
    for layer in range(0, Nlay):
        for row_col in range(0, N):
            chd_rec.append(((layer, row_col, 0), h1, 0.0))
            chd_rec.append(((layer, row_col, N - 1), h2, 0.0))

    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        stress_period_data=chd_rec,
        save_flows=True,
        **chd_kwargs
    )

    # Create well package
    # wel location at (1,
    wel_kwargs = {
        'auxiliary': 'CONCENTRATION',
        'pname': 'WEL-1'
        }

    wel_rec = [
        ((0, int(N / 2.5), int(N / 4)), q, 0),
        ((0, int(N / 3), int(N / 4)), q, 0),
        ((0, int(N / 4), int(N / 4)), q, 0)
    ]
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wel_rec,
        save_flows=True,
        **wel_kwargs
    )

    # Create the output control package
    headfile = f"{gwfmodel_name}.hds"
    budgetfile = f"{gwfmodel_name}.bud"
    saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]

    oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
        gwf,
        pname="gwfgwt",
        budget_filerecord=budgetfile,
        head_filerecord=headfile,
        saverecord=saverecord,
        )

    # create gwt name
    gwtname = "gwt_pumptreat"
    # create transport model
    gwt = flopy.mf6.MFModel(
            sim,
            model_type='gwt6',
            modelname=gwtname,
            model_nam_file=f'{gwtname}.nam'
        )


    # Solver for GWT (register this one SECOND)
    imsgwt = flopy.mf6.ModflowIms(
            sim,
            print_option="SUMMARY",
            # outer_dvclose=hclose,
            # outer_maximum=nouter,
            under_relaxation="NONE",
            # inner_maximum=ninner,
            # inner_dvclose=hclose,
            # rcloserecord=rclose,
            linear_acceleration="BICGSTAB",
            scaling_method="NONE",
            reordering_method="NONE",
            # relaxation_factor=relax,
            filename=f'{gwtname}.ims',
            pname='ims'
        )

    sim.register_ims_package(imsgwt, [gwt.name])  # Register only to GWT

    # 2. Add DIS package to GWT (MUST match GWF dimensions)
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=Nlay,
        nrow=N,
        ncol=N,
        delr=delrow,
        delc=delcol,
        top=0.0,
        botm=bot,
        filename=f'{gwtname}.dis',
        pname='dis'
    )
    # 3. Initial concentration
    flopy.mf6.ModflowGwtic(
        gwt,
        strt=0.0,
        filename=f'{gwtname}.ic',
        pname='ic')  # Initial concentration = 0 everywhere

    # advection package
    adv = flopy.mf6.ModflowGwtadv(
        gwt,
        scheme="UPSTREAM",
        filename=f'{gwtname}.adv',
        pname='adv')

    # dispersion package
    dsp = flopy.mf6.ModflowGwtdsp(gwt,
                                alh=1.0,
                                ath1=0.5,
                                filename=f'{gwtname}.dsp',
                                pname='dsp')

    # Add MST package (Mass Storage and Transfer) — REQUIRED!
    flopy.mf6.ModflowGwtmst(gwt,
                            porosity=0.2,
                            first_order_decay=False,
                            decay=None,
                            decay_sorbed=None,
                            sorption=None,
                            bulk_density=None,
                            distcoef=None,
                            filename=f'{gwtname}.mst',
                            pname='mst'
                            )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    sourcerecarray = [
        ('CHD-1', 'AUX', 'CONCENTRATION'),
    ]
    sourcerecarray.append(('WEL-1', 'AUX', 'CONCENTRATION'))

    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=sourcerecarray,
        filename=f'{gwtname}.ssm',
        pname='ssm'
        )


    # CNC package
    cnc_rec = []
    # Add contamination sources
    cnc_rec.append(((0, int(N / 2.5), int(N / 4)), con_max))  # Single source
    cnc_rec.append(((0, int(N / 2.8), int(N / 4)), con_max))  # Second source
    cnc_rec.append(((0, int(N / 2.4), int(N / 4)), con_max))  # Single source
    cnc_rec.append(((0, int(N / 3), int(N / 4)), con_max))
    cnc_rec.append(((0, int(N / 2.4), int(N / 4.3)), con_max))  # Single source
    cnc_rec.append(((0, int(N / 3), int(N / 4.4)), con_max))


    # Create CNC package
    cnc = flopy.mf6.ModflowGwtcnc(
        gwt,
        stress_period_data=cnc_rec,
        pname="cnc"
    )

    # Output control for GWT
    gwt_oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "LAST")],   # REQUIRED to save .UCN
        printrecord=[("CONCENTRATION", "LAST")],
        filename=f"{gwtname}.oc",
    )

    # GWF-GWT exchange
    flopy.mf6.ModflowGwfgwt(sim,
                            exgtype="GWF6-GWT6",
                            exgmnamea=gwfmodel_name,
                            exgmnameb=gwtname,
                            filename=f"{model_name}.gwfgwt")

    # Write and run
    sim.write_simulation()
    success, msg = sim.run_simulation()
    print("Success:", success)


    # Create modflow 6 input files
    # Write and run
    sim.write_simulation()
    success, buff = sim.run_simulation()
    print("Success:", success)
    if not success:
        print("\n".join(buff))

        # visualize -----------------------------------------------------------------------------------
        # Create ibound array to identify boundary condition locations
    ibd = np.ones((Nlay, N, N), dtype=int)  # Start with all cells active (1)

    # Mark constant head cells with -1
    for record in chd_rec:
        cellid = record[0]
        k, i, j = cellid
        ibd[k, i, j] = -1

    # Mark contamination sources with 3
    for record in cnc_rec:
        cellid = record[0]
        k, i, j = cellid
        ibd[k, i, j] = 3

    iper = 0
    ra = chd.stress_period_data.get_data(key=iper)
    ra

    # Read the binary head file and plot the results
    # We can use the existing Flopy HeadFile class because
    # the format of the headfile for MODFLOW 6 is the same
    # as for previous MODFLOW verions
    fname = os.path.join(workspace, headfile)
    hds = flopy.utils.binaryfile.HeadFile(fname)
    h = hds.get_data(kstpkper=(0, 0))
    x = y = np.linspace(0, L, N)
    y = y[::-1]
    c = plt.contour(x, y, h[0])
    plt.clabel(c, fmt="%2.1f")
    plt.title('Head file')
    plt.axis("scaled")
    #plt.show()

    # We can also use the Flopy PlotMapView capabilities for MODFLOW 6
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax)

    # Then we can use the plot_grid() method to draw the grid
    # The return value for this function is a matplotlib LineCollection object,
    # which could be manipulated (or used) later if necessary.
    quadmesh = modelmap.plot_ibound(ibound=ibd)
    linecollection = modelmap.plot_grid()
    contours = modelmap.contour_array(h[0], levels=np.arange(90, 100.1, 0.2))

    # We can also use the Flopy PlotMapView capabilities for MODFLOW 6
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")

    # Next we create an instance of the ModelMap class
    modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax)

    # Then we can use the plot_grid() method to draw the grid
    # The return value for this function is a matplotlib LineCollection object,
    # which could be manipulated (or used) later if necessary.
    quadmesh = modelmap.plot_ibound(ibound=ibd)
    linecollection = modelmap.plot_grid()
    pa = modelmap.plot_array(h[0])
    cb = plt.colorbar(pa, shrink=0.5)

    # plot well head layer 1
    # read binary file
    h = gwf.output.head().get_data(kstpkper=(0, 0))
    x = y = np.linspace(0, L, N)
    y = y[::-1]
    vmin, vmax = 90.0, 100.0
    contour_intervals = np.arange(90, 100.1, 1.0)

    # plot map layer 1
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    c = ax.contour(x, y, h[0], contour_intervals, colors="black")
    plt.clabel(c, fmt="%2.1f")

    # plot ap using modflow 6 capabilities
    fig, axes = plt.subplots(2, 1, figsize=(6, 12), constrained_layout=True)
    # first subplot
    ax = axes[0]
    ax.set_title("Model Layer 1")
    modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax)
    pa = modelmap.plot_array(h, vmin=vmin, vmax=vmax)
    quadmesh = modelmap.plot_bc("CHD")
    linecollection = modelmap.plot_grid(lw=0.5, color="0.5")
    contours = modelmap.contour_array(h, levels=contour_intervals, colors="black")
    ax.clabel(contours, fmt="%2.1f")
    cb = plt.colorbar(pa, shrink=0.5, ax=ax)
    # second subplot
    ax = axes[1]
    ax.set_title(f"Model Layer {Nlay}")
    modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=Nlay - 1)
    linecollection = modelmap.plot_grid(lw=0.5, color="0.5")
    pa = modelmap.plot_array(h, vmin=vmin, vmax=vmax)
    quadmesh = modelmap.plot_bc("CHD")
    contours = modelmap.contour_array(h, levels=contour_intervals, colors="black")
    ax.clabel(contours, fmt="%2.1f")
    cb = plt.colorbar(pa, shrink=0.5, ax=ax)
    #plt.show()

    # visualize plume contamination
    conc = gwt.output.concentration().get_data()[-1]

    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(1, 1, 1, aspect="equal")
    ax.set_title('Contamination Plume')
    pmv = flopy.plot.PlotMapView(gwf)
    c = pmv.plot_array(conc, cmap="jet")
    pmv.contour_array(conc, levels=(0.0001, 0.001, 0.01, 0.1), colors="y")
    plt.colorbar(c, shrink=0.5)
    plt.savefig(os.path.join(workspace, "plume.png"))
    #plt.show()

    # Save all plots

    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap
    import numpy as np

    # Create ibound array to identify boundary condition locations
    ibd = np.ones((Nlay, N, N), dtype=int)  # Start with all cells active (1)

    # Mark constant head cells with -1
    for record in chd_rec:
        cellid = record[0]
        k, i, j = cellid
        ibd[k, i, j] = -1

    # Mark contamination sources with 3
    for record in cnc_rec:
        cellid = record[0]
        k, i, j = cellid
        ibd[k, i, j] = 3

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
    #plt.show()

    # =======================================================================
    # 2. Plot Head Distribution
    # =======================================================================
    print("Plotting head distribution...")

    # Load head data
    head_file = os.path.join(workspace, f"{model_name}.hds")
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
    #plt.show()

    # Plots if the model has more layers

    if Nlay > 1:
        print("Plotting cross-section...")
        fig, ax = plt.subplots(figsize=(12, 6))

        # Create cross-section along row midpoint
        row = N // 2
        xsect = flopy.plot.PlotCrossSection(
            model=gwf,
            ax=ax,
            line={"row": row},
            geographic_coords=True
        )

        # Plot concentration in cross-section
        conc_cross = xsect.plot_array(conc, cmap='viridis', alpha=0.8, vmin=0, vmax=550)

        # Plot grid lines
        xsect.plot_grid(alpha=0.3)

        # Add colorbar
        plt.colorbar(conc_cross, shrink=0.5, label="Concentration")
        ax.set_title(f"Concentration Cross-Section (Row {row})")
        plt.savefig(os.path.join(workspace, "concentration_cross_section.png"))
        plt.show()

    print("All plots saved to:", workspace)


if __name__ == '__main__':
    example_2_wells()