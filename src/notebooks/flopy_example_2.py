import os
import sys
from fileinput import filename
from pprint import pformat
from tempfile import TemporaryDirectory
from matplotlib.colors import LinearSegmentedColormap
from flopy.utils.binaryfile import HeadFile
from flopy.utils import UcnFile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import flopy

# For this example, we will set up a temporary workspace.
# Model input files and output files will reside here.
dir_name = os.getcwd()
workspace = os.path.join(dir_name + '/models', "pumptreat")

# Set up workspace
model_name = "pumptreat"

# create flopy objects
h1 = 25
h2 = 24.5
Nlay = 1
N = 101
L = 100.0
H = 30.0
k = 1.0
k33 = 0.3
q = -100.0
times = (100.0, 200, 1.0)

# Create the Flopy simulation object
sim = flopy.mf6.MFSimulation(
    sim_name=model_name, exe_name="C:/Users/lucialabarca/mf6.6.2_win64/bin/mf6.exe", version="mf6", sim_ws=workspace
)

# Create the Flopy temporal discretization object
tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
    sim, pname="tdis", time_units="DAYS", nper=2, perioddata=[(1.0, 1, 1.0), (100.0, 200, 1.0)]
)

# Create the Flopy groundwater flow (gwf) model object
model_nam_file = f"{model_name}.nam"
gwf = flopy.mf6.ModflowGwf(sim, modelname=model_name, model_nam_file=model_nam_file, save_flows=True)

# Solver for GWF (register this one FIRST)
ims_gwf = flopy.mf6.ModflowIms(
    sim,
    pname="ims_gwf",
    filename="gwf_" + model_name + ".ims",
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
    ((0, int(N / 2), int(N / 4)), q, 0),
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
headfile = f"{model_name}.hds"
budgetfile = f"{model_name}.bud"
saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]

oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
    gwf,
    pname="gwfgwt",
    budget_filerecord=budgetfile,
    head_filerecord=headfile,
    saverecord=saverecord,
    )

# create gwt name
gwtname = 'gwt_' + model_name
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
cnc_rec.append(((0, int(N/2), int(N/3)), 550.0))  # Single source
cnc_rec.append(((0, int(N/3), int(N/3)), 550.0))  # Second source

# Create CNC package
cnc = flopy.mf6.ModflowGwtcnc(
    gwt,
    stress_period_data=cnc_rec,
    pname="cnc"
)

# Output control for GWT
gwt_oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f'{gwtname}.cbc',
        concentration_filerecord=f'{gwtname}.ucn',
        concentrationprintrecord=[
            ('COLUMNS', 10, 'WIDTH', 15, 'DIGITS', 6, 'GENERAL')
        ],
        saverecord=[('CONCENTRATION', 'LAST'), ('BUDGET', 'LAST')],
        printrecord=[('CONCENTRATION', 'LAST'), ('BUDGET', 'LAST')],
    pname='oc'
    )

# GWF-GWT exchange
flopy.mf6.ModflowGwfgwt(sim,
                        exgtype="GWF6-GWT6",
                        exgmnamea=model_name,
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

# =======================================================================
# 1. Load Model Results
# =======================================================================
print("Loading model results...")

# Head file
head_file = os.path.join(dir_name + '/models', f"{model_name}.hds")
hds = HeadFile(head_file)
head = hds.get_data(kstpkper=(0, 0))  # First stress period

# Concentration file
ucn = flopy.utils.UcnFile(f"{gwtname}.ucn", model=gwt)
conc_obj = ucn.get_data()
conc = conc_obj.get_data(totim=100.0)  # Concentration at end of simulation

# Create grid coordinates
x = y = np.linspace(0, L, N)
y = y[::-1]  # Flip y-axis for proper orientation

# =======================================================================
# 2. Plot Boundary Conditions
# =======================================================================
print("Plotting boundary conditions...")

chd_rec = []
# chd_rec.append(((0, int(N), int(N)), h2))
for layer in range(0, Nlay):
    for row_col in range(0, N):
        chd_rec.append(((layer, row_col, 0), h1, 0.0))
        chd_rec.append(((layer, row_col, N - 1), h2, 0.0))

wel_rec = [
    ((0, int(N / 2), int(N / 4)), q, 0),
    ((0, int(N / 3), int(N / 4)), q, 0),
    ((0, int(N / 4), int(N / 4)), q, 0)
]

# CNC package
cnc_rec = []
# Add contamination sources
cnc_rec.append(((0, int(N / 2), int(N / 3)), 550.0))  # Single source
cnc_rec.append(((0, int(N / 3), int(N / 3)), 550.0))  # Second source

# Create CHD and WEL location arrays
chd_locs = np.zeros((Nlay, N, N), dtype=bool)
wel_locs = np.zeros((Nlay, N, N), dtype=bool)
cnc_locs = np.zeros((Nlay, N, N), dtype=bool)

# Mark CHD locations
for record in chd_rec:
    cellid = record[0]
    k, i, j = cellid
    chd_locs[k, i, j] = True

# Mark WEL locations
for record in wel_rec:
    cellid = record[0]
    k, i, j = cellid
    wel_locs[k, i, j] = True

# Mark CNC locations
for record in cnc_rec:
    cellid = record[0]
    k, i, j = cellid
    cnc_locs[k, i, j] = True

# Create boundary condition plot
fig, ax = plt.subplots(figsize=(10, 8))
modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)

# Plot grid
grid = modelmap.plot_grid(alpha=0.3)

# Plot boundary conditions
chd_plot = modelmap.plot_bc("CHD-1", color="blue", alpha=0.7, label="Constant Head")
wel_plot = modelmap.plot_bc("WEL-1", color="red", alpha=0.7, label="Well")
cnc_plot = ax.scatter(
    [x[j] for (_, _, j) in np.argwhere(cnc_locs[0])],
    [y[i] for (_, i, _) in np.argwhere(cnc_locs[0])],
    color="purple", s=50, marker="s", label="Contamination Source"
)

# Add legend and title
ax.legend(loc="upper right")
ax.set_title("Boundary Conditions")
plt.savefig(os.path.join(workspace, "boundary_conditions.png"))
plt.show()

# =======================================================================
# 3. Plot Head Distribution
# =======================================================================
print("Plotting head distribution...")

# Create head contour plot
fig, ax = plt.subplots(figsize=(10, 8))
modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)

# Plot grid
modelmap.plot_grid(alpha=0.3)

# Plot head contours
levels = np.linspace(head.min(), head.max(), 20)
contours = modelmap.contour_array(
    head[0],
    levels=levels,
    colors="black",
    linewidths=0.5
)
plt.clabel(contours, fmt="%.1f", fontsize=8)

# Plot filled head contours
head_plot = modelmap.plot_array(
    head[0],
    cmap="viridis",
    alpha=0.7
)
plt.colorbar(head_plot, shrink=0.5, label="Head (m)")

# Plot boundary conditions
modelmap.plot_bc("CHD-1", color="blue", alpha=0.5)
modelmap.plot_bc("WEL-1", color="red", alpha=0.7)

ax.set_title("Head Distribution")
plt.savefig(os.path.join(workspace, "head_distribution.png"))
plt.show()

# =======================================================================
# 4. Plot Concentration Distribution
# =======================================================================
print("Plotting concentration distribution...")

# Create custom colormap for concentrations
colors = [(0, "white"), (0.2, "cyan"), (0.5, "green"), (0.7, "yellow"), (1, "red")]
cmap = LinearSegmentedColormap.from_list("contam_cmap", colors)

fig, ax = plt.subplots(figsize=(10, 8))
modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)

# Plot concentration array
conc_plot = modelmap.plot_array(
    conc[0],
    cmap=cmap,
    alpha=0.8,
    vmin=0,
    vmax=550
)
plt.colorbar(conc_plot, shrink=0.5, label="Concentration")

# Plot concentration contours
conc_levels = [1, 10, 50, 100, 200, 300, 400, 500]
conc_contours = modelmap.contour_array(
    conc[0],
    levels=conc_levels,
    colors="black",
    linewidths=0.5
)
plt.clabel(conc_contours, fmt="%d", fontsize=8)

# Plot contamination sources
ax.scatter(
    [x[j] for (_, _, j) in np.argwhere(cnc_locs[0])],
    [y[i] for (_, i, _) in np.argwhere(cnc_locs[0])],
    color="purple", s=50, marker="s", label="Source"
)

ax.set_title("Concentration Distribution at End of Simulation")
plt.savefig(os.path.join(workspace, "concentration_distribution.png"))
plt.show()

# =======================================================================
# 5. Plot Cross-Section (If multi-layer)
# =======================================================================
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
    conc_cross = xsect.plot_array(
        conc,
        cmap=cmap,
        alpha=0.8,
        vmin=0,
        vmax=550
    )

    # Plot grid lines
    xsect.plot_grid(alpha=0.3)

    # Add colorbar
    plt.colorbar(conc_cross, shrink=0.5, label="Concentration")
    ax.set_title(f"Concentration Cross-Section (Row {row})")
    plt.savefig(os.path.join(workspace, "concentration_cross_section.png"))
    plt.show()

# =======================================================================
# 6. Create Animation of Contaminant Plume Development
# =======================================================================
print("Creating concentration animation...")

# Get all concentration times
times = conc_obj.get_times()
nsteps = min(10, len(times))  # Limit to 10 frames for animation
selected_times = np.linspace(0, len(times) - 1, nsteps, dtype=int)

# Create animation frames
animation_dir = os.path.join(workspace, "animation")
os.makedirs(animation_dir, exist_ok=True)

for idx, time_idx in enumerate(selected_times):
    t = times[time_idx]
    conc_frame = conc_obj.get_data(totim=t)

    fig, ax = plt.subplots(figsize=(10, 8))
    modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)

    # Plot concentration
    conc_plot = modelmap.plot_array(
        conc_frame[0],
        cmap=cmap,
        alpha=0.8,
        vmin=0,
        vmax=550
    )

    # Plot concentration contours
    conc_contours = modelmap.contour_array(
        conc_frame[0],
        levels=conc_levels,
        colors="black",
        linewidths=0.5
    )

    # Plot sources
    ax.scatter(
        [x[j] for (_, _, j) in np.argwhere(cnc_locs[0])],
        [y[i] for (_, i, _) in np.argwhere(cnc_locs[0])],
        color="purple", s=50, marker="s"
    )

    plt.colorbar(conc_plot, shrink=0.5, label="Concentration")
    ax.set_title(f"Concentration at Time = {t:.1f} days")
    plt.savefig(os.path.join(animation_dir, f"conc_frame_{idx:03d}.png"))
    plt.close()

print("Visualization complete! Check the output directory for results.")
print(f"All plots saved to: {workspace}")
print(f"Animation frames saved to: {animation_dir}")































# Read the binary head file and plot the results
# We can use the existing Flopy HeadFile class because
# the format of the headfile for MODFLOW 6 is the same
# as for previous MODFLOW verions
fname = os.path.join(workspace, headfile)
hds = flopy.utils.binaryfile.HeadFile(fname)
h = hds.get_data(kstpkper=(0, 0))
x = y = np.linspace(0, L, N)
y = y[::-1]
c = plt.contour(x, y, h[0], np.arange(90, 100.1, 0.2))
plt.clabel(c, fmt="%2.1f")
plt.axis("scaled")
plt.show()


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
plt.show()

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
plt.show()

# visualize plume contamination
conc = gwt.output.concentration().get_data()

fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(1, 1, 1, aspect="equal")
pmv = flopy.plot.PlotMapView(gwf)
c = pmv.plot_array(conc, cmap="jet")
pmv.contour_array(conc, levels=(0.0001, 0.001, 0.01, 0.1), colors="y")
plt.colorbar(c, shrink=0.5)