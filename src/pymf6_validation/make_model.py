"""Create and run a MODFLOW 6 model with flopy.
"""

from pathlib import Path
import shutil
from itertools import zip_longest

import flopy

CONFIG = {
    'internal_dir': '.internal',
    'model_files': 'model_files',
}


def _get_mf6_exe(exe_name):
    """Get name of MODFLOW6 executable"""
    if exe_name is None:
        # make pymf6 optional to reduce dependencies
        import pymf6  # pylint: disable-msg=import-outside-toplevel
        exe_name = pymf6.__mf6_exe__
    return exe_name


def _save_model_file_names(
        model_path, model_file_names, append=False, config=CONFIG):
    """Save names of input model files to file"""
    internal = Path(model_path) / config['internal_dir']
    internal.mkdir(exist_ok=True)
    files = internal / config['model_files']
    old_content = ''
    if append:
        old_content = files.read_text() + '\n'
    files.write_text(old_content + '\n'.join(sorted(model_file_names)))


def make_input(
        model_data,
        exe_name=None,
        verbosity_level=0):
    """Create MODFLOW 6 input"""
    # pylint: disable-msg=too-many-locals
    exe_name = _get_mf6_exe(exe_name)
    model_path = model_data['model_path']
    model_name = model_data['name']
    file_extensions = ['nam']
    sim = flopy.mf6.MFSimulation(
        sim_name=model_data['name'],
        sim_ws=model_path,
        exe_name=exe_name,
        verbosity_level=verbosity_level,
    )
    times = model_data['times']
    repeat_times = model_data['repeat_times']
    tdis_rc = [(1.0, 1, 1.0)] + [times] * repeat_times
    pname = 'tdis'

    # Instantiating time discretization package
    flopy.mf6.ModflowTdis(
        sim, pname=pname,
        time_units=model_data['time_units'],
        nper=repeat_times + 1,
        perioddata=tdis_rc,
    )
    file_extensions.append(pname)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=model_name,
        save_flows=True)
    file_extensions.append('ims')
    dim_kwargs = {name: model_data[name] for name in
                  ['nrow', 'ncol', 'nlay', 'delr', 'delc', 'top', 'botm',
                   'length_units']
                  }
    #im_kwargs['lenuni'] = model_data['length_units']
    model_data['dim_kwargs'] = dim_kwargs

    # Instantiating spatial discretization package
    flopy.mf6.ModflowGwfdis(gwf, **dim_kwargs)
    file_extensions.append('dis')

    # Instantiating initial conditions package
    flopy.mf6.ModflowGwfic(gwf, strt=model_data['initial_head'])
    file_extensions.append('ic')
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=[0],
        k=model_data['k'],
        k33=model_data['k33'],
    )
    file_extensions.append('npf')
    sy = flopy.mf6.ModflowGwfsto.sy.empty(  # pylint: disable-msg=invalid-name
        gwf,
        default_value=model_data['sy']
    )
    ss = flopy.mf6.ModflowGwfsto.ss.empty(  # pylint: disable-msg=invalid-name
        gwf, default_value=model_data['ss']
    )
    pname = 'sto'

    # Instantiating storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        pname=pname,
        save_flows=True,
        iconvert=1,
        ss=ss,
        sy=sy,
        steady_state={0: True},
        transient={index: True for index in range(1, len(times))},
        )
    file_extensions.append(pname)

    if model_data['wells_active']:
        # Stress period data for the well
        stress_period_data = {}
        for index in range(len(times)):
            entry = []
            for well in model_data['wells'].values():
                value = [well['coords'], well['q'][index]]
                if model_data['transport']:
                    value.append(0)
                entry.append(tuple(value))
            stress_period_data[index + 1] = entry
        wel_kwargs = {}
        if model_data['transport']:
            wel_kwargs.update({
                'auxiliary': 'CONCENTRATION',
                'pname': 'WEL-1'})
        # Instantiating well package
        flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data=stress_period_data,
            **wel_kwargs,
        )
        file_extensions.append('wel')

    chd_kwargs = {}
    if model_data['transport']:
        chd_kwargs.update({
            'auxiliary': 'CONCENTRATION',
            'pname': 'CHD-1'})

    # Instantiating constant head package
    flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=model_data['chd'],
    **chd_kwargs
    )

    file_extensions.append('chd')

    budget_file = model_data['name'] + '.bud'
    head_file = model_data['name'] + '.hds'

    # Instantiating output control package
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budget_file,
        head_filerecord=head_file,
        saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])
    file_extensions.append('oc')

    if model_data['transport']:
        file_extensions.append('gwfgwt')
    model_file_names = set(f'{model_name}.{ext}' for ext in file_extensions)
    model_file_names.add('mfsim.nam')
    _save_model_file_names(model_path, model_file_names)
    if model_data['transport']:
        make_transport_model(sim, model_data)
    if model_data['river_active']:
        make_river(model_data=model_data, gwf=gwf, file_extensions=file_extensions)
    sim.write_simulation()


def make_transport_model(sim, model_data):
    """Create MODFLOW 6 input for transport model"""
    # Instantiating MODFLOW 6 groundwater transport package
    model_path = model_data['model_path']
    gwtname = 'gwt_' + model_data['name']
    gwt = flopy.mf6.MFModel(
        sim,
        model_type='gwt6',
        modelname=gwtname,
        model_nam_file=f'{gwtname}.nam'
    )
    file_extensions = ['nam']
    gwt.name_file.save_flows = True

    # create iterative model solution and register the gwt model with it
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
        filename=f'{gwtname}.ims'
    )
    file_extensions.append('ims')
    sim.register_ims_package(imsgwt, [gwt.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
        **model_data['dim_kwargs'],
        filename=f'{gwtname}.dis'
    )
    file_extensions.append('dis')

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGwtic(
        gwt, strt=model_data['initial_concentration'], filename=f'{gwtname}.ic'
    )
    file_extensions.append('ic')

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGwtadv(
        gwt, scheme=model_data['scheme'], filename=f'{gwtname}.adv'
    )
    file_extensions.append('adv')

    # Instantiating MODFLOW 6 transport dispersion package
    long_disp = model_data['longitudinal_dispersivity']
    ratio = model_data['dispersivity_ratio']
    if long_disp != 0:
        flopy.mf6.ModflowGwtdsp(
            gwt,
            xt3d_off=True,
            alh=long_disp,
            ath1=long_disp * ratio,
            filename=f'{gwtname}.dsp'
        )
    file_extensions.append('dsp')

    # Instantiating MODFLOW 6 transport mass storage package
    # (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=model_data['porosity'],
        first_order_decay=False,
        decay=None,
        decay_sorbed=None,
        sorption=None,
        bulk_density=None,
        distcoef=None,
        filename=f'{gwtname}.mst'
    )
    file_extensions.append('mst')

    # Instantiating MODFLOW 6 transport source-sink mixing package
    sourcerecarray = [
        ('CHD-1', 'AUX', 'CONCENTRATION'),
    ]

    if model_data['wells_active']:
        sourcerecarray.append(('WEL-1', 'AUX', 'CONCENTRATION'))

    flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f'{gwtname}.ssm'
    )
    file_extensions.append('ssm')
    if 'cnc' in model_data:
        flopy.mf6.ModflowGwtcnc(
            gwt,
            stress_period_data=model_data['cnc']
        )
        file_extensions.append('cnc')
    # Instantiating MODFLOW 6 transport output control package
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f'{gwtname}.cbc',
        concentration_filerecord=f'{gwtname}.ucn',
        concentrationprintrecord=[
            ('COLUMNS', 10, 'WIDTH', 15, 'DIGITS', 6, 'GENERAL')
        ],
        saverecord=[('CONCENTRATION', 'LAST'), ('BUDGET', 'LAST')],
        printrecord=[('CONCENTRATION', 'LAST'), ('BUDGET', 'LAST')],
    )
    file_extensions.append('oc')

    # Instantiate observation package (for transport)
    obs = model_data['obs']
    if obs:
        obslist = []
        for name, obs_coords in obs:
            obslist.append([name, 'concentration', obs_coords])
        obsdict = {f'{gwtname}.obs.csv': obslist}
        flopy.mf6.ModflowUtlobs(
            gwt, print_input=False, continuous=obsdict
        )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype='GWF6-GWT6',
        exgmnamea=model_data['name'],
        exgmnameb=gwtname,
        filename=f'{model_data["name"]}.gwfgwt'
    )
    model_file_names = set(f'{gwtname}.{ext}' for ext in file_extensions)
    _save_model_file_names(model_path, model_file_names, append=True)


def get_simulation(model_path, exe_name=None, verbosity_level=0):
    """Get simulation for a model."""
    exe_name = _get_mf6_exe(exe_name)
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=model_path,
        exe_name=exe_name,
        verbosity_level=verbosity_level,
    )
    return sim


def run_simulation(model_path, verbosity_level=0):
    """Run a MODFLOW 6 model"""
    sim = get_simulation(
        model_path,
        verbosity_level=verbosity_level)
    sim.run_simulation()


def clone_model(src, dst=None, config=CONFIG):
    """Clone a MODFLOW6 model.

    Copy all input files listed in `.internal/model_files from `src` to
    `dst`."""
    src = Path(src)
    if dst is None:
        dst = src.parent / (src.name + '_controlled')
    else:
        dst = Path(dst)
    dst.mkdir(exist_ok=True)
    model_files_path = src / config['internal_dir'] / config['model_files']
    model_files = model_files_path.read_text().split('\n')
    for file_name in model_files:
        shutil.copy(src / file_name, dst / file_name)

def make_river( model_data, file_extensions, gwf):
    entry = []
    riv = model_data['river_spd']
    for counter in range(len(riv['rivlay'])):
            entry.append([
                riv['rivlay'][counter],
                riv['rivrow'][counter],
                riv['rivcol'][counter],
                riv['rivstg'][counter],
                riv['rivcnd'][counter],
                riv['rivbot'][counter]])
    riv_kwargs = {}
    riv_kwargs.update({'pname':'RIV-1'})
    flopy.mf6.ModflowGwfriv(
            gwf,
            stress_period_data=entry,
            #boundnames=model_data['river_boundnames'], # boolean to indicate that boundary names may be in river list cells
            #observations= model_data['obs_dict'], # dictionary or data containing data for the observation package
           # timeseries= model_data['tsdict'], # dictionary or data for the time-series
            **riv_kwargs
            )
    file_extensions.append('riv')
    model_path=model_data['model_path']
    gwfname=model_data['name']
    model_file_names = set(f'{gwfname}.{ext}' for ext in file_extensions)
    _save_model_file_names(model_path, model_file_names, append=True)
