import numpy as np
import pickle

import openmc


def build_inputs(n_mesh_bins, energy_bins, n_azim_bins, **kwargs):
    """Build OpenMC input XML files for a pincell with detailed flux tallying"""
    if n_azim_bins % 8 != 0: raise ValueError("The number of azimuthal bins "
         "must be divisible by eight")
    pitch = 0.62992*2

    ####################
    # Define materials
    ####################

    uo2 = openmc.Material(material_id=9201,
                          name='UO2 fuel at 2.4% wt enrichment')
    uo2.set_density('g/cm3', 10.29769)
    uo2.add_nuclide('U-234', 4.4842e-06)
    uo2.add_nuclide('U-235', 5.5814e-04)
    uo2.add_nuclide('U-238', 2.2407e-02)
    uo2.add_nuclide('O-16', 4.5828e-02)
    uo2.add_nuclide('O-17', 1.7457e-05 + 9.4176e-05)

    borated_water = openmc.Material(material_id=101,
                                    name='Borated water at 975 ppm')
    borated_water.set_density('g/cm3', 0.740582)
    borated_water.add_nuclide('B-10', 8.0042e-6)
    borated_water.add_nuclide('B-11', 3.2218e-5)
    borated_water.add_nuclide('H-1', 4.9457e-2)
    borated_water.add_nuclide('H-2', 7.4196e-6)
    borated_water.add_nuclide('O-16', 2.4672e-2)
    borated_water.add_nuclide('O-17', 9.3982e-06 + 5.0701e-05)
    borated_water.add_s_alpha_beta('HH2O', '71t')

    helium = openmc.Material(material_id=201, name='Helium for gap')
    helium.set_density('g/cm3', 0.001598)
    helium.add_nuclide('He-4', 2.4044e-4)

    zircaloy = openmc.Material(material_id=4001, name='Zircaloy 4')
    zircaloy.set_density('g/cm3', 6.55)
    zircaloy.add_nuclide('O-16', 3.0743e-04)
    zircaloy.add_nuclide('O-17', 1.1711e-07 + 6.3176e-07)
    zircaloy.add_nuclide('Cr-50', 3.2962e-06)
    zircaloy.add_nuclide('Cr-52', 6.3564e-05)
    zircaloy.add_nuclide('Cr-53', 7.2076e-06)
    zircaloy.add_nuclide('Cr-54', 1.7941e-06)
    zircaloy.add_nuclide('Fe-54', 8.6699e-06)
    zircaloy.add_nuclide('Fe-56', 1.3610e-04)
    zircaloy.add_nuclide('Fe-57', 3.1431e-06)
    zircaloy.add_nuclide('Fe-58', 4.1829e-07)
    zircaloy.add_nuclide('Zr-90', 2.1827e-02)
    zircaloy.add_nuclide('Zr-91', 4.7600e-03)
    zircaloy.add_nuclide('Zr-92', 7.2758e-03)
    zircaloy.add_nuclide('Zr-94', 7.3734e-03)
    zircaloy.add_nuclide('Zr-96', 1.1879e-03)
    zircaloy.add_nuclide('Sn-112', 4.6735e-06)
    zircaloy.add_nuclide('Sn-114', 3.1799e-06)
    zircaloy.add_nuclide('Sn-115', 1.6381e-06)
    zircaloy.add_nuclide('Sn-116', 7.0055e-05)
    zircaloy.add_nuclide('Sn-117', 3.7003e-05)
    zircaloy.add_nuclide('Sn-118', 1.1669e-04)
    zircaloy.add_nuclide('Sn-119', 4.1387e-05)
    zircaloy.add_nuclide('Sn-120', 1.5697e-04)
    zircaloy.add_nuclide('Sn-122', 2.2308e-05)
    zircaloy.add_nuclide('Sn-124', 2.7897e-05)

    materials_file = openmc.MaterialsFile()
    materials_file.default_xs = '71c'
    materials_file.add_materials([uo2, helium, zircaloy, borated_water])
    materials_file.export_to_xml()

    ####################
    # Define geometry
    ####################

    # Surfaces.
    fuel_or = openmc.ZCylinder(R=0.39218, name='Fuel OR')
    clad_ir = openmc.ZCylinder(R=0.40005, name='Clad IR')
    clad_or = openmc.ZCylinder(R=0.45720, name='Clad OR')

    bottom = openmc.YPlane(y0=0.0, name='bottom')
    right = openmc.XPlane(x0=pitch/2.0, name='right')
    top = openmc.Plane(A=1, B=-1, name='top') # 45 degree angle
    lower = openmc.ZPlane(z0=-10, name='lower')
    upper = openmc.ZPlane(z0=10, name='upper')
    for s in (bottom, right, top, lower, upper): s.boundary_type = 'reflective'

    # Cells.
    fuel_cell = openmc.Cell()
    gap_cell = openmc.Cell()
    clad_cell = openmc.Cell()
    water_cell = openmc.Cell()

    # Cell regions.
    fuel_cell.region = -fuel_or
    gap_cell.region = +fuel_or & -clad_ir
    clad_cell.region = +clad_ir & -clad_or
    water_cell.region = +clad_or & -right
    for c in (fuel_cell, gap_cell, clad_cell, water_cell):
        c.region = c.region & +bottom & +top & +lower & -upper

    # Cell fills.
    fuel_cell.fill = uo2
    gap_cell.fill = helium
    clad_cell.fill = zircaloy
    water_cell.fill = borated_water

    # Universe, geometry, and XML.
    root = openmc.Universe(universe_id=0, name='Root universe')
    root.add_cells([fuel_cell, gap_cell, clad_cell, water_cell])

    geometry = openmc.Geometry()
    geometry.root_universe = root

    geometry_file = openmc.GeometryFile()
    geometry_file.geometry = geometry
    geometry_file.export_to_xml()

    ####################
    # Define settings
    ####################

    settings_file = openmc.SettingsFile()
    settings_file.batches = kwargs.setdefault('batches', 25)
    settings_file.inactive = kwargs.setdefault('inactive', 5)
    settings_file.particles = kwargs.setdefault('particles', 1000)
    settings_file.source = openmc.source.Source(
                                            openmc.stats.Point((0.2, 0.1, 0.0)))
    settings_file.export_to_xml()

    ####################
    # Define tallies
    ####################

    mesh = openmc.Mesh(mesh_id=1)
    delta = pitch / 2.0 / (n_mesh_bins + 1)
    mesh.dimension = [n_mesh_bins, n_mesh_bins, 1]
    mesh.lower_left = [-delta/2.0, -delta/2.0, -1.e50]
    mesh.upper_right = [pitch/2.0 + delta/2.0, pitch/2.0 + delta/2.0, 1.e50]

    energy_filter = openmc.Filter(type='energy', bins=energy_bins)
    mesh_filter = openmc.Filter()
    mesh_filter.mesh = mesh
    azim_filter = openmc.Filter(type='azimuthal', bins=n_azim_bins)

    tally = openmc.Tally(tally_id=1)
    tally.add_filter(energy_filter)
    tally.add_filter(mesh_filter)
    tally.add_filter(azim_filter)
    tally.add_score('flux')

    tallies_file = openmc.TalliesFile()
    tallies_file.add_mesh(mesh)
    tallies_file.add_tally(tally)
    tallies_file.export_to_xml()

    ####################
    # Define plots
    ####################

    plotfile = openmc.PlotsFile()

    plot = openmc.Plot()
    plot.filename = 'matplot'
    plot.origin = (pitch/4.0, pitch/4.0, 0)
    plot.width = (0.5*pitch, 0.5*pitch)
    plot.pixels = (400, 400)
    plot.color = 'mat'
    plot.col_spec = {
         101: (100, 200, 200),
         201: (220, 220, 220),
         4001: (150, 150, 150),
         9201: (255, 50, 50)}
    plotfile.add_plot(plot)

    plot = openmc.Plot()
    plot.filename = 'cellplot'
    plot.origin = (pitch/4.0, pitch/4.0, 0)
    plot.width = (0.5*pitch, 0.5*pitch)
    plot.pixels = (400, 400)
    plot.color = 'cell'
    plotfile.add_plot(plot)

    plotfile.export_to_xml()


def extract_data(statepoint_name, output_name):
    """Save the tally data from a statepoint as a pickled DataFrame."""
    sp = openmc.StatePoint(statepoint_name)
    with open(output_name, 'wb') as fh:
        pickle.dump(sp.tallies[1].get_pandas_dataframe(), fh)


def small_example():
    n_mesh_bins = 10
    energy_bins = [0.0, 0.625, 20.0]
    n_azim_bins = 16
    build_inputs(n_mesh_bins, energy_bins, n_azim_bins, 
         inactive=5, batches=25, particles=1000)


def medium_example():
    n_mesh_bins = 10
    energy_bins = np.logspace(-8, 1, 9+1)
    n_azim_bins = 16
    build_inputs(n_mesh_bins, energy_bins, n_azim_bins, 
         inactive=5, batches=55, particles=10000)


def big_example():
    n_mesh_bins = 10
    energy_bins = np.logspace(-8, 1, 9+1)
    n_azim_bins = 200
    build_inputs(n_mesh_bins, energy_bins, n_azim_bins, 
         inactive=5, batches=105, particles=int(1e6))


if __name__ == '__main__':
    small_example()
    #extract_data('statepoint.25.h5', 'temp.p')
