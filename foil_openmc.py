import openmc
import matplotlib.pyplot as plt
import numpy as np
import openmc.data.data as openmcdata
import pandas

def foil_simulation(df_gamma, distance_foil, energy_bins):
    #Material
    #lead
    pb = openmc.Material(name='lead')
    pb.add_element('Pb', 1, 'ao')
    pb.set_density('g/cm3', 11.34)
    #copper
    #cu = openmc.Material(name='copper')
    #cu.add_element('Cu', 1, 'ao')
    #cu.set_density('g/cm3', 8.96)
    #aluminium
    al = openmc.Material(name='aluminium')
    al.add_element('Al', 1, 'ao')
    al.set_density('g/cm3', 2.7)
    #al.add_s_alpha_beta("c_Al27")
    #beryllium
    be = openmc.Material(name='beryllium')
    be.add_element('Be', 1, 'ao')
    be.set_density('g/cm3', 1.848)
    #germanium
    ge = openmc.Material(name="germanium")
    ge.add_element('Ge', 1, 'ao')
    ge.set_density('g/cm3', 5.323)
    #polythene terephelate (mylar)
    mylar = openmc.Material(name="mylar")
    mylar.add_element('H', 0.041959, 'wo')
    mylar.add_nuclide('C0', 0.625017, 'wo')
    mylar.add_element('O', 0.333025, 'wo')
    mylar.set_density('g/cm3', 1.4)

    materials = openmc.Materials([al, be, ge, mylar, pb])
    materials.cross_sections = "nuclearcx_lib/cross_sections.xml" #endfb-vii.1
    materials.export_to_xml(f"materials.xml")

    #Geometry (cm)
    detector_radius = 2.465
    detector_hollow_radius = 0.595
    detector_hollow_height = 4.08
    detector_hollow_nominal_radius = 0.5
    detector_height = 5.29
    al_case_o_cyl_o_radius = 3.5
    al_case_o_cyl_i_radius = 3.4
    al_case_o_cyl_height = 9.5
    al_case_o_cyl_pos = 0.3
    al_case_o_cap_cyl_radius = 3.5
    al_case_o_cap_cyl_height = 0.1
    al_case_o_cap_cyl_pos = 0.3
    al_case_i_cyl_o_radius = 2.88 #uncertain variablele
    al_case_i_cyl_i_radius = 2.8 #uncertain variable
    al_case_i_cyl_height = 9.4
    al_end_cyl_radius = al_case_i_cyl_o_radius
    al_end_cyl_height = 0.3
    al_end_cyl_pos = -9.1
    almy_cap_cyl_radius = al_case_i_cyl_o_radius
    almy_cap_cyl_height = 0.03
    almy_cap_cyl_pos = 0
    almy_cap_2_cyl_radius = al_case_i_cyl_o_radius
    almy_cap_2_cyl_height = 0.03
    almy_cap_2_cyl_pos = almy_cap_cyl_pos + almy_cap_cyl_height
    be_cap_cyl_radius = al_case_i_cyl_o_radius
    be_cap_cyl_height = 0.05
    be_cap_cyl_pos = 0.3
    leadshield_distance = 5
    leadshield_position_v = 1
    leadshield_position_h = -5
    leadshield_length = 10
    leadshield_width = 1
    leadshield_height = 9
    sourcedistance = distance_foil + be_cap_cyl_pos + be_cap_cyl_height

    world_boundary = openmc.Sphere(r=21, boundary_type='vacuum')
    detector_cyl = openmc.ZCylinder(r=detector_radius, boundary_type='transmission')
    detector_bot = openmc.ZPlane(z0=-detector_height, boundary_type='transmission')
    detector_top = openmc.ZPlane(z0=0, boundary_type='transmission')
    detector_hollow_cyl = openmc.ZCylinder(r=detector_hollow_radius, boundary_type='transmission')
    detector_hollow_top = openmc.ZPlane(z0=-detector_height+detector_hollow_height, boundary_type='transmission')
    detector_hollow_bot = openmc.ZPlane(z0=-detector_height, boundary_type='transmission')
    detector_hollow_nominal_sphere = openmc.Sphere(r=0.5,  x0=0, y0=0, z0=-detector_height+detector_hollow_height, boundary_type='transmission')
    al_case_o_cyl_o = openmc.ZCylinder(r=al_case_o_cyl_o_radius, boundary_type='transmission')
    al_case_o_cyl_i = openmc.ZCylinder(r=al_case_o_cyl_i_radius, boundary_type='transmission')
    al_case_o_bot = openmc.ZPlane(z0=al_case_o_cyl_pos-al_case_o_cyl_height, boundary_type='transmission')
    al_case_o_top = openmc.ZPlane(z0=al_case_o_cyl_pos, boundary_type='transmission')
    al_case_o_cap_cyl = openmc.ZCylinder(r=al_case_o_cap_cyl_radius, boundary_type='transmission')
    al_case_o_cap_bot = openmc.ZPlane(z0=al_case_o_cap_cyl_pos, boundary_type='transmission')
    al_case_o_cap_top = openmc.ZPlane(z0=al_case_o_cap_cyl_pos+al_case_o_cap_cyl_height, boundary_type='transmission')
    al_case_i_cyl_o = openmc.ZCylinder(r=al_case_i_cyl_o_radius, boundary_type='transmission')
    al_case_i_cyl_i = openmc.ZCylinder(r=al_case_i_cyl_i_radius, boundary_type='transmission')
    al_case_i_bot = openmc.ZPlane(z0=-al_case_i_cyl_height, boundary_type='transmission')
    al_case_i_top = openmc.ZPlane(z0=0, boundary_type='transmission')
    al_end_cyl = openmc.ZCylinder(r=al_end_cyl_radius, boundary_type='transmission')
    al_end_bot = openmc.ZPlane(z0=al_end_cyl_pos-al_end_cyl_height, boundary_type='transmission')
    al_end_top = openmc.ZPlane(z0=al_end_cyl_pos, boundary_type='transmission')
    almy_cap_cyl = openmc.ZCylinder(r=almy_cap_cyl_radius, boundary_type='transmission')
    almy_cap_bot = openmc.ZPlane(z0=almy_cap_cyl_pos, boundary_type='transmission')
    almy_cap_top = openmc.ZPlane(z0=almy_cap_cyl_pos+almy_cap_cyl_height, boundary_type='transmission')
    almy_cap_2_cyl = openmc.ZCylinder(r=almy_cap_2_cyl_radius, boundary_type='transmission')
    almy_cap_2_bot = openmc.ZPlane(z0=almy_cap_2_cyl_pos, boundary_type='transmission')
    almy_cap_2_top = openmc.ZPlane(z0=almy_cap_2_cyl_pos+almy_cap_2_cyl_height, boundary_type='transmission')
    be_cap_cyl = openmc.ZCylinder(r=be_cap_cyl_radius, boundary_type='transmission')
    be_cap_bot = openmc.ZPlane(z0=be_cap_cyl_pos, boundary_type='transmission')
    be_cap_top = openmc.ZPlane(z0=be_cap_cyl_pos+be_cap_cyl_height, boundary_type='transmission')
    lead_shield_top = openmc.ZPlane(z0=leadshield_position_v, boundary_type='transmission')
    lead_shield_bottom = openmc.ZPlane(z0=leadshield_position_v-leadshield_height, boundary_type='transmission')
    lead_shield_1_side_inner = openmc.XPlane(x0=leadshield_distance, boundary_type='transmission')
    lead_shield_1_side_outer = openmc.XPlane(x0=leadshield_distance+leadshield_width, boundary_type='transmission')
    lead_shield_2_side_inner = openmc.XPlane(x0=-leadshield_distance, boundary_type='transmission')
    lead_shield_2_side_outer = openmc.XPlane(x0=-leadshield_distance-leadshield_width, boundary_type='transmission')
    lead_shield_side_b = openmc.YPlane(y0=leadshield_position_h+leadshield_length, boundary_type='transmission')
    lead_shield_side_a = openmc.YPlane(y0=leadshield_position_h, boundary_type='transmission')

    detector_hollow_region = -detector_hollow_top & +detector_hollow_bot & -detector_hollow_cyl
    detector_hollow_nominal_sphere_region = -detector_hollow_nominal_sphere
    detector_region = -detector_top & +detector_bot & -detector_cyl & ~detector_hollow_region & ~detector_hollow_nominal_sphere_region
    al_case_o_region = -al_case_o_top & +al_case_o_bot & -al_case_o_cyl_o & +al_case_o_cyl_i | -al_case_o_cap_top & +al_case_o_cap_bot & -al_case_o_cap_cyl & +be_cap_cyl
    al_case_i_region = -al_case_i_top & +al_case_i_bot & -al_case_i_cyl_o & +al_case_i_cyl_i | -al_end_top & +al_end_bot & -al_end_cyl & +detector_hollow_cyl
    almy_cap_region = -almy_cap_top & +almy_cap_bot & -almy_cap_cyl
    almy_cap_2_region = -almy_cap_2_top & +almy_cap_2_bot & -almy_cap_2_cyl
    be_cap_region = -be_cap_top & +be_cap_bot & -be_cap_cyl
    lead_shield_1 = -lead_shield_top & +lead_shield_bottom & +lead_shield_side_a & -lead_shield_side_b & +lead_shield_1_side_inner & -lead_shield_1_side_outer
    lead_shield_2 = -lead_shield_top & +lead_shield_bottom & +lead_shield_side_a & -lead_shield_side_b & -lead_shield_2_side_inner & +lead_shield_2_side_outer
    vacuum_void_region = -world_boundary & ~detector_region & ~al_case_i_region & ~al_case_o_region & ~almy_cap_region & ~almy_cap_2_region & ~be_cap_region & ~lead_shield_1 & ~lead_shield_2

    detector_cell = openmc.Cell(region=detector_region, fill=ge)
    al_case_o_cell = openmc.Cell(region=al_case_o_region, fill=al)
    al_case_i_cell = openmc.Cell(region=al_case_i_region, fill=al)
    almy_cap_cell = openmc.Cell(region=almy_cap_region, fill=al)
    almy_cap_2_cell = openmc.Cell(region=almy_cap_2_region, fill=mylar)
    be_cap_cell = openmc.Cell(region=be_cap_region, fill=be)
    lead_shield_cell_1 = openmc.Cell(region=lead_shield_1, fill=pb)
    lead_shield_cell_2 = openmc.Cell(region=lead_shield_2, fill=pb)
    vacuum_void_cell = openmc.Cell(region=vacuum_void_region)

    #universe
    universe = openmc.universe = openmc.Universe(cells=[detector_cell, al_case_o_cell, al_case_i_cell, almy_cap_cell, almy_cap_2_cell, be_cap_cell, vacuum_void_cell, lead_shield_cell_1, lead_shield_cell_2])
    geometry = openmc.Geometry(universe)
    geometry.export_to_xml(f"geometry.xml")

    #plot
    color_assignment = {detector_cell : 'steelblue',
                        al_case_o_cell : 'tan',
                        al_case_i_cell : 'tan',
                        be_cap_cell : 'indigo',
                        almy_cap_cell : 'gold',
                        almy_cap_2_cell : 'orange',
                        lead_shield_cell_1 : 'dimgrey',
                        lead_shield_cell_2 : 'dimgrey',
                        vacuum_void_cell : 'whitesmoke',
                        }

    #plotxz = geometry.plot(pixels=(1000,1000), basis='xz', color_by='cell' ,colors=color_assignment)
    #plotxz.figure.savefig("xz-cell.png")
    #plotxy = geometry.plot(pixels=(1000,1000), origin=(0, 0, be_cap_cyl_pos + be_cap_cyl_height), basis='xy', color_by='cell', colors=color_assignment)
    #plotxy.figure.savefig("xy-cell.png")
    #plotyz = geometry.plot(pixels=(1000,1000), basis='yz', colors=color_assignment)
    #plotyz.figure.savefig('yz-cell.png')

    #Source
    source = openmc.IndependentSource()
    source.particle = "photon"
    source.space = openmc.stats.Point((0,0, sourcedistance)) #source location
    source.energy = openmc.stats.Discrete(df_gamma.energy.to_numpy() * 1000, df_gamma.intensity_real.to_numpy()) #source energy profile
    source.angle = openmc.stats.Isotropic()
    source.strength = df_gamma.decay.sum()
    
    #Settings
    sim_batch = 20
    sim_particle = 200000

    settings = openmc.Settings()
    settings.particles = sim_particle 
    settings.batches = sim_batch
    settings.photon_transport = True
    settings.run_mode = "fixed source"
    settings.source = source
    settings.verbosity = 1
    settings.material_cell_offsets = False
    settings.threads = 1
    settings.export_to_xml(f"settings.xml")

    #getting eneregy bins from experimental data (convert kev to ev)
    #first_energy_bin = energy_data[0] * 1000
    #last_energy_bin = energy_data[bin_no - 1] * 1000
    #Tallies
    tallies = openmc.Tallies()
    photon_particle_filter = openmc.ParticleFilter(['photon'])
    energy_filter = openmc.EnergyFilter(energy_bins * 1000)
    detector_cell_filter = openmc.CellFilter(detector_cell)

    #tally_1 (detector pulse-height)
    tally_1 = "pulse-height"
    tally_1_type = "pulse-height"

    tallies_tally_1 = openmc.Tally(name=tally_1)
    tallies_tally_1.scores = [tally_1_type]
    tallies_tally_1.filters = [detector_cell_filter, energy_filter]
    tallies.append(tallies_tally_1)
    tallies.export_to_xml(f"tallies.xml")

    openmc.run()