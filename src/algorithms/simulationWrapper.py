from .flowSolver.simulator import structure_solver
from .mesher.structureMesh_class import mark_inside ,mark_boundaries
from .mesher.close_holes import close_holes
from .plot.plot_class import plot_fenics_object, plot_subdomain

import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np

from skimage.io import imread
from dolfin import UnitSquareMesh, File
from pathlib import Path
from zipfile import ZipFile


# l_s = ['/home/rafael/Documents/web-app/python-code/Fenics-Simulation/simu_wrapper.py', '/home/rafael/Desktop/imagetest.png', 'coarseMesh', 'black', 'Stokes', '10', '-10', '0.001']
# subprocess.Popen(l_s)
# python3 simu_wrapper rock_img_path (nx,ny) free_flow_pixel flow_equation boundary_condition (pin,pout) mu k
# python3 simu_wrapper.py '/home/rafael/Desktop/web_dev/Fenics-Simulation/test-images/image2.png' [100,100] 0 'Brinkman' 'Slip' [1,-1] 1e-3 1e-11


def simu_wrapper(*args, k = None, del_cells = None):

    ROOT_PATH = os.path.dirname(os.path.abspath(__file__))
    ROOT_PATH = os.path.dirname(ROOT_PATH)

    for i_input, val in enumerate(args):
        print('Input {} has value {}'.format(i_input, val))

    rock_image_path = args[0]  
    mesh_refinement = args[1]
    free_flow_pixel = args[2]
    flow_equation = args[3]
    pin = eval(args[4])  # 1, -1 #input 6
    pout = eval(args[5])  # 1, -1 #input 6
    mu = eval(args[6])  # 1e-3 #input 7

    if free_flow_pixel == "black":
        no_flow_pixel_n = 255
        free_flow_pixel_n = 0  # 0 #input 3
    elif free_flow_pixel == "white":
        no_flow_pixel_n = 0
        free_flow_pixel_n = 255  # 0 #input 3

    if flow_equation == "Brinkman":
        k = eval(k)  # 1e-10 #input 8


    #####################################################################
    ###################### -- GENERAL STUFF -- ##########################
    #####################################################################
    # rock-image and the pixels array of the image
    results_path = rock_image_path[:rock_image_path.rfind('/')]

    # "/home/rafael/Desktop/web_dev/Fenics-Simulation/test-images/aspect-test.png" #input 1
    rock_array = imread(rock_image_path, as_gray=True)
    if len(rock_array.shape) == 3:
        rock_array = rock_array[:,:,0]

    if flow_equation == 'Stokes':
        rock_array = close_holes(rock_array, free_flow_pixel_n, no_flow_pixel_n)

    #####################################################################
    ######################### -- MESH STUFF -- ##########################
    #####################################################################
    # number of elements for each direction and structured mesh
    if mesh_refinement == "coarseMesh":
        nx_elements, ny_elements = 20, 20  # 30, 30 # #input 2
    elif mesh_refinement == "fineMesh":
        nx_elements, ny_elements = 50, 50  # 30, 30 # #input 2
    elif mesh_refinement == "evenFinerMesh":
        nx_elements, ny_elements = 80, 80  # 30, 30 # #input 2
    else:
        nx_elements, ny_elements = 3, 3  # 30, 30 # #input 2

    square_mesh = UnitSquareMesh(nx_elements, ny_elements, diagonal="crossed")

    # n_iterations set the number of pixel points to be analysed inside of each triangle
    # n_iteration = 1 -> center point + each vertex
    # n_iteartion = 2 -> same of n_iteration = 1 + the mean of each combinations
    # marker -> all elements are marked with 1, mark_inside marks the other subdomains with marker int
    n_iteration = 1


    square_mesh, subdomains = mark_inside(
        square_mesh, rock_array, n_iteration, free_flow_pixel_n, flow_equation, del_cells
    )

    boundaries = mark_boundaries(square_mesh, flow_equation)

    # mesh and subdomains for simulation
    mesh_folder_name = "meshes"
    mesh_dir_path = "/".join([results_path, mesh_folder_name])
    if not os.path.isdir(mesh_dir_path):
        os.mkdir(mesh_dir_path)
    
    mesh_xml_path = "/".join([mesh_dir_path, "mesh.xml"])
    boundaries_xml_path = "/".join([mesh_dir_path, "boundaries.xml"])
    subdomain_xml_path = "/".join([mesh_dir_path, "subdomain.xml"])

    File(mesh_xml_path) << square_mesh
    File(boundaries_xml_path) << boundaries
    File(subdomain_xml_path) << subdomains
    #####################################################################
    ################### -- SIMULATION STUFF -- ##########################
    #####################################################################
    # simulation inputs
    u, p = structure_solver(
        mesh_xml_path, boundaries_xml_path, subdomain_xml_path, flow_equation, pin, pout, mu, k=k
    )

    #####################################################################
    ###################### -- EXTRA SAVINGS -- ##########################
    #####################################################################
    # mesh - image configurations
    edgecolor = "blue"  # or 'none'
    linewidth = 0.05  # or 'none'
    fig = plot_subdomain(subdomains, edgecolor, linewidth)

    # save figure
    img_folder_name = "images"
    img_dir_path = "/".join([results_path, img_folder_name])
    if not os.path.isdir(img_dir_path):
        os.mkdir(img_dir_path)

    fig_path = "/".join([img_dir_path, "mesh.png"])
    plt.savefig(fig_path)

    # mesh and subdomains for paraviewimg_dir_path = "/".join([ROOT_PATH, results_path, "images"])
    mesh_vizu_folder_name = "mesh_visualization"
    paraview_mesh_dir_path = "/".join([results_path, mesh_vizu_folder_name])
    if not os.path.isdir(paraview_mesh_dir_path):
        os.mkdir(paraview_mesh_dir_path)

    # mesh_paraview_path = "/".join([paraview_mesh_dir_path, "mesh.pvd"])
    subdomain_paraview_path = "/".join([paraview_mesh_dir_path, "subdomain.pvd"])

    # File(mesh_paraview_path) << square_mesh
    File(subdomain_paraview_path) << subdomains

    # save csv of results
    results_folder_name = "results_csv"
    csv_dir_path = "/".join([results_path, results_folder_name])
    if not os.path.isdir(csv_dir_path):
        os.mkdir(csv_dir_path)

    mesh_coords = u.function_space().mesh().coordinates()

    u_vertex = u.compute_vertex_values()
    u_x_vertex = u_vertex[: int(len(u_vertex) / 2)]
    u_y_vertex = u_vertex[int(len(u_vertex) / 2) :]

    p_vertex = p.compute_vertex_values()

    header = "x,y,velocity_x(pixel/s),velocity_y(pixel/s),pressure(Pa)"
    csv_path = "/".join([csv_dir_path, "vertex_pressure_velocity.csv"])
    np.savetxt(
        csv_path,
        np.c_[mesh_coords[:, 0], mesh_coords[:, 1], u_x_vertex, u_y_vertex, p_vertex],
        header=header,
        delimiter=",",
    )

    # save pvd of results
    results_vizu_folder_name = "results_visualization"
    paraview_dir_path = "/".join([results_path, results_vizu_folder_name])
    if not os.path.isdir(paraview_dir_path):
        os.mkdir(paraview_dir_path)

    velocity_pvd_path = "/".join([paraview_dir_path, "velocity.pvd"])
    pressure_pvd_path = "/".join([paraview_dir_path, "pressure.pvd"])

    File(velocity_pvd_path) << u
    File(pressure_pvd_path) << p

    # save velocity and pressure field figures
    velocity_field_fig_path = "/".join([img_dir_path, "velocity_field.png"])
    fig_u = plot_fenics_object(u)
    plt.savefig(velocity_field_fig_path, bbox_inches="tight", pad_inches=0)

    pressure_field_fig_path = "/".join([img_dir_path, "pressure_field.png"])
    fig_p = plot_fenics_object(p)
    plt.savefig(pressure_field_fig_path, bbox_inches="tight", pad_inches=0)

    # zip for user download
    mesh_dir_path = "/".join([results_path, "meshes"])
    if not os.path.isdir(mesh_dir_path):
        os.mkdir(mesh_dir_path)

    abs_folders_in_zip = [
        img_dir_path,
        paraview_dir_path,
        paraview_mesh_dir_path,
        csv_dir_path,
    ]
    rel_folders_in_zip = [
        img_folder_name,
        results_vizu_folder_name,
        mesh_vizu_folder_name,
        results_folder_name,
    ]

    zipObj_path = "/".join([results_path, "simulation_results.zip"])
    with ZipFile(zipObj_path, "w") as zipObj:
        for abs_folderName, rel_folderName in zip(abs_folders_in_zip, rel_folders_in_zip):
            for filename in os.listdir(abs_folderName):
                # create complete filepath of file in directory
                filePath = os.path.join(abs_folderName, filename)

                # Add file to zip
                zipObj.write(filePath, os.path.join(rel_folderName, filename))
