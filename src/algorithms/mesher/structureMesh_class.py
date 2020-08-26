from dolfin import *
import numpy as np
import skimage.io as sk_io
import sys
import matplotlib.pyplot as plt
import matplotlib.tri.triangulation as triangulation

# sys.path.insert(1, '/home/rafael/Documents/web-app/python-code/Fenics-Simulation/plot')
# from plot_class import *

def recursive_sum(res, coords, loc, len):
    for i in range(len - 1):
        res = np.vstack((res, (coords[loc] + coords[i + 1]) / 2))
    len = len - 1
    loc = loc + 1
    if len == 0:
        return res
    else:
        return recursive_sum(res, coords, loc, len)

def image_indexes(cell, n_iteration):
    coords = []
    coords.append([cell.midpoint().x(), cell.midpoint().y()])
    coords.append(cell.get_vertex_coordinates()[:2])
    coords.append(cell.get_vertex_coordinates()[2:4])
    coords.append(cell.get_vertex_coordinates()[4:6])
    coords1 = np.array(coords)

    if n_iteration == 1:
        return coords1
        # return np.delete(coords1, index_1, axis=0)  # return middle point

    coords2 = recursive_sum(coords1, coords1, 0, coords1.shape[0])
    if n_iteration == 2:
        return coords2

    coords3 = recursive_sum(coords2, coords2, 0, coords2.shape[0])
    if n_iteration == 3:
        return coords3

    coords4 = recursive_sum(coords3, coords3, 0, coords3.shape[0])
    if n_iteration == 4:
        return coords4

def sum_coords_in(coords, coquina_array, pixels_len, free_flow_pixel):
    count = 0
    for coord in coords:
        point_image = (
            int(coord[0] * pixels_len[0]),
            int(coord[1] * pixels_len[1]),
        )
        if coquina_array[point_image] == free_flow_pixel:
            count += 1
    if (count / coords.shape[0]) > 0.5:
        return 1
    else:
        return 0

def get_neighbors_cells(cell):
    tdim = 2
    face_it = [facet.entities(tdim) for facet in facets(cell)]
    face_it = np.unique(np.concatenate(face_it)).tolist()
    return list(filter(lambda ci: ci != cell.index(), face_it))
    
def indice_of_isolated_cells(mesh, subdomains):
    #indexes of no flow type cells
    pixel_flow_indexes = np.concatenate(np.argwhere(subdomains.array()==2)).tolist()

    mesh.init(1,2)

    cell_i_del = []
    check_cell = []

    iterator = 0

    for iterator in range(len(pixel_flow_indexes)):
        index = pixel_flow_indexes[iterator]
                
        cell = Cell(mesh, index)
        cell_neighbors_indexes = get_neighbors_cells(cell)
        cell_neighbors_markers = subdomains.array()[cell_neighbors_indexes]

        if all(cell_neighbors_markers == 1):

            cell_i_del.append(index)
            # pixel_flow_indexes.remove(index)

        elif len(cell_neighbors_markers) == 3 and (cell_neighbors_markers == 1).sum() == 2:

            neighbor_freeflow_local = np.where(cell_neighbors_markers == 2)[0][0]
            neighbor_freeflow_index = cell_neighbors_indexes[neighbor_freeflow_local]

            cell = Cell(mesh, neighbor_freeflow_index)
            cell_neighbors_indexes = get_neighbors_cells(cell)

            cell_neighbors_indexes.remove(index)
            cell_neighbors_markers = subdomains.array()[cell_neighbors_indexes]

            if all(cell_neighbors_markers == 1):

                cell_i_del.append(index)
                # pixel_flow_indexes.remove(index)
                cell_i_del.append(neighbor_freeflow_index)
                # pixel_flow_indexes.remove(neighbor_freeflow_index)

        # iterator += 1

    if len(cell_i_del) == 0:
        return False
    else:
        return cell_i_del

def delete_cell(rect_mesh, subdomains, del_cells = False):
    old_cells = rect_mesh.cells()
    old_coors = rect_mesh.coordinates()
    # File('/home/rafael/Desktop/mesh.pvd') << subdomains
    if del_cells:
        cell_i_del_new = indice_of_isolated_cells(rect_mesh, subdomains)
        subdomains.array()[cell_i_del_new] = 1

    del_i_cells = np.where(subdomains.array() == 1)[0] # indice das celulas que serao deletadas
    keep_i_cells = np.where(subdomains.array() == 2)[0] # indice das celulas que serao mantidas
     
    del_i_coords = old_cells[del_i_cells] 
    del_i_coords_unique = np.unique(del_i_coords.flatten())

    keep_i_coords = old_cells[keep_i_cells] 
    keep_i_coords_unique = np.unique(keep_i_coords.flatten())

    intersect_i_coords = np.intersect1d(del_i_coords_unique,keep_i_coords_unique)

    inside_i_coords_boolean = np.invert(np.isin(del_i_coords_unique, intersect_i_coords, assume_unique=True))
    
    inside_i_coords = del_i_coords_unique[inside_i_coords_boolean] # indice das coordenadas que serao deletadas

    intersect_del_i_coords = np.invert(np.isin(np.array(list(range(len(old_coors)))),inside_i_coords))

    key_transform_dict = np.array(list(range(len(old_coors))))[intersect_del_i_coords]
    vals_transform_dict = np.array(list(range(len(key_transform_dict))))

    transform_dict = {key:val for key,val in zip(key_transform_dict, vals_transform_dict)} #dicionario o indice da coord que era aponta pro que vai

    new_cells = np.delete(old_cells, del_i_cells, axis=0)

    for i_cell in range(len(new_cells)):
        new_cells[i_cell,0] = transform_dict[new_cells[i_cell,0]]
        new_cells[i_cell,1] = transform_dict[new_cells[i_cell,1]]
        new_cells[i_cell,2] = transform_dict[new_cells[i_cell,2]]
    

    new_coords = np.delete(old_coors, inside_i_coords, axis=0)

    new_mesh = Mesh() 
    editor = MeshEditor()

    cell_type = 'triangle'
    editor.open(new_mesh,cell_type, 2, 2)
    editor.init_vertices(len(new_coords))
    editor.init_cells(len(new_cells))

    for vert_id, vert in enumerate(new_coords):
        editor.add_vertex(vert_id, vert)

    cell_id = 0
    for c in range(len(new_cells)): 
        editor.add_cell(cell_id, new_cells[c])
        cell_id += 1 

    editor.close()

    new_subdomains = MeshFunction("size_t", new_mesh, 2)
    new_subdomains.set_all(2)

    return new_mesh, new_subdomains

def mark_inside(square_mesh, coquina_array, n_iteration, free_flow_pixel, flowMethod, del_cells = False):
    pixels_len = np.zeros(2)

    # coquina_array = np.transpose(coquina_array)

    pixels_len[0] = coquina_array.shape[0]-1
    pixels_len[1] = coquina_array.shape[1]-1

    aspect_image = (pixels_len / pixels_len.max()).round(10)

    # square_mesh.coordinates()[:,0] *= aspect_image[1]
    # square_mesh.coordinates()[:,1] *= aspect_image[0]


    debug_text = "Mesh Size = {} x {}".format(square_mesh.coordinates()[:,0].max(), square_mesh.coordinates()[:,1].max())

    subdomains = MeshFunction("size_t", square_mesh, 2)
    subdomains.set_all(1)

    for c in cells(square_mesh):
        index = c.index()
        coords_in = image_indexes(c, n_iteration)
        color_cell = sum_coords_in(coords_in, coquina_array, pixels_len, free_flow_pixel)
        if color_cell == 1:
            subdomains.array()[index] = 2


    new_coordinates = np.zeros(square_mesh.coordinates().shape)
    new_coordinates[:, 0] = square_mesh.coordinates()[:, 1]
    new_coordinates[:, 1] = (
        square_mesh.coordinates()[:, 0].max() - square_mesh.coordinates()[:, 0]
    )
    square_mesh.coordinates()[:, :] = new_coordinates
    square_mesh.coordinates()[:,0] *= aspect_image[1]
    square_mesh.coordinates()[:,1] *= aspect_image[0]

    if flowMethod == 'Stokes':
        square_mesh, subdomains = delete_cell(square_mesh, subdomains, del_cells)

    debug1_text = "Mesh Size = {} x {}".format(square_mesh.coordinates()[:,0].max(), square_mesh.coordinates()[:,1].max())
    debug2_text = "SubMesh Size = {} x {}".format(subdomains.mesh().coordinates()[:,0].max(), subdomains.mesh().coordinates()[:,1].max())

    # print(debug1_text)
    # print(debug2_text)

    return square_mesh, subdomains

def mark_boundaries(rect_mesh, flowEquation):
 
    mesh_width = rect_mesh.coordinates()[:,0].max()
    mesh_height = rect_mesh.coordinates()[:,1].max()

    boundaries = MeshFunction('size_t', rect_mesh, 1)
    boundaries.set_all(0)

    # Sub domain for inflow (right)
    class OutflowBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] > mesh_width - DOLFIN_EPS and on_boundary

    # Sub domain for outflow (left)
    class InflowBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] < DOLFIN_EPS and on_boundary

    # Sub Y 0
    class TopBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] < DOLFIN_EPS and on_boundary

    # Sub Y L
    class BotBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] > mesh_height- DOLFIN_EPS and on_boundary

    inflow, outflow = InflowBoundary(), OutflowBoundary()
    botbound, topbound = BotBoundary(), TopBoundary()

    boundaries_draw = """  
                        __4___
       Fluid Flow      |     |
      ===========>   1 |     | 2
                       |_____|
                           3
    """

    # print('---------------------------------------')
    # text = "Total facets: {}\nmarked0: {}\nmarked1: {}\nmarked2: {}\nmarked3: {}\nmarked4: {}\nmarked5 {}".format(boundaries.array().size,np.where(boundaries.array() == 0)[0].size, np.where(boundaries.array() == 1)[0].size, np.where(boundaries.array() == 2)[0].size, np.where(boundaries.array() == 3)[0].size, np.where(boundaries.array() == 4)[0].size, np.where(boundaries.array() == 5)[0].size)
    # print(text)
    # print(boundaries_draw)
    # inflow.mark(boundaries, 1)
    # print('----------1-----------------------------')
    # outflow.mark(boundaries, 2)
    # print('---------2------------------------------')
    # botbound.mark(boundaries, 3)
    # print('-----------3----------------------------')
    # topbound.mark(boundaries, 4)
    # print('-------------4--------------------------')

    for facet in facets(rect_mesh):
        i_facet = facet.index()
        cell_is = [cell.index() for cell in cells(facet)]
        if len(cell_is) == 1:
            if facet.midpoint().x() < DOLFIN_EPS:
                boundaries.array()[i_facet] = 1

            elif facet.midpoint().x() > mesh_width - DOLFIN_EPS:
                boundaries.array()[i_facet] = 2

            elif facet.midpoint().y() < DOLFIN_EPS:
                boundaries.array()[i_facet] = 3
                
            elif facet.midpoint().y() > mesh_height - DOLFIN_EPS:
                boundaries.array()[i_facet] = 4
            else:
                boundaries.array()[i_facet] = 5

    # text = "Total facets: {}\nmarked1: {}\nmarked2: {}\nmarked3: {}\nmarked4: {}\nmarked5 {}".format(boundaries.array().size, np.where(boundaries.array() == 1)[0].size, np.where(boundaries.array() == 2)[0].size, np.where(boundaries.array() == 3)[0].size, np.where(boundaries.array() == 4)[0].size, np.where(boundaries.array() == 5)[0].size)
    # print(text)
    # print(boundaries_draw)
    # print('---------------------------------------')

    return boundaries
















# nx_elements, ny_elements = 100, 100
# square_mesh = UnitSquareMesh(nx_elements, ny_elements, diagonal="crossed")
# rock_image_path = "/home/rafael/Desktop/test.png" #input 1
# rock_array = sk_io.imread(rock_image_path)
# n_iteration = 1
# free_flow_pixel = 255
# flowEq = 'Stokes'

# square_mesh, subdomains = mark_inside(square_mesh, rock_array, n_iteration, free_flow_pixel, flowEq)

# edgecolor = "black"  # or 'none'
# linewidth = None # or 'none'
# fig = plot_subdomain(subdomains, edgecolor, linewidth)
# plt.savefig('/home/rafael/Desktop/mesh_to_fig.png', dpi=300)

# edgecolor = "blue"  # or 'none'
# linewidth = 0.05  # or 'none'
# fig = plot_subdomain(subdomains, edgecolor, linewidth)
# plt.savefig('/home/rafael/Desktop/mesh.png')

# File('/home/rafael/Desktop/mesh.pvd') << subdomains