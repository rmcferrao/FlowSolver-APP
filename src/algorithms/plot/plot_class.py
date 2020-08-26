from dolfin import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri

def plot_subdomain(subdomains, edgecolor='none', linewidth='none'):
    x = subdomains.mesh().coordinates()[:, 0]
    y = subdomains.mesh().coordinates()[:, 1]
    t = subdomains.mesh().cells()
    i = subdomains.array()

    fig = plt.figure()

    plt.tripcolor(x, y, t, facecolors=i, edgecolors=edgecolor, cmap="binary_r", vmax=3, linewidth=linewidth)
    plt.axis("off")
    plt.axis('equal')
    fig.tight_layout()

    return fig

def mesh2triang(mesh):
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

def plot_fenics_object(obj):

    fig = plt.figure()

    plt.gca().set_aspect('equal')

    # pressure
    mesh = obj.function_space().mesh()
    C = obj.compute_vertex_values(mesh)
    cmap="plasma"

    plt.xlabel('x direction', fontsize=16)
    plt.ylabel('y direction', fontsize=16)
    plt.title('Pressure Field', fontsize=22)
    color_bar_label = 'Pressure (Pa)'

    if not C.size == mesh.coordinates().shape[0]:
        #half/half -> velocity
        cmap='jet'
        C1 = C[:int(len(C)/2)]
        C2 = C[int(len(C)/2):]

        color_bar_label = 'Velocity Magnitude\n(pixels/s)'
        plt.title('Velocity Field', fontsize=22)

        C = (C1 **2 + C2 ** 2) ** (1/2)
    
    plt.tripcolor(mesh2triang(mesh), C, cmap=cmap, shading='gouraud')
    plt.colorbar(shrink=0.6).set_label(color_bar_label)
    
    plt.axis('equal')
    fig.tight_layout()
    return fig
