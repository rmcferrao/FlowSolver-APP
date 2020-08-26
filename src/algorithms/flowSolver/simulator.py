from dolfin import *
import numpy as np
from skimage import io
import sys
import matplotlib.pyplot as plt


# insert at 1, 0 is the script path (or '' in REPL)
# sys.path.insert(1, '/home/rafael/Documents/web-app/python-code/Fenics-Simulation/plot')
# from plot_class import *

def structure_solver(rect_mesh_path, boundaries_path, subdomains_path, flowEquation, pin, pout, mu, f = (0, 0), k = None):
    parameters['dof_ordering_library'] = 'Boost'

    # print(parameters['dof_ordering_library'])
    rect_mesh = Mesh(rect_mesh_path)

    mesh_width = rect_mesh.coordinates()[:,0].max()
    mesh_height = rect_mesh.coordinates()[:,1].max()

    subdomains = MeshFunction('size_t', rect_mesh, subdomains_path)
    boundaries = MeshFunction('size_t', rect_mesh, boundaries_path)
    # boundaries.set_all(0)

    mu = Constant(mu)
    pin = Constant(pin)
    pout = Constant(pout)
    f = Constant(f)
    if k:
        k = Constant(k)

    # # Sub domain for inflow (right)
    # class Outflow(SubDomain):
    #     def inside(self, x, on_boundary):
    #         return x[0] > mesh_width - DOLFIN_EPS and on_boundary

    # # Sub domain for outflow (left)
    # class Inflow(SubDomain):
    #     def inside(self, x, on_boundary):
    #         return x[0] < DOLFIN_EPS and on_boundary

    # # Sub Y 0
    # class Y0(SubDomain):
    #     def inside(self, x, on_boundary):
    #         return x[1] < DOLFIN_EPS and on_boundary

    # # Sub Y L
    # class Yl(SubDomain):
    #     def inside(self, x, on_boundary):
    #         return x[1] > mesh_height- DOLFIN_EPS and on_boundary

    # inflow, outflow = Inflow(), Outflow()
    # y0, yl = Y0(), Yl()

    # inflow.mark(boundaries, 1)
    # outflow.mark(boundaries, 2)
    # y0.mark(boundaries, 3)
    # yl.mark(boundaries, 4)
    
    # if flowEquation == 'Stokes':
    #     for facet in facets(rect_mesh):
    #         i_facet = facet.index()
    #         cell_is = [cell.index() for cell in cells(facet)]
    #         if len(cell_is) == 1:
    #             # print(boundaries.array()[i_facet])
    #             if boundaries.array()[i_facet] == 0:
    #                 boundaries.array()[i_facet] = 5
                
    # File('/home/rafael/Desktop/boudaries.pvd') <<  boundaries

    V = VectorElement('P', rect_mesh.ufl_cell(), 2)
    Q = FiniteElement('P', rect_mesh.ufl_cell(), 1)
    Element = MixedElement([V, Q])

    W = FunctionSpace(rect_mesh, Element)
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)

    # boundary conditions
    if flowEquation == 'Stokes':
        bc1 = DirichletBC(W.sub(0),Constant((0.0,0.0)),boundaries,3)
        bc2 = DirichletBC(W.sub(0),Constant((0.0,0.0)),boundaries,4)
        bc3 = DirichletBC(W.sub(0),Constant((0.0,0.0)),boundaries,5)
    
    elif flowEquation == 'Brinkman':
        bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,3)
        bc2 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,4)
        

    # Demarcacao dos indices das regioes
    ds = Measure('ds', domain=rect_mesh, subdomain_data=boundaries)
    dx = Measure('dx', domain=rect_mesh, subdomain_data=subdomains)

    # vetor normal a malha
    n = FacetNormal(rect_mesh)

    # Variational Forumalation and Results
    aD = lambda  u, v: inner(u, v)
    aS = lambda u, v: inner(grad(u), grad(v))
    b = lambda u, q: div(u) * q
    lf = lambda f, v: inner(f, v)

    # 1 rock matrix - darcy or null
    # 2 free flow - stokes
    if flowEquation == 'Stokes':
        a = mu * aS(u,v)*dx - b(v,p)*dx - b(u,q)*dx
        l = lf(f,v)*dx - pin*dot(v,n)*ds(1) - pout*dot(v,n)*ds(2)

    elif flowEquation == 'Brinkman':
        a = mu * aS(u,v)*dx(2) +(mu/k)*aD(u,v)*dx(1) - b(v,p)*dx - b(u,q)*dx
        l = lf(f,v)*dx - pin*dot(v,n)*ds(1) - pout*dot(v,n)*ds(2)

    A = assemble(a)
    L = assemble(l)

    bc1.apply(A,L)
    bc2.apply(A,L)
    if flowEquation == 'Stokes':
        bc3.apply(A,L)

    w = Function(W)

    try:
        solve(A, w.vector(), L, 'mumps', 'none')
    except:
        solve(A, w.vector(), L)

    u,p = w.split()

    return u, p#, boundaries, subdomains
    

# rect_mesh_path = '/home/rafael/Documents/web-app/python-code/Fenics-Simulation/temp/meshes/meshFine.xml'
# subdomains_path = '/home/rafael/Documents/web-app/python-code/Fenics-Simulation/temp/meshes/subdomainFine.xml'
# flowEquation = 'Stokes'
# pin = 1
# pout = -1
# mu = 1e-3

# u, p, boundaries, subdomains = stucture_solver(rect_mesh_path, subdomains_path, flowEquation, pin, pout, mu)

# fig_u = plot_fenics_object(u)
# plt.savefig('/home/rafael/Documents/web-app/python-code/Fenics-Simulation/temp/images/velocity_field.png', bbox_inches="tight", pad_inches=0)

# fig_p = plot_fenics_object(p)
# plt.savefig('/home/rafael/Documents/web-app/python-code/Fenics-Simulation/temp/images/pressure_field.png', bbox_inches="tight", pad_inches=0)

# File('boundaries.pvd') << boundaries
# File('subdomains.pvd') << subdomains
