# ------------------------------------------------------------------------------
#
#  Gmsh Python 
#
# ------------------------------------------------------------------------------

import gmsh
import sys
import json
from subprocess import call

gmsh.initialize()

xlen = 20
ylen = 3
zlen = 1
meshname = 'x_35'
meshname_withformat = meshname + '.msh2'
lc = 1e-2
p1 = gmsh.model.geo.addPoint(0, 0, zlen)
p2 = gmsh.model.geo.addPoint(0, 0, 0)
p3 = gmsh.model.geo.addPoint(0, ylen, zlen)
p4 = gmsh.model.geo.addPoint(0, ylen, 0)

p5 = gmsh.model.geo.addPoint(xlen, 0, zlen)
p6 = gmsh.model.geo.addPoint(xlen, 0, 0)
p7 = gmsh.model.geo.addPoint(xlen, ylen, zlen)
p8 = gmsh.model.geo.addPoint(xlen, ylen, 0)



gmsh.model.geo.addLine(p2, p1, 1)
gmsh.model.geo.addLine(p1, p3, 2)
gmsh.model.geo.addLine(p4, p3, 3)
gmsh.model.geo.addLine(p2, p4, 4)

gmsh.model.geo.addLine(p6, p5, 5)
gmsh.model.geo.addLine(p5, p7, 6)
gmsh.model.geo.addLine(p8, p7, 7)
gmsh.model.geo.addLine(p6, p8, 8)

gmsh.model.geo.addLine(p2, p6, 9)
gmsh.model.geo.addLine(p1, p5, 10)
gmsh.model.geo.addLine(p4, p8, 11)
gmsh.model.geo.addLine(p3, p7, 12)

gmsh.model.geo.addCurveLoop([1, 2, -3, -4], 111)
gmsh.model.geo.addCurveLoop([5, 6, -7, -8], 112)
gmsh.model.geo.addCurveLoop([9, 5, -10, -1], 113)
gmsh.model.geo.addCurveLoop([11, 7, -12, -3], 114)
gmsh.model.geo.addCurveLoop([4, 11, -8, -9], 115)
gmsh.model.geo.addCurveLoop([2, 12, -6, -10], 116)

s1 = gmsh.model.geo.addPlaneSurface([111])
s2 = gmsh.model.geo.addPlaneSurface([112])
s3 = gmsh.model.geo.addPlaneSurface([113])
s4 = gmsh.model.geo.addPlaneSurface([114])
s5 = gmsh.model.geo.addPlaneSurface([115])
s6 = gmsh.model.geo.addPlaneSurface([116])


gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6], 1)
v1=gmsh.model.geo.addVolume([1])
gmsh.model.geo.synchronize()

pl1 =gmsh.model.addPhysicalGroup(2, [s1])
gmsh.model.setPhysicalName(2, pl1, "dirichlet")

pl2 =gmsh.model.addPhysicalGroup(2, [s2])
gmsh.model.setPhysicalName(2, pl2, "neumann")

ps = gmsh.model.addPhysicalGroup(3, [1])
gmsh.model.setPhysicalName(3, ps, "material")

# Set this to True to build a fully hex mesh:
transfinite = True
#transfinite = False
transfiniteAuto = False

if transfinite:
    NN = 9
    for c in gmsh.model.getEntities(1):
        gmsh.model.mesh.setTransfiniteCurve(c[1], NN)

    for s in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(s[1])
        gmsh.model.mesh.setRecombine(s[0], s[1])
        gmsh.model.mesh.setSmoothing(s[0], s[1], 100)
    gmsh.model.mesh.setTransfiniteVolume(v1)
elif transfiniteAuto:
    gmsh.option.setNumber('Mesh.MeshSizeMin', 0.5)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.5)
    # setTransfiniteAutomatic() uses the sizing constraints to set the number
    # of points
    gmsh.model.mesh.setTransfiniteAutomatic()
else:
    gmsh.option.setNumber('Mesh.MeshSizeMin', 0.05*100000)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.05*100000)

gmsh.model.mesh.generate(3)
gmsh.write(meshname_withformat)

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

# ------------------------------------------------------------------------------
#
#  AMfe Python 
#
# ------------------------------------------------------------------------------

import matplotlib.pyplot as plt
from amfe.component import StructuralComponent
from amfe.material import KirchhoffMaterial
from amfe.forces import constant_force
import numpy as np
from amfe import ui
import scipy.sparse as sparse
from subprocess import call
import os
import time

E_alu = 210e9
nu_alu = 0.3
rho_alu = 2.7e-3
my_material = KirchhoffMaterial(E_alu, nu_alu, rho_alu)
force           = constant_force(10000)
force_direction = np.array([0.0, 1.0, 0.0])
alpha = 1e-6
beta  = 1e-3
print('############################################################')
print('################### CREATING COMPONENTS ####################')
print('############################################################')
components_list = []
components_ids_list  = []
input_file = meshname +  '.msh2'
input_file_rename = meshname +  '.msh'
os.rename(input_file,input_file_rename)

my_mesh = ui.import_mesh_from_file(input_file_rename)
my_component = StructuralComponent(my_mesh)
ui.assign_material_by_group(my_component, my_material, 'material')
if 'dirichlet' in my_mesh.groups.keys():
    ui.set_dirichlet_by_group(my_component, 'dirichlet', ('ux'), 'Dirichlet_x')
    ui.set_dirichlet_by_group(my_component, 'dirichlet', ('uy'), 'Dirichlet_y')
    ui.set_dirichlet_by_group(my_component, 'dirichlet', ('uz'), 'Dirichlet_z')
if 'neumann' in my_mesh.groups.keys():
    ui.set_neumann_by_group(my_component, 'neumann', force_direction, False, 'Load', force)
system, formulation = ui.create_constrained_mechanical_system_from_component(my_component, all_linear=True)

no_of_dofs = system.dimension
q0 = np.zeros(no_of_dofs)
dq0 = q0
ddq0 = dq0
K1 = system.K(q0, dq0, 0)
M1 = system.M(q0, dq0, 0)
fext = system.f_ext(q0, dq0, 0)


nodecoord_reference = [8000,3000,1000]
nodesdf_forx= my_mesh.nodes_df[my_mesh.nodes_df['x'] == nodecoord_reference[0]]
nodes_forxy = nodesdf_forx[nodesdf_forx['y'] == nodecoord_reference[1]]
nodes_forxyz = nodes_forxy[nodes_forxy['z'] == nodecoord_reference[2]].index.tolist()
reference_xdof = my_component.mapping.nodal2global.loc[nodes_forxyz].values[0][0]
reference_ydof = my_component.mapping.nodal2global.loc[nodes_forxyz].values[0][1]
reference_zdof = my_component.mapping.nodal2global.loc[nodes_forxyz].values[0][2]
print('############################################################')
print('##################### Getting into FRF #####################')
print('############################################################')

class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
    def __call__(self, rk=None):
        self.niter += 1
        #if self._disp:
            #print('iter %3i\trk = %s' % (self.niter, str(rk)))
def cm2inch(value):
    return value/2.54
def inch2cm(value):
    return value*2.54
################################
################################
count = 0
update_ilu = True

time_taken_ilu_solver = 0.0
time_taken_ilu_precon = 0.0
time_taken_ilu_precon_list = []
time_taken_ilu_solver_list = []
time_taken_ilu_diff_list = []
time_taken_ilu = 0.0

time_taken_ilu_full_solver = 0.0
time_taken_ilu_full_precon = 0.0
time_taken_ilu_full_precon_list = []
time_taken_ilu_full_solver_list = []
time_taken_ilu_full_diff_list = []
time_taken_ilu_full = 0.0

number_iterations_ilu_list = []
number_iterations_full_ilu_list = []
number_update_ilu = 0.0
u_solution = dict()
u_solution_dof = []
freq_list = []
################################
fc = fext + 0.0 * 1J * fext
buildPrec = lambda w, M1, K1, alpha=1e-6, beta=1e-3: -w ** 2 * M1 + K1 + 1J * w * (alpha * K1 + beta * M1)
buildZ = lambda w, M, K, alpha=1e-6, beta=1e-3: -w ** 2 * M + K + 1J * w * (alpha * K + beta * M)
tolere = 1e-3
w_list = np.linspace(100.0,10000, 150)
initial_guess = fc * 0.0



for i in range(2 * len(w_list)):
    w = w_list[count]
    Z = buildZ(w, M1, K1, alpha, beta)
    Z1 = buildPrec(w, M1, K1, alpha, beta)

    #################################################################3
    begin_precon = time.time()
    M_ilu_full = sparse.linalg.spilu(Z1, drop_tol=tolere)
    Mx_ilu_full = lambda x: M_ilu_full.solve(x)
    M3_prec_full = sparse.linalg.LinearOperator((K1.shape[0], K1.shape[0]), Mx_ilu_full)
    end_precon = time.time()

    begin_solver = time.time()
    counter_full_ilu = gmres_counter()
    usol_ilu_full, info_ilu_full = sparse.linalg.bicgstab(Z, fc, x0=initial_guess, M=M3_prec_full, tol=1e-07,
                                                          callback=counter_full_ilu)
    end_solver = time.time()
    time_taken_ilu_full_precon = end_precon - begin_precon
    time_taken_ilu_full_precon_list.append(time_taken_ilu_full_precon)

    time_taken_ilu_full_solver = end_solver - begin_solver
    time_taken_ilu_full_solver_list.append(time_taken_ilu_full_solver)

    number_iterations_full_ilu_list.append(counter_full_ilu.niter)
    time_taken_ilu_full = time_taken_ilu_full + time_taken_ilu_full_solver + time_taken_ilu_full_precon
    time_taken_ilu_full_diff_list.append(
        end_solver - begin_precon - time_taken_ilu_full_solver - time_taken_ilu_full_precon)
    #################################################################3
    usol_recovered, du, ddu = formulation.recover(usol_ilu_full, dq0, ddq0, 0)
    u_solution.update({w:usol_recovered})
    u_solution_dof.append(np.abs(usol_recovered[reference_ydof]))
    freq_list.append(w/(2*np.pi))
    count += 1

    if w == w_list[-1]:
        break

plt.figure()
plt.plot(freq_list,u_solution_dof,'-')
plt.yscale('log')
plt.show()
