import fenics as fn
import matplotlib.pyplot as plt

# Import the mesh
mesh = fn.Mesh('/home/sabiqislam/Freecad_files/Assembly_DW-Compound.xml')

# Import subdomains and boundaries
markers = fn.MeshFunction("size_t", mesh, '/home/sabiqislam/Freecad_files/Assembly_DW-Compound_gmsh_physical.xml')
boundaries = fn.MeshFunction("size_t", mesh, '/home/sabiqislam/Freecad_files/Assembly_DW-Compound_gmsh_geometrical.xml')

# Define function space
V = fn.FunctionSpace(mesh, 'P', 1)

# Silicon nitride permittivity
epsilon_r_dielectric = 7.5
epsilon_0 = 8.854e-12  # Vacuum permittivity (F/m)

# Define permittivity function
class Permittivity(fn.UserExpression):
    def __init__(self, markers, **kwargs):
        super().__init__(**kwargs)
        self.markers = markers

    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 61:  # ID for silicon nitride
            values[0] = epsilon_r_dielectric * epsilon_0
        else:
            values[0] = epsilon_0  

    def value_shape(self):
        return ()

eps = Permittivity(markers, degree=2)

# Boundary conditions
top_wire_bc = fn.DirichletBC(V, fn.Constant(1.0), boundaries, 63)
bottom_wire_bc = fn.DirichletBC(V, fn.Constant(0.0), boundaries, 62)
bcs = [top_wire_bc, bottom_wire_bc]


u = fn.TrialFunction(V)
v = fn.TestFunction(V)
a = fn.dot(fn.grad(u), fn.grad(v)) * eps * dx
L = fn.Constant('0') * v * dx


u_solution = fn.Function(V)
fn.solve(a == L, u_solution, bcs)

# Plot
plt.figure()
p = fn.plot(u_solution)
plt.colorbar(p)
plt.title('Electric Potential')
plt.show()

# Calculate and plot the electric field
E = fn.project(-fn.grad(u_solution), V)
plt.figure()
p = fn.plot(E)
plt.colorbar(p)
plt.title('Electric Field')
plt.show()
