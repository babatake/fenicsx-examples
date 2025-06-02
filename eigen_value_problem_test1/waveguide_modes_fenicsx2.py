import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import mesh, fem
from dolfinx.fem.petsc import assemble_matrix
import dolfinx.plot
import ufl
from ufl import grad, dx
from slepc4py import SLEPc
import pyvista

# ---- パラメータ定義 ----
Lx, Ly = 1.0, 0.5        # 導波路のサイズ
nx, ny = 40, 20          # メッシュ分割数
num_modes = 5            # 計算するモード数

# ---- メッシュ作成 ----
domain = mesh.create_rectangle(MPI.COMM_WORLD,
                               points=[[0, 0], [Lx, Ly]],
                               n=[nx, ny],
                               cell_type=mesh.CellType.triangle)

# ---- 関数空間定義 ----
V = fem.functionspace(domain, ("CG", 2))

# ---- 異方性誘電率定義（例：定数） ----
eps_xx = fem.Function(V)
eps_yy = fem.Function(V)
eps_xx.x.array[:] = 4.0  # x方向の誘電率
eps_yy.x.array[:] = 2.0  # y方向の誘電率

# ---- 試行関数とテスト関数 ----
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
grad_u = grad(u)
grad_v = grad(v)

# ---- 剛性行列（異方性考慮） ----
a = eps_xx * grad_u[0] * grad_v[0] * dx + eps_yy * grad_u[1] * grad_v[1] * dx

# ---- 質量行列（通常） ----
b = u * v * dx

# ---- PEC境界条件の適用 ----
def pec_boundary(x):
    return np.isclose(x[0], 0) | np.isclose(x[0], Lx) | np.isclose(x[1], 0) | np.isclose(x[1], Ly)

bc = fem.dirichletbc(fem.Function(V), fem.locate_dofs_geometrical(V, pec_boundary))

# ---- 行列アセンブル ----
A = assemble_matrix(fem.form(a), bcs=[bc])
A.assemble()

B = assemble_matrix(fem.form(b), bcs=[bc])
B.assemble()

# ---- SLEPc 固有値解析 ----
eps = SLEPc.EPS().create()
eps.setOperators(A, B)
eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
eps.setDimensions(num_modes)
eps.setFromOptions()
eps.solve()

# ---- 結果取得 ----
nconv = eps.getConverged()
print(f"Number of converged eigenpairs: {nconv}")

Hz_modes = []
for i in range(min(nconv, num_modes)):
    r = eps.getEigenvalue(i)
    eigvec = A.getVecLeft()
    eps.getEigenvector(i, eigvec)

    Hz = fem.Function(V)
    Hz.x.array[:] = eigvec[:]
    Hz_modes.append((r, Hz))

    print(f"Mode {i}: k_t^2 = {r:.5f}")

# ---- モードの可視化 ----
plotter = pyvista.Plotter(shape=(1, num_modes))
for i in range(num_modes):
    r, Hz = Hz_modes[i]
    grid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(V))
    grid.point_data["Hz"] = Hz.x.array.real

    plotter.subplot(0, i)
    plotter.add_mesh(grid, scalars="Hz", cmap="coolwarm", show_edges=False)
    plotter.add_text(f"Mode {i}\nk_t²={r:.2f}", font_size=10)
    plotter.view_xy()

plotter.link_views()
plotter.show()
