import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import mesh, fem
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, create_vector
import dolfinx.plot 
import ufl
from slepc4py import SLEPc
import pyvista

# ----- パラメータ -----
Lx, Ly = 1.0, 0.5         # 導波路の幅と高さ
nx, ny = 40, 20           # メッシュ分割数
num_modes = 3             # 計算するモード数

# ----- メッシュ作成 -----３
domain = mesh.create_rectangle(MPI.COMM_WORLD,
                               points=[[0, 0], [Lx, Ly]],
                               n=[nx, ny],
                               cell_type=mesh.CellType.triangle)

# ----- 関数空間の定義 -----
V = fem.functionspace(domain, ("CG", 2))

# ----- 試行関数・テスト関数 -----
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# ----- バイリニア形式（剛性行列と質量行列） -----
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
b = ufl.inner(u, v) * ufl.dx

# ----- PEC境界条件 (H_z = 0) -----
def pec_boundary(x):
    return np.isclose(x[0], 0) | np.isclose(x[0], Lx) | np.isclose(x[1], 0) | np.isclose(x[1], Ly)

bc = fem.dirichletbc(value=fem.Function(V), dofs=fem.locate_dofs_geometrical(V, pec_boundary))

# ----- 行列アセンブル -----
A = assemble_matrix(fem.form(a), bcs=[bc])
A.assemble()

B = assemble_matrix(fem.form(b), bcs=[bc])
B.assemble()

# ----- 固有値問題の設定 -----
eps = SLEPc.EPS().create()
eps.setOperators(A, B)
eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)  # 一般化固有値問題
eps.setDimensions(num_modes)
eps.setFromOptions()
eps.solve()

# ----- 結果の表示 -----
nconv = eps.getConverged()
print(f"Number of converged eigenpairs: {nconv}")

Hz_modes = []
for i in range(min(nconv, num_modes)):
    r = eps.getEigenvalue(i)  # 固有値（k_t^2）

    eigvec_petsc = A.getVecLeft()
    eps.getEigenvector(i, eigvec_petsc)

    Hz = fem.Function(V)
    Hz.x.array[:] = eigvec_petsc[:]
    Hz_modes.append((r, Hz))


    print(f"Mode {i}: k_t^2 = {r:.5f}")


# ----- 可視化（最初のモード） -----
# 可視化用データ準備
grid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(V))
grid.point_data["Hz"] = Hz_modes[0][1].x.array.real

# 描画
plotter = pyvista.Plotter()
plotter.add_mesh(grid, show_edges=False, scalars="Hz", cmap="coolwarm")
plotter.add_title("TE Mode 0")
plotter.show()

plotter = pyvista.Plotter(shape=(1, num_modes), border=False)
for i in range(num_modes):
    r, Hz = Hz_modes[i]

    # グリッド再生成（各モードごとにコピー）
    grid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(V))
    grid.point_data["Hz"] = Hz.x.array.real

    # サブプロット位置を指定して追加
    plotter.subplot(0, i)
    plotter.add_mesh(grid,
                     scalars="Hz",
                     cmap="coolwarm",
                     show_edges=False)
    plotter.add_text(f"Mode {i}\nk_t² = {r:.3f}", font_size=10)
    plotter.view_xy()  # 真上視点で固定

plotter.link_views()  # カメラ同期（便利）
plotter.show()

