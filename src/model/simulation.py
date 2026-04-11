"""
Module for running the advanced finite element analysis solver for the furnace simulation.
Sets up the mesh, physics, executes the solver, and outputs visualizations and XDMF data.
"""

import numpy as np
from mpi4py import MPI
from dolfinx import io, fem, mesh, default_scalar_type, geometry
from dolfinx.fem.petsc import LinearProblem
import ufl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
from pathlib import Path

def run_advanced_solver() -> None:
    """
    Executes the main finite element simulation.

    Reads the provided XDMF domain and boundary tags, assembles the steady-state
    thermal conduction and convection problem, solves it, and exports
    2D heatmap visualizations as well as full 3D XDMF outputs.
    """
    print("0. Architecting Directory Structure...")
    base_dir = Path.cwd() / "data"
    
    dirs = {
        "meshes": base_dir / "meshes",
        "xdmf": base_dir / "xdmf-s",
        "graphs": base_dir / "graphs",
        "visuals": base_dir / "visuals"
    }
    
    for folder in dirs.values():
        folder.mkdir(parents=True, exist_ok=True)
    
    print("1. Loading XDMF Mesh and Material Tags...")
    mesh_path = dirs["xdmf"] / "furnace_domain.xdmf"
    
    with io.XDMFFile(MPI.COMM_WORLD, str(mesh_path), "r") as xdmf:
        domain = xdmf.read_mesh(name="Grid")
        cell_tags = xdmf.read_meshtags(domain, name="Grid")

    domain.geometry.x[:, :] /= 1000.0

    print("2. Configuring Physics...")
    V = fem.functionspace(domain, ("Lagrange", 1))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    dx = ufl.Measure("dx", domain=domain, subdomain_data=cell_tags)

    print("3. Isolating Boundaries (Box vs. Tube)...")
    box_facets = mesh.locate_entities_boundary(
        domain, domain.topology.dim - 1, 
        lambda x: np.logical_or(np.abs(x[0]) >= 0.23, np.abs(x[1]) >= 0.23)
    )
    tube_facets = mesh.locate_entities_boundary(
        domain, domain.topology.dim - 1, 
        lambda x: np.logical_and(np.abs(x[0]) < 0.23, np.abs(x[1]) < 0.23)
    )

    facet_indices = np.concatenate([box_facets, tube_facets])
    facet_values = np.concatenate([np.full(len(box_facets), 1, dtype=np.int32), 
                                   np.full(len(tube_facets), 2, dtype=np.int32)])
    
    sort_order = np.argsort(facet_indices)
    facet_tags = mesh.meshtags(domain, domain.topology.dim - 1, 
                               facet_indices[sort_order], facet_values[sort_order])
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

    print("4. Building Effective Property Matrix...")
    k_values = {
        1: 150.0, 2: 1.4, 3: 0.1, 4: 0.026, 5: 1.0, 6: 5.0, 7: 0.08, 8: 15.0
    }

    print("5. Defining Spatial Buoyancy Gradients...")
    x = ufl.SpatialCoordinate(domain)
    y = x[1] 
    
    h_box = 8.5 - 6.5 * (y / 0.24) 
    h_tube = default_scalar_type(2.0) 
    T_amb = default_scalar_type(20.0)

    a = 0
    for tag, k in k_values.items():
        a += k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx(tag)
    
    a += h_box * u * v * ds(1)
    a += h_tube * u * v * ds(2)

    print("6. Applying Volumetric Heat Generation to Coiled Coil (Tag 6)...")
    Q = fem.Constant(domain, default_scalar_type(168960.0))
    L = Q * v * dx(6)
    L += h_box * T_amb * v * ds(1)
    L += h_tube * T_amb * v * ds(2)

    print("7. Executing Matrix Solver...")
    problem = LinearProblem(a, L, bcs=[], petsc_options={"ksp_type": "cg", "pc_type": "gamg"}, petsc_options_prefix="furnace_solve_")
    T_solution = problem.solve()
    T_solution.name = "Temperature"

    print(f"   -> Max Temp: {T_solution.x.array.max():.1f} °C")

    np.save(dirs["xdmf"] / "T_array.npy", T_solution.x.array)

    print("8. Probing Mesh and Generating 2D Annotated Heatmaps...")
    plt.style.use('dark_background')

    bb_tree = geometry.bb_tree(domain, domain.topology.dim)

    print("   -> Rendering Camera 1: Longitudinal (Length) Slice...")
    z_long = np.linspace(-0.42, 0.42, 400) # Tube is 800mm long
    y_long = np.linspace(-0.26, 0.26, 250) # Box is 480mm high
    Z_grid_L, Y_grid_L = np.meshgrid(z_long, y_long)
    X_fixed_L = np.zeros_like(Z_grid_L)
    
    points_L = np.vstack((X_fixed_L.flatten(), Y_grid_L.flatten(), Z_grid_L.flatten()))
    temp_L = np.full(points_L.shape[1], np.nan)
    
    cells_L, points_on_proc_L, valid_L = [], [], []
    colliding_cells_L = geometry.compute_colliding_cells(domain, geometry.compute_collisions_points(bb_tree, points_L.T), points_L.T)
    for i, p in enumerate(points_L.T):
        if len(colliding_cells_L.links(i)) > 0:
            cells_L.append(colliding_cells_L.links(i)[0])
            points_on_proc_L.append(p)
            valid_L.append(i)
            
    if len(points_on_proc_L) > 0:
        temp_L[valid_L] = T_solution.eval(points_on_proc_L, cells_L).flatten()
        
    fig1, ax1 = plt.subplots(figsize=(14, 8))
    heatmap1 = ax1.pcolormesh(Z_grid_L * 1000, Y_grid_L * 1000, temp_L.reshape(Z_grid_L.shape), cmap='inferno', shading='auto')
    cbar1 = plt.colorbar(heatmap1, ax=ax1, pad=0.02)
    cbar1.set_label('Temperature (°C)', fontsize=14, fontweight='bold')
    
    ax1.add_patch(patches.Rectangle((-240, -240), 480, 480, linewidth=2, edgecolor='white', facecolor='none', alpha=0.5, linestyle='--'))
    ax1.text(-230, 220, "Steel Casing", color='white', fontsize=12, fontweight='bold')
    
    ax1.add_patch(patches.Rectangle((-400, -35), 800, 70, linewidth=1.5, edgecolor='cyan', facecolor='none', alpha=0.6, linestyle='-'))
    ax1.text(260, 45, "Quartz Tube", color='cyan', fontsize=12, fontweight='bold')
    
    ax1.add_patch(patches.Rectangle((-225, -58), 450, 116, linewidth=1.5, edgecolor='lime', facecolor='none', alpha=0.8, linestyle=':'))
    ax1.text(-70, 70, "Heating Zone", color='lime', fontsize=12, fontweight='bold')
    
    ax1.set_xlabel("Length along Z-Axis (mm)", fontsize=12)
    ax1.set_ylabel("Height along Y-Axis (mm)", fontsize=12)
    ax1.set_title("Longitudinal Cross-Section (YZ Plane at X=0)", fontsize=16, fontweight='bold', pad=15)
    ax1.set_aspect('equal')
    plt.tight_layout()
    fig1.savefig(str(dirs["visuals"] / "furnace_lengthwise_slice.png"), dpi=300, facecolor='black')
    plt.close(fig1)

    print("   -> Rendering Camera 2: Perpendicular (Barrel) Slice...")
    x_perp = np.linspace(-0.26, 0.26, 250) 
    y_perp = np.linspace(-0.26, 0.26, 250) 
    X_grid_P, Y_grid_P = np.meshgrid(x_perp, y_perp)
    Z_fixed_P = np.zeros_like(X_grid_P)
    
    points_P = np.vstack((X_grid_P.flatten(), Y_grid_P.flatten(), Z_fixed_P.flatten()))
    temp_P = np.full(points_P.shape[1], np.nan)
    
    cells_P, points_on_proc_P, valid_P = [], [], []
    colliding_cells_P = geometry.compute_colliding_cells(domain, geometry.compute_collisions_points(bb_tree, points_P.T), points_P.T)
    for i, p in enumerate(points_P.T):
        if len(colliding_cells_P.links(i)) > 0:
            cells_P.append(colliding_cells_P.links(i)[0])
            points_on_proc_P.append(p)
            valid_P.append(i)
            
    if len(points_on_proc_P) > 0:
        temp_P[valid_P] = T_solution.eval(points_on_proc_P, cells_P).flatten()
        
    fig2, ax2 = plt.subplots(figsize=(10, 10))
    heatmap2 = ax2.pcolormesh(X_grid_P * 1000, Y_grid_P * 1000, temp_P.reshape(X_grid_P.shape), cmap='inferno', shading='auto')
    cbar2 = plt.colorbar(heatmap2, ax=ax2, pad=0.02, fraction=0.046)
    cbar2.set_label('Temperature (°C)', fontsize=14, fontweight='bold')
    
    ax2.add_patch(patches.Rectangle((-240, -240), 480, 480, linewidth=2, edgecolor='white', facecolor='none', alpha=0.5, linestyle='--'))
    ax2.text(-230, 220, "Steel Casing", color='white', fontsize=12, fontweight='bold')
    
    ax2.add_patch(patches.Circle((0, 0), radius=58, linewidth=1.5, edgecolor='lime', facecolor='none', alpha=0.8, linestyle=':'))
    ax2.text(65, 55, "Heating Zone", color='lime', fontsize=12, fontweight='bold')
    
    ax2.add_patch(patches.Circle((0, 0), radius=35, linewidth=1.5, edgecolor='cyan', facecolor='none', alpha=0.8, linestyle='-'))
    ax2.text(-40, -10, "Quartz Tube", color='cyan', fontsize=12, fontweight='bold')
    
    ax2.set_xlabel("Width along X-Axis (mm)", fontsize=12)
    ax2.set_ylabel("Height along Y-Axis (mm)", fontsize=12)
    ax2.set_title("Perpendicular Cross-Section (XY Plane at Z=0)", fontsize=16, fontweight='bold', pad=15)
    ax2.set_aspect('equal')
    plt.tight_layout()
    fig2.savefig(str(dirs["visuals"] / "furnace_barrel_slice.png"), dpi=300, facecolor='black')
    plt.close(fig2)

    plt.style.use('default')
    print("-> Both cross-sections saved successfully to the 'visuals' folder!")

    print("9. Exporting XDMF for ParaView...")
    output_path = dirs["xdmf"] / "furnace_temperature.xdmf"
    with io.XDMFFile(MPI.COMM_WORLD, str(output_path), "w") as out_file:
        out_file.write_mesh(domain)
        out_file.write_function(T_solution)

    print(f"Done! Find your results in {base_dir}")

if __name__ == "__main__":
    run_advanced_solver()