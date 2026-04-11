import dolfinx
from dolfinx import mesh, fem, io, geometry
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

def analyze_furnace():
    print("1. Loading FEniCSx Thermal Data (Numpy Injector)...")
    comm = MPI.COMM_WORLD
    
    with io.XDMFFile(comm, "data/xdmf-s/furnace_domain.xdmf", "r") as xdmf:
        domain = xdmf.read_mesh(name="Grid")
        
    bounds_max = domain.geometry.x.max(axis=0)
    r_max = float(bounds_max[0])
    z_max = float(bounds_max[2])
    
    is_meters = r_max < 1.0
    scale_to_mm = 1000.0 if is_meters else 1.0
    
    print(f"   -> Detected Mesh Bounds: Radius={r_max:.3f}, Length={z_max:.3f}")
    print(f"   -> Units detected as: {'Meters' if is_meters else 'Millimeters'}")

    V = fem.functionspace(domain, ("Lagrange", 1))
    T_solution = fem.Function(V)
    T_solution.x.array[:] = np.load("data/xdmf-s/T_array.npy")
    
    print("2. Initializing Probing Engine...")
    bb_tree = geometry.bb_tree(domain, domain.topology.dim)
    plt.style.use('dark_background')
    
    graphs_dir = Path("data/graphs")
    graphs_dir.mkdir(parents=True, exist_ok=True)

    def probe_mesh(points_array):
        temp_vals = np.full(points_array.shape[1], np.nan)
        cell_candidates = geometry.compute_collisions_points(bb_tree, points_array.T)
        colliding_cells = geometry.compute_colliding_cells(domain, cell_candidates, points_array.T)
        
        cells, pts_on_proc, valid = [], [], []
        for i, p in enumerate(points_array.T):
            if len(colliding_cells.links(i)) > 0:
                cells.append(colliding_cells.links(i)[0])
                pts_on_proc.append(p)
                valid.append(i)
                
        if len(pts_on_proc) > 0:
            temp_vals[valid] = T_solution.eval(pts_on_proc, cells).flatten()
        return temp_vals

    print("3. Generating Longitudinal Profile (Z-Axis)...")
    z_vals = np.linspace(-z_max, z_max, 500) 
    points_L = np.vstack((np.zeros_like(z_vals), np.zeros_like(z_vals), z_vals))
    T_long = probe_mesh(points_L)

    core_max_T = np.nanmax(T_long)
    tolerance = 10.0
    working_indices = np.where(T_long >= (core_max_T - tolerance))[0]
    
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(z_vals * scale_to_mm, T_long, color='magenta', linewidth=2, label="Centerline Temp")
    
    if len(working_indices) > 0:
        z_work_min = z_vals[working_indices[0]] * scale_to_mm
        z_work_max = z_vals[working_indices[-1]] * scale_to_mm
        ax1.axvspan(z_work_min, z_work_max, color='lime', alpha=0.2, label=f"Working Volume (±{tolerance}°C)")
        ax1.text(0, core_max_T - 50, f"Sweet Spot: {z_work_max - z_work_min:.1f} mm wide", color='lime', ha='center')

    ax1.axvspan(-400, -240, color='red', alpha=0.1, label="Exposed Tube (Dead Zone)")
    ax1.axvspan(240, 400, color='red', alpha=0.1)

    ax1.set_title("Longitudinal Temperature Profile (X=0, Y=0)", fontweight='bold')
    ax1.set_xlabel("Z-Axis Position (mm)")
    ax1.set_ylabel("Temperature (°C)")
    ax1.legend(loc="lower center")
    ax1.grid(color='gray', linestyle='--', alpha=0.3)
    fig1.savefig(graphs_dir / "analysis_longitudinal_profile.png", dpi=300)
    fig1.savefig(graphs_dir / "analysis_longitudinal_profile.png", dpi=300)
    plt.close(fig1)

    print("4. Generating Radial Profile (X-Axis)...")
    x_vals = np.linspace(0, r_max, 500) 
    points_R = np.vstack((x_vals, np.zeros_like(x_vals), np.zeros_like(x_vals)))
    T_rad = probe_mesh(points_R)

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.plot(x_vals * scale_to_mm, T_rad, color='cyan', linewidth=2)
    
    layers = [
        (0, 33, "Air Core", "gray"),
        (33, 33.5, "Quartz", "white"),
        (33.5, 37.5, "Air Gap", "gray"),
        (37.5, 80, "Cement", "orange"),
        (80, 240, "Ceramic Wool", "brown")
    ]
    
    for r_min_layer, r_max_layer, name, col in layers:
        ax2.axvspan(r_min_layer, r_max_layer, color=col, alpha=0.15)
        ax2.text((r_min_layer + r_max_layer)/2, np.nanmin(T_rad) + 50, name, color=col, rotation=90, va='bottom', ha='center')

    ax2.axvline(45.0, color='red', linestyle='--', label="Coiled Coil Location")

    ax2.set_title("Radial Thermal Drop-Off (Z=0)", fontweight='bold')
    ax2.set_xlabel("Distance from Center (mm)")
    ax2.set_ylabel("Temperature (°C)")
    ax2.legend()
    ax2.grid(color='gray', linestyle='--', alpha=0.3)
    fig2.savefig(graphs_dir / "analysis_radial_dropoff.png", dpi=300)
    fig2.savefig(graphs_dir / "analysis_radial_dropoff.png", dpi=300)
    plt.close(fig2)

    print("5. Generating Parallel Core Heatmap...")
    z_para = np.linspace(-z_max*0.375, z_max*0.375, 200) 
    y_para = np.linspace(-r_max*0.25, r_max*0.25, 200) 
    Z_grid, Y_grid = np.meshgrid(z_para, y_para)
    X_fixed = np.zeros_like(Z_grid)

    points_P = np.vstack((X_fixed.flatten(), Y_grid.flatten(), Z_grid.flatten()))
    temp_P = probe_mesh(points_P)

    fig3, ax3 = plt.subplots(figsize=(12, 6))
    hm3 = ax3.pcolormesh(Z_grid * scale_to_mm, Y_grid * scale_to_mm, temp_P.reshape(Z_grid.shape), cmap='inferno', shading='auto')
    cb3 = plt.colorbar(hm3, ax=ax3)
    cb3.set_label('Temperature (°C)', fontweight='bold')
    
    ax3.add_patch(patches.Rectangle((-58, -45), 116, 90, fill=False, edgecolor='cyan', linestyle=':', linewidth=1.5, alpha=0.8))
    ax3.text(0, 48, "Heater Core Boundary", color='cyan', ha='center')

    ax3.set_title("Core Longitudinal Heatmap (X=0)", fontweight='bold')
    ax3.set_xlabel("Z-Axis Length (mm)")
    ax3.set_ylabel("Y-Axis Height (mm)")
    ax3.set_aspect('equal')
    fig3.savefig(graphs_dir / "analysis_heatmap_parallel.png", dpi=300)
    fig3.savefig(graphs_dir / "analysis_heatmap_parallel.png", dpi=300)
    plt.close(fig3)

    print("6. Generating Full Perpendicular Heatmap...")
    x_perp = np.linspace(-r_max, r_max, 200) 
    y_perp = np.linspace(-r_max, r_max, 200) 
    X_grid2, Y_grid2 = np.meshgrid(x_perp, y_perp)
    Z_fixed2 = np.zeros_like(X_grid2)

    points_P2 = np.vstack((X_grid2.flatten(), Y_grid2.flatten(), Z_fixed2.flatten()))
    temp_P2 = probe_mesh(points_P2)

    fig4, ax4 = plt.subplots(figsize=(8, 8))
    hm4 = ax4.pcolormesh(X_grid2 * scale_to_mm, Y_grid2 * scale_to_mm, temp_P2.reshape(X_grid2.shape), cmap='inferno', shading='auto')
    cb4 = plt.colorbar(hm4, ax=ax4, fraction=0.046, pad=0.04)
    cb4.set_label('Temperature (°C)', fontweight='bold')

    ax4.add_patch(patches.Rectangle((-240, -240), 480, 480, fill=False, edgecolor='white', linewidth=2, alpha=0.5))
    ax4.text(-230, 220, "Outer Steel Casing", color='white')

    ax4.set_title("Full Furnace Cross-Section (Z=0)", fontweight='bold')
    ax4.set_xlabel("X (mm)")
    ax4.set_ylabel("Y (mm)")
    ax4.set_aspect('equal')
    fig4.savefig(graphs_dir / "analysis_heatmap_full_perp.png", dpi=300)
    plt.close(fig4)

    print("Analysis Complete! 4 High-Resolution Engineering Graphs saved to data/graphs/")

if __name__ == "__main__":
    analyze_furnace()