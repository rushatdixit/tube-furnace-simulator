import gmsh
import math
import sys
from pathlib import Path

def generate_embedded_1d_coil():
    print("1. Initializing Gmsh...")
    gmsh.initialize(sys.argv)
    gmsh.model.add("coiled_coil_1d")
    occ = gmsh.model.occ

    N_major = 13                 
    minor_turns_per_major = 50   
    R_major = 45.0               
    R_minor = 5.0                
    pitch_major = 8.9            

    total_minor_turns = N_major * minor_turns_per_major
    total_length = N_major * pitch_major

    print("2. Computing 1D Double Helix Skeleton (40 pts/turn)...")
    points = []
    points_per_minor_turn = 40  # Your high-fidelity upgrade!
    total_points = total_minor_turns * points_per_minor_turn

    for i in range(total_points + 1):
        t = (i / total_points) * (2 * math.pi * N_major)
        phi = t * minor_turns_per_major

        x = (R_major + R_minor * math.cos(phi)) * math.cos(t)
        y = (R_major + R_minor * math.cos(phi)) * math.sin(t)
        z = (pitch_major * (t / (2 * math.pi))) - (total_length / 2) + R_minor * math.sin(phi)
        points.append(occ.addPoint(x, y, z))

    spline = occ.addSpline(points)
    wire = occ.addWire([spline])

    print("3. Generating 3D Cement Block...")
    cement_block = occ.addCylinder(0, 0, -(total_length/2 + 10), 0, 0, total_length + 20, R_major + R_minor + 10)

    print("4. Embedding 1D Wire into 3D Block (Fragmentation)...")
    occ.fragment([(3, cement_block)], [(1, wire)])
    occ.synchronize()

    print("5. Assigning Physical Groups for FEniCSx...")
    # Tag 1: The 3D Cement Volume
    cement_vols = [v[1] for v in gmsh.model.getEntities(3)]
    gmsh.model.addPhysicalGroup(3, cement_vols, 1)
    gmsh.model.setPhysicalName(3, 1, "Cement_Block")

    # Tag 2: The 1D Wire (We grab all 1D curves)
    wire_curves = [c[1] for c in gmsh.model.getEntities(1)]
    gmsh.model.addPhysicalGroup(1, wire_curves, 2)
    gmsh.model.setPhysicalName(1, 2, "Coiled_Coil_Wire")

    print("6. Configuring Smart Mesh Density...")
    gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 15.0)
    
    curves = gmsh.model.getEntities(1)
    gmsh.model.mesh.setSize(gmsh.model.getBoundary(curves, False, False, True), 1.5)

    print("7. Generating Mesh...")
    gmsh.model.mesh.generate(3)

    print("8. Saving Tagged Mesh...")
    base_dir = Path.cwd() / "data" / "meshes"
    base_dir.mkdir(parents=True, exist_ok=True)
    
    output_path = base_dir / "embedded_1d_coiled_coil.msh"
    gmsh.write(str(output_path))
    
    print(f"Success! Tagged 1D embedded mesh saved to {output_path}")
    #gmsh.fltk.run()
    gmsh.finalize()

if __name__ == "__main__":
    generate_embedded_1d_coil()