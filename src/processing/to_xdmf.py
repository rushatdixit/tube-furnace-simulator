"""
Processing module for converting Gmsh mesh formats to XDMF format suitable for FEniCS.
"""

import meshio
import numpy as np

def extract_mesh_data() -> None:
    """
    Reads a gmsh file, extracts all 3D tetrahedron blocks and physical tags,
    assembles them into continuous arrays, and writes a pristine FEniCS
    compatible mesh in XDMF format.
    """
    print("Reading Gmsh file...")
    import os
    os.makedirs("data/xdmf-s", exist_ok=True)
    msh = meshio.read("data/meshes/flawless_homogenized_furnace.msh")

    print("Extracting ALL 3D Tetrahedron blocks and Physical Tags...")
    tetra_blocks = []
    tetra_data = []

    for i, block in enumerate(msh.cells):
        if block.type == "tetra":
            tetra_blocks.append(block.data)
            tetra_data.append(msh.cell_data["gmsh:physical"][i])

    if not tetra_blocks:
        raise ValueError("No tetrahedrons found! Make sure you meshed 3D.")

    tetra_cells_array = np.vstack(tetra_blocks)
    tetra_data_array = np.concatenate(tetra_data)

    clean_mesh = meshio.Mesh(
        points=msh.points,
        cells=[("tetra", tetra_cells_array)],
        cell_data={"subdomains": [tetra_data_array]}
    )

    print(f"Total elements converted: {len(tetra_cells_array)}")
    print("Writing to XDMF format for FEniCS...")
    meshio.write("data/xdmf-s/furnace_domain.xdmf", clean_mesh)
    print("Success! Created data/xdmf-s/furnace_domain.xdmf and data/xdmf-s/furnace_domain.h5")

if __name__ == "__main__":
    extract_mesh_data()