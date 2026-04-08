"""
Main entry point for the furnace thermal simulation pipeline.
This script coordinates mesh generation, format conversion, and FEA solving.
"""

from src.geometry.mesher import run_mesher
from src.processing.to_xdmf import extract_mesh_data
from src.model.simulation import run_advanced_solver

def main() -> None:
    print("1. Generating Furnace Mesh...")
    run_mesher()
    
    print("\n2. Converting Mesh to XDMF...")
    extract_mesh_data()
    
    print("\n3. Running FEA Solver & Exporting...")
    run_advanced_solver()

    print("\n Simulation complete. Please check furnace/data.")

if __name__ == "__main__":
    main()
