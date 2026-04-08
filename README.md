# Furnace Thermal Simulation

A Finite Element Analysis (FEA) based thermal simulation of an electric heating furnace. This project models the spatial buoyancy gradients, heat conduction, and convective heat transfer of a tubular furnace apparatus. Built on FEniCSx and Gmsh.

## Prerequisites

To run this project, we rely on advanced FEA scientific computing libraries which are best installed via `conda-forge`. FEniCSx and Gmsh have complex C++ bindings that Conda handles automatically.

### System Setup
If you don't already have Conda installed, you can quickly install `Miniforge` (a lightweight Conda installer prioritizing `conda-forge`) entirely from your terminal:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p $HOME/miniforge3

source $HOME/miniforge3/etc/profile.d/conda.sh
conda init
```
*(Note: You may need to restart your terminal or open a new window afterward for the `conda` command to be fully active.)*

## Installation Instructions

1. **Clone/Download the Repository**
Navigate to the desired location on your machine and ensure you are in the project root:
```bash
cd /path/to/furnace
```

2. **Create a Conda Environment**
We will create an isolated Conda environment and install Python.
```bash
conda create -n furnace_env -c conda-forge python=3.11
conda activate furnace_env
```

3. **Install Dependencies**
Using the provided `requirements.txt`, install all required mathematical and FEA packages simultaneously. By using `conda-forge`, it ensures `dolfinx`, `mpi4py`, and `gmsh` link properly.
```bash
conda install -c conda-forge --file requirements.txt
```

*(Note: Depending on your exact architecture (e.g., Apple Silicon M-series vs Intel), Conda will automatically fetch the right compatible MPICH/OpenMPI libraries.)*

---

## Running the Simulation Pipeline

The entire system is orchestrated by a single script that runs from the **root directory** (`/furnace`). It sequentially generates the mesh, converts the geometry into FEA data, and executes the thermal solver.

```bash
python simulate.py
```

## Outputs

The simulation script will automatically construct any missing subdirectories inside the `data` folder. You can expect:
- **`data/visuals/`**: High-quality 2D annotated heatmaps of both longitudinal and perpendicular cross-sections.
- **`data/xdmf-s/`**: The complete 3D simulation results capable of being opened natively in applications like **ParaView**.
- **`data/meshes/`**: The compiled Gmsh tetrahedron files.

---

## Exploring with ParaView

Because this project solves complex 3D thermodynamics, the primary raw outputs are `.xdmf` matrices. To view these interactively:

### 1. Install ParaView
You can easily install ParaView directly from the terminal (if you are on macOS with Homebrew):
```bash
brew install --cask paraview
```
*(On Linux/Windows, or without Homebrew, download it directly from [paraview.org](https://www.paraview.org/download/)).*

### 2. Loading the Model
1. Open the **ParaView** application.
2. Go to **File > Open** and locate your `data/xdmf-s/furnace_temperature.xdmf` file.
3. In the left-hand **Properties** panel, click the green **Apply** button to render the physical domain block.

### 3. Rendering Features
- **Color by Temperature**: In the top toolbar, find the dropdown menu that typically defaults to `Solid Color` and switch it to **`Temperature`**. The heat gradient will reflect our FEA solutions immediately.
- **Cross-section Slicing**: To look inside at the quartz tube and heating zone, click the **Slice** tool (icon of a plane cutting a cube) in the top menu. In its properties panel, choose the normal axis you wish to slice across (e.g. `X Normal`), uncheck the 'Show Plane' visualizer, and hit **Apply**.
