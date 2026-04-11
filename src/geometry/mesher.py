"""
Geometry and meshing module for the furnace using Gmsh.
Generates a 3D geometry of a tagged furnace with different physical volumes, 
fragments them, and outputs a .msh file.
"""

import gmsh
import os

def run_mesher() -> None:
    """
    Executes the gmsh geometry creation and meshing routine.
    Outputs the file to data/meshes/flawless_homogenized_furnace.msh
    """
    os.makedirs("data/meshes", exist_ok=True)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("Tagged_Furnace")
    occ = gmsh.model.occ

    R_AIR = 33e-3
    R_QUARTZ = 35e-3
    R_GAP = 37.5e-3
    R_INNER_CEM = 40e-3
    R_HEATING = 58e-3
    R_OUTER_CEM = 80e-3
    BOX_HALF = 240e-3
    STEEL_THICK = 2e-3
    Z_PLANES = [-400e-3, -290e-3, -240e-3, -225e-3, 225e-3, 240e-3, 290e-3, 400e-3]

    def make_annulus(z_start, z_end, r_in, r_out):
        """
        Creates an annular solid cylinder between z_start and z_end with inner radius r_in and outer radius r_out.
        """
        length = z_end - z_start
        solid = occ.addCylinder(0, 0, z_start, 0, 0, length, r_out)
        if r_in > 0:
            void = occ.addCylinder(0, 0, z_start, 0, 0, length, r_in)
            return occ.cut([(3, solid)], [(3, void)], removeTool=True)[0][0][1]
        return solid

    tags = {
        "air": [], "quartz": [], "plugs": [], "gap": [], 
        "cement": [], "heating_zone": [], "wool": [], "steel": []
    }

    print("1. Building and Categorizing Slices...")

    for z in [(Z_PLANES[0], Z_PLANES[1]), (Z_PLANES[6], Z_PLANES[7])]:
        tags["air"].append(make_annulus(z[0], z[1], 0, R_AIR))
        tags["quartz"].append(make_annulus(z[0], z[1], R_AIR, R_QUARTZ))

    for z in [(Z_PLANES[1], Z_PLANES[2]), (Z_PLANES[5], Z_PLANES[6])]:
        tags["plugs"].append(make_annulus(z[0], z[1], 0, R_AIR))
        tags["quartz"].append(make_annulus(z[0], z[1], R_AIR, R_QUARTZ))

    for z in [(Z_PLANES[2], Z_PLANES[3]), (Z_PLANES[4], Z_PLANES[5])]:
        tags["air"].append(make_annulus(z[0], z[1], 0, R_AIR))
        tags["quartz"].append(make_annulus(z[0], z[1], R_AIR, R_QUARTZ))
        tags["gap"].append(make_annulus(z[0], z[1], R_QUARTZ, R_GAP))
        tags["cement"].append(make_annulus(z[0], z[1], R_GAP, R_OUTER_CEM))

    z = (Z_PLANES[3], Z_PLANES[4])
    tags["air"].append(make_annulus(z[0], z[1], 0, R_AIR))
    tags["quartz"].append(make_annulus(z[0], z[1], R_AIR, R_QUARTZ))
    tags["gap"].append(make_annulus(z[0], z[1], R_QUARTZ, R_GAP))
    tags["cement"].append(make_annulus(z[0], z[1], R_GAP, R_INNER_CEM))
    tags["heating_zone"].append(make_annulus(z[0], z[1], R_INNER_CEM, R_HEATING))
    tags["cement"].append(make_annulus(z[0], z[1], R_HEATING, R_OUTER_CEM))

    for z in [(Z_PLANES[2], Z_PLANES[3]), (Z_PLANES[3], Z_PLANES[4]), (Z_PLANES[4], Z_PLANES[5])]:
        L = z[1] - z[0]
        w_s = occ.addBox(-BOX_HALF, -BOX_HALF, z[0], BOX_HALF*2, BOX_HALF*2, L)
        w_v = occ.addCylinder(0, 0, z[0], 0, 0, L, R_OUTER_CEM)
        tags["wool"].append(occ.cut([(3, w_s)], [(3, w_v)])[0][0][1])
        
        sh, ss = BOX_HALF + STEEL_THICK, (BOX_HALF + STEEL_THICK)*2
        s_s = occ.addBox(-sh, -sh, z[0], ss, ss, L)
        s_v = occ.addBox(-BOX_HALF, -BOX_HALF, z[0], BOX_HALF*2, BOX_HALF*2, L)
        tags["steel"].append(occ.cut([(3, s_s)], [(3, s_v)])[0][0][1])

    print("2. Fragmenting and Tracking Topology...")
    all_vols_tuples = []
    vol_to_name = {}
    for name, vol_list in tags.items():
        for v in vol_list:
            all_vols_tuples.append((3, v))
            vol_to_name[v] = name

    out, out_map = occ.fragment(all_vols_tuples, [])
    occ.removeAllDuplicates()
    occ.synchronize()

    new_tags = {name: [] for name in tags.keys()}
    for old_tuple, new_entities in zip(all_vols_tuples, out_map):
        old_tag = old_tuple[1]
        name = vol_to_name[old_tag]
        for entity in new_entities:
            new_tags[name].append(entity[1])

    new_tags = {k: list(set(v)) for k, v in new_tags.items()}

    labels = {"air": 1, "quartz": 2, "plugs": 3, "gap": 4, 
              "cement": 5, "heating_zone": 6, "wool": 7, "steel": 8}

    for name, label in labels.items():
        gmsh.model.addPhysicalGroup(3, new_tags[name], label, name)
        print(f"Tagged {name} as Physical Group {label} (Volumes: {new_tags[name]})")

    print("3. Meshing...")
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.003)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.02)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1) 
    gmsh.model.mesh.generate(3)
    gmsh.fltk.run() # Disable GUI for automated pipeline 
    gmsh.write("data/meshes/flawless_homogenized_furnace.msh")
    gmsh.finalize()

if __name__ == "__main__":
    run_mesher()