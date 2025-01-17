import mdtraj as md
from tqdm import tqdm

topology_file = "cht_box.pdb"
trajectory_file = "md_new.dcd"
output_file = "noHOH.dcd"
chunk_size = 1

topology = md.load(topology_file).topology
non_water_indices = topology.select("not water")

with md.formats.DCDTrajectoryFile(output_file, "w") as outfile:
    frame_count = 0
    for chunk in tqdm(md.iterload(trajectory_file, top=topology_file, chunk=chunk_size), desc="Processing chunks"):
        filtered_positions = chunk.xyz[:, non_water_indices, :]
        outfile.write(filtered_positions)
        frame_count += chunk.n_frames

print(f"\nFiltered trajectory saved as '{output_file}'")
print(f"Processed {frame_count} frames.")
