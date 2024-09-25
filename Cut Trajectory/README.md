# DCD Frame Editor

## Overview

`dcd_frame_editor.py` is a Python script designed to efficiently manipulate and extract frames from Molecular Dynamics (MD) trajectory files in the DCD format. Using the powerful [MDAnalysis](https://www.mdanalysis.org/) library, this tool allows users to easily skip frames, select specific ranges, and generate new DCD files.

This tool is ideal for researchers working in computational chemistry, bioinformatics, and related fields who need to modify large trajectory files for further analysis or visualization.

## Features

- Select specific ranges of frames to process (start, end, step)
- Efficiently skip frames to reduce output file size
- Supports input PDB files for topology and DCD files for trajectory
- Generates new DCD files with selected frames
- Easy to use with command-line arguments

## Requirements

- Python 3.6+
- [MDAnalysis](https://www.mdanalysis.org/) (install using `pip` or `conda`)

### Install MDAnalysis

You can install MDAnalysis using `pip`:

```
pip install MDAnalysis
```

Or using `conda`:

```
conda install -c conda-forge mdanalysis
```

## Usage

```
python dcd_frame_editor.py -p <pdb_file> -d <dcd_file> -o <output_file> [-b <begin>] [-e <end>] [-s <step>]
```

### Command-line Arguments:

- `-p`, `--pdbfile`: Input PDB file for topology (default: md.pdb)
- `-d`, `--dcdfile`: Input DCD file for trajectory (default: md.dcd)
- `-o`, `--outfile`: Output DCD file (default: md_new.dcd)
- `-b`, `--begin`: Starting frame (default: 0)
- `-e`, `--end`: Ending frame (default: the last frame in the trajectory)
- `-s`, `--step`: Step size for skipping frames (default: 1)

## Example

To extract every 10th frame between frame 100 and frame 500 from a trajectory:

```bash
python dcd_frame_editor.py -p topology.pdb -d trajectory.dcd -o filtered.dcd -b 100 -e 500 -s 10
```

This will output a new DCD file (`filtered.dcd`) containing every 10th frame from the specified range.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! If you have suggestions for improvements or bug fixes, feel free to open an issue or submit a pull request.

## Acknowledgements

- MDAnalysis for providing powerful tools for molecular dynamics analysis.
