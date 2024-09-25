#!/usr/bin/env python
# Frame editor for DCD file
# Inputs: PDB file, trajectory DCD file
# Outputs: target DCD file

import sys
import MDAnalysis as mda
import argparse

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        description='DCD frame skipper')
    
    # Argument definitions
    parser.add_argument('-p', '--pdbfile', default='md.pdb', 
                        help='Input PDB file for topology')
    parser.add_argument('-d', '--dcdfile', default='md.dcd', 
                        help='Input DCD file')
    parser.add_argument('-o', '--outfile', default='md_new.dcd', 
                        help='Output DCD file')
    parser.add_argument('-b', '--begin', type=int, default=0, 
                        help='Starting frame (after equilibration)')
    parser.add_argument('-e', '--end', type=int, 
                        help='Ending frame')
    parser.add_argument('-s', '--step', type=int, default=1, 
                        help='Step size between frames')
    args = parser.parse_args()

    # Load universe
    u = mda.Universe(args.pdbfile, args.dcdfile)
    real = u.select_atoms("all")
    
    # Define the frame range
    end_frame = args.end if args.end is not None else len(u.trajectory)
    begin_frame = max(0, args.begin)
    
    if begin_frame >= end_frame:
        print("Error: begin frame is greater than or equal to end frame.")
        sys.exit(1)

    # Open writer
    with mda.Writer(args.outfile, len(real)) as ww:
        print(f"Processing {len(u.trajectory)} total frames...")
        frame_count = 0
        
        # Loop through the trajectory and write frames
        for ts in u.trajectory[begin_frame:end_frame:args.step]:
            ww.write(real)
            frame_count += 1
        
        print(f"New file: wrote {frame_count} frames.")

if __name__ == '__main__':
    main()
