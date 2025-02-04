#!/usr/bin/env python3
import argparse
from ase import io
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description=("Read two multi-frame xyz files (one with positions and energy and one with forces), "
                     "convert energies and forces from hartree to eV, and write out an extxyz file "
                     "with the energy stored under the key 'energy' and forces stored under the key 'frc'.")
    )
    parser.add_argument("--positions", 
                        help="XYZ file with atomic positions (in Å) and energy (in hartree) stored in the second line with key 'E'")
    parser.add_argument("--forces", 
                        help="XYZ file with forces (in hartree/Å) stored as coordinates")
    parser.add_argument("--output", 
                        help="Output filename (with .xyz extension) in extxyz format")
    args = parser.parse_args()

    # Conversion factor from hartree to eV
    hartree_to_ev = 27.211386024367243

    # Read all frames from the positions and forces files.
    pos_frames = io.read(args.positions, index=":")
    frc_frames = io.read(args.forces, index=":")

    if len(pos_frames) != len(frc_frames):
        raise ValueError(f"Mismatch in number of frames: positions file has {len(pos_frames)} "
                         f"while forces file has {len(frc_frames)}.")

    # Process each frame
    frames = []
    for i, (atoms_pos, atoms_frc) in enumerate(zip(pos_frames, frc_frames)):
        # Get the energy from the positions file.
        # We expect that each frame has an 'E' entry in the info dictionary.
        if 'E' not in atoms_pos.info:
            raise KeyError(f"Frame {i} in positions file does not contain key 'E' in its comment line.")
        # Remove the original key so that we only have the converted value in the output.
        try:
            energy_hartree = float(atoms_pos.info.pop('E'))
        except ValueError:
            raise ValueError(f"Frame {i}: Unable to convert energy '{atoms_pos.info['E']}' to float.")
        energy_ev = energy_hartree * hartree_to_ev

        # Store the converted energy under the key 'energy' (this will appear in the extxyz comment line)
        atoms_pos.info['energy'] = energy_ev

        # In the forces file, the per-atom forces are stored in the positions attribute.
        # Multiply by the same conversion factor (hartree/Å to eV/Å).
        forces_hartree = atoms_frc.positions.copy()  # copy to avoid modifying the original
        forces_ev = forces_hartree * hartree_to_ev

        # Save the converted forces under the key 'frc'. ASE will include array properties
        # in the extxyz comment line.
        atoms_pos.arrays['frc'] = forces_ev

        frames.append(atoms_pos)

    # Write out all frames in extxyz format.
    # The extra info in atoms.info and arrays in atoms.arrays will be included
    # in the second line of each frame.
    io.write(args.output, frames, format="extxyz")
    print(f"Converted {len(frames)} frames written to {args.output}")

if __name__ == '__main__':
    main()

