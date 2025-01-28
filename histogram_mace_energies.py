import argparse
import pandas as pd
import numpy as np
import sys

import matplotlib.pyplot as plt

from ase.io import read, write
from ase.units import fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import Langevin
from ase.md import MDLogger

from mace.calculators import MACECalculator

def parse():
    parser = argparse.ArgumentParser(description="This script will take a xyz trajectory and calculate the MACE energy for each frame,\
                                                  then for a list of nodes it will calculate the node_energy and plot a histogram of the energies\
                                                  for each element.")
    parser.add_argument("-i","--input", help="XYZ Trajectory File", required=True)
    parser.add_argument("-n","--nframes", help="Number of frames to run", type=int, default=None)
    parser.add_argument("-m","--model", help="MACE model file", required=True)
    parser.add_argument("-e","--elements", required=True, type=str, nargs="+", help="Element (s) to print histogram for")
    parser.add_argument("-d","--device", type=str, help="Device to run on", default="cuda")
    return parser.parse_args()

def get_node_energies(xyz,calculator,nframes):
    indices=list(range(len(xyz[0])))*nframes
    elements=xyz[0].get_chemical_symbols()*nframes
    node_energies=None
    for i in range(args.nframes):
        print(f"processing frame {i}")
        xyz[i].calc=calculator
        xyz[i].get_potential_energy()
        if node_energies is None:
            node_energies=xyz[i].calc.results['node_energy']
        else: 
            node_energies=np.append(node_energies,xyz[i].calc.results['node_energy']) 
    
    df=pd.DataFrame()
    df["index"]=indices
    df["element"]=elements
    df["node_energy"]=node_energies
    return df

if __name__=="__main__":
    args=parse()
    calculator=MACECalculator(model_paths=args.model,device=args.device)
    xyz=read(args.input,":")
    if args.nframes is None:
        args.nframes=len(xyz)
    
    node_energies_df=get_node_energies(xyz,calculator,args.nframes)
    
    for element in args.elements:
        element_df=node_energies_df[node_energies_df["element"]==element]
        plt.hist(element_df["node_energy"],bins=100,label=element)
        plt.legend()
        plt.xlabel("Node Energy")
        plt.ylabel("Frequency")
        plt.title("Histogram of Node Energies")
        plt.show()
    
    


