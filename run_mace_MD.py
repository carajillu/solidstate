import argparse

from ase.io import read, write
from ase.units import fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import Langevin
from ase.md import MDLogger

from mace.calculators import MACECalculator

def parse():
    parser = argparse.ArgumentParser(description="This script will take an xyz file and run a molecular dynamics simulation using the MACE calculator.")
    parser.add_argument("-i","--input", help="Initial Structure File", required=True)
    parser.add_argument("-m","--model", help="MACE model file", required=True)
    parser.add_argument("-d","--device", help="Device to run on", default="cuda")
    parser.add_argument("-t","--temperature", type=float, help="Temperature in K", default=300)
    parser.add_argument("-s","--stride",type=int,help="Stride to print output",default=1000)
    parser.add_argument("-n","--nsteps",type=int,help="Number of steps to run",default=1000000)
    return parser.parse_args()

def return_print_snapshot(atoms):
    def print_snapshot():
        #atoms.arrays["forces"]=atoms.get_forces()
        atoms.arrays['node_energy']=atoms.calc.results['node_energy']
        basename=args.input.split(".xyz")[0]
        write(f"{basename}_trj.xyz",atoms,format="extxyz",append=True)
    return print_snapshot

if __name__=="__main__":
    args=parse()

    xyz = read(args.input)
    #xyz.set_pbc([True,True,True])
    calc = MACECalculator(model_paths=args.model,device=args.device)
    xyz.set_calculator(calc)
    MaxwellBoltzmannDistribution(atoms=xyz, temperature_K=args.temperature)
    dyn=Langevin(atoms=xyz, timestep=1*fs, temperature_K=args.temperature, friction=0.1/fs)
    logger=MDLogger(dyn,atoms=xyz,logfile=f"{args.input}.log",header=True)
    dyn.attach(return_print_snapshot(xyz),interval=args.stride)
    dyn.attach(logger, interval=args.stride)
    dyn.run(args.nsteps)