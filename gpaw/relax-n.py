from readst import readst
import os
import numpy as np
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter
from gpaw import * #imports evertyhing from gpaw
from ase.dft.bee import BEEFEnsemble
from ase.parallel import parprint
from ase.io import read,Trajectory
from ase.io import write


atoms = read('POSCAR')
sym = atoms.get_chemical_symbols()

#get information from previous calculations

#-----------------kp2,nat
nat=len(atoms.get_atomic_numbers())
kp2=[1 for i in range(3)]
for i in range(3):
    kp2[i]=int(np.ceil(20/np.linalg.norm(atoms.cell[i,:])))

vdW='BEEF-vdW'

calc = GPAW(xc=vdW,mode=PW(700, dedecut='estimate'),kpts={'size': kp2, 'gamma': True},occupations=FermiDirac(width=0.05))

atoms.set_calculator(calc)
name = atoms.get_chemical_formula(mode='hill')

atoms.calc.set(txt=name+'_final.txt')
atoms.calc.attach(atoms.calc.write, 5, name+'_final.gpw')
uf = ExpCellFilter(atoms,constant_volume=True)
relax = BFGS(uf)
traj = Trajectory(name+'_final.traj', 'w', atoms)
relax.attach(traj)
relax.run(fmax=0.03)  
atoms.calc.write(name+'_final.gpw')
#Writing ensemble energies to file
if  atoms.calc.get_xc_functional()=='BEEF-vdW':
    ens = BEEFEnsemble(atoms.calc)
    ens_material = ens.get_ensemble_energies()
    np.savetxt(str(name)+'_Ensemble_Energies.txt',ens_material)
