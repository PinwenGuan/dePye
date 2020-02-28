# dePye
A python code for EOS fitting and calculations of thermodynamics properties based on the Debye model with DFT uncertainty quantification functions.

## Install dePye

Suppose you download the files and put them in the folder \~/bin. Then add the related paths by adding in your .bashrc:<br>
```
export PATH=~/bin/depye:$PATH
export PYTHONPATH=~/bin/depye:$PYTHONPATH
```
The EOS fitting is dependent on pymatgen, which should be installed first.<br>

## Prepare the input files

Since you have finished the energy calculations, you should already have the structure file in your working folder (if not you need to do so). Currently dePye supports the structure files of VASP and Quantum ESPRESSO (QE), i.e., POSCAR and the .in file.<br> 
Normally, the only input file you need to prepare is e-v.dat (can be a name ended with .dat) such as:<br>
```
# V(angstrom^3) E(eV) vol_0.98 vol_0.99 vol_1.00 vol_1.01 vol_1.02<br>

69.131955   -7728.27656469   -7728.27686878   …
71.217852   -7728.28518104   -7728.29116264   …
73.345300   -7728.25348998   -7728.25450442   …
75.514684   -7728.18594466   -7728.17339096   …
77.726434   -7728.08681822   -7728.05316957   …
```
The comment lines do not matter. The first column is the volume (in angstrom^3), and the second, third, … is the energies (in eV). At least 5 V-E data points are needed. In the case of multiple sets of energies, each energy set will generate a result. For example, one can get two E-V datasets with and without relaxation, and get two sets of properties.<br>
Note: the input data should be corresponded to the number of atoms in the structure file. 

## Run dePye

First activate the virtual environment where pymatgen is installed:<br>
```
source activate (pymatgen environment)
```
In your working folder, type:<br>
```
depye -tmax=1000 –tstep=10 POSCAR
```
Here,  tmax means the maximum temperature you want to calculate, tstep is the temperature step in unit of 10 K in the output file, and POSCAR is the structure file. If not specified, the default setting will be used, i.e., tmax=1000, tstep=1 and POSCAR. The order of these flags does not matter. Another example:<br>
```
depye -tmax=1000 cu.in
```
where the default tstep=1 (10 K) and the QE input file cu.in will be used.<br>
You can also use an input file whose name is not e-v.dat (but should end with .dat):<br>
```
depye POSCAR e-v2.dat
```
You can adjust the scaling factor in the Debye temperature (default is 0.617):<br>
```
depye POSCAR s=0.8
```
You can put a 'poisson' file which only contains the value of the poisson ratio of the material, then the code will calculate the scaling factor from it and override the default or the assigned s.<br>
You can turn off showing figures by adding show='F':<br>
```
depye POSCAR show='F'
```
You can determine how many atoms the output quantities are corresponded to:<br>
```
depye POSCAR nat=1
```
The above outputs the quantities per atom.<br>
Note: the default setting is per formula.<br>
You can decide which properties are kept in the outputed figures:<br>
```
depye POSCAR prop=[0,3]
```
The above shows thermodynamic quantities for only Gibbs energy and entropy. The quantity codes are<br>
```
G V B S H TEC Cp Cv Bp
0 1 2 3 4  5  6  7  8
```
You can append experimental data (or data you get by other methods, e.g., phonon method):<br>
```
depye POSCAR data.expt
```
Note: the default name of the experimental file is expt. If you want to use another name, the file should end with .expt.<br> 
Note: for VASP, the structure file should be named “POSCAR”. For QE, it should be ended with “.in”. The expt file should look like:<br>
```
# T (K) G (kJ/mol) V (cm^3/mol) B (GPa) S (J/mol/K) H (kJ/mol) TEC (1e-6/K) Cp (J/mol/K) Cv (J/mol/K)
T B ref
298 100 1
298 95 2 
500  90 2

T TEC ref
298 40 3

#You can name the references for literature data shown in the figures by alias below. The names are in the order as 1,2,3
#alias: author A, author B, author C
```

## Get the results

After dePye running, a series of figures about EOS and thermodynamic properties will be prompted out.  The output file “Debye-thermo” contains the EOS parameters and the data of Gibbs energy, volume, bulk modulus, entropy, enthalpy, thermal expansion and heat capacity as functions of temperature, which can be further as inputs for thermodynamic modelling.<br> 
