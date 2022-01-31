# Molcular-Dynamics-Simulation
## Introduction
Molecular dynamics (MD) simulations predict how every atom in a protein or other molecular system will move over time, based on a general model of the physics governing interatomic interactions. These simulations can capture a wide variety of important biomolecular processes, including conformational change, ligand binding, and protein folding, revealing the positions of all the atoms at femtosecond temporal resolution. Importantly, such simulations can also predict how biomolecules will respond—at an atomic level—to perturbations such as mutation, phosphorylation, protonation, or the addition or removal of a ligand.
## Running MD by Gromacs 
1. We must download the protein structure file with which we will be working. For this tutorial, we will utilize hen egg white lysozyme (PDB code 1AKI). Go to the RCSB website and download the PDB text for the crystal structure.

Once you have downloaded the structure, you can visualize the structure using a viewing program such as VMD, Chimera, PyMOL, etc. Once you've had a look at the molecule, you are going to want to strip out the crystal waters. To delete the water molecules (residue "HOH" in the PDB file), either use a plain text editor like vi, emacs (Linux/Mac), or Notepad (Windows). Do not use word processing software! Alternatively, you can use grep to delete these lines very easily:
```rupy
grep -v HOH complex2.pdb > complex2_clean.pdb
gmx pdb2gmx -ignh -f complex2-clean.pdb -o complex2_processed.gro -water spce
select15, total charge 6
gmx editconf -f complex2_processed.gro -o complex2_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp complex2_newbox.gro -cs spc216.gro -o complex2_solv.gro -p topol.top
gmx grompp -f ions.mdp -c complex2_solv.gro -p topol.top -o ions.tpr #need ions.mdp file
gmx genion -s ions.tpr -o complex2_solv_ions.gro -p topol.top -nname CL -nn 6 #to neutralize +6 charge of system
# for neutralize negetive charge:
gmx genion -s ions.tpr -o target_solv_ions.gro -p target.top -pname NA -np 1 
gmx grompp -f minim.mdp -c complex2_solv_ions.gro -p topol.top -o em.tpr  #need minim.mdp file
gmx mdrun -v -deffnm em &   # htop -u "ghaedi"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1

#For MD:
gmx grompp -f md.mdp -c target-npt.gro -t target-npt.cpt -p target.top -o target-md.tpr
gmx mdrun -v -nt 20 -deffnm target-md
```
for REMD:

# calculate T distribution via webpage then replace ref_t with resulted value from webpage. Now you have N remd_N.mdp files which are same but they have different ref_t. such as; remd_0.mdp , remd_1.mdp , remd_N.mdp

for each remd_N.mdp file, run this script (you should make n(7) remd_.tpr file same as your remd_.mdp files 0, 1, ..., 6): 
grompp -f remd_n.mdp -c target-npt.gro -p target.top -o remd_n.tpr -maxwarn 1
or just use :for i in {0..5}; do grompp -f remd_${i}.mdp -c mix10-npt.gro -p mix10.top -o remd_${i}.tpr -maxwarn 1 ;done

mpirun -np 7 mdrun_mpi -s remd_${i}.tpr -multi 7 -replex 1000 -reseed -1

*There are two very important factors to evaluate to determine if EM was successful.
The first is the potential energy (printed at the end of the EM process, even without -v).
Epot should be negative, and (for a simple protein in water) on the order of 105-106, 
depending on the system size and number of water molecules. The second important feature is the maximum force, Fmax,
the target for which was set in minim.mdp - "emtol = 1000.0" - indicating a target Fmax of no greater than 1000 kJ mol-1 nm-1.
It is possible to arrive at a reasonable Epot with Fmax > emtol. If this happens, your system may not be stable enough for simulation.
Evaluate why it may be happening, and perhaps change your minimization parameters (integrator, emstep, etc).

## Analysis 
