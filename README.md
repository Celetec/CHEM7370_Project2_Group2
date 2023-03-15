# CHEM7370 Spring2023 Project 2 Group 2
import os
from pathlib import Path
import numpy as np
import scipy
Downloads2 = str(os.path.join(Path.home(), 'downloads','reactions.txt'))
reactions = np.genfromtxt(fname=Downloads2, delimiter='/n', dtype='unicode')
print(reactions)

#--IMPORTED DATA FROM FILE--#

def stoichiometric_coefficient(molecule):
    position = 0
    coefficient = ''
    while position < len(molecule):
        if molecule[position].isdigit():
            coefficient += molecule[position]
        else:
            break
        position += 1
    if coefficient == '':
        just_molecule = molecule
        return (1,just_molecule)
    else:
        just_molecule = molecule[len(coefficient):]
        return (int(coefficient),just_molecule)
def count_atoms(molecule):
    charge = 0
    poscharge = molecule.find('+')
    negcharge = molecule.find('-')
    if poscharge != -1:
        charge = molecule.split('+')[1]
        if charge == '':
            charge = 1
        else:
#--DEFINED IMPORTANT VARIABLES --#
#--BEGIN MAJOR FOR LOOP--#

for x in reactions:
    all_reactants = []
    all_products = []
    reactants_together = x.split(" -> ")[0]
    products_together = x.split(" -> ")[1]
    reactants = reactants_together.split(" + ")
    products = products_together.split(" + ")

    #--SPLITS EACH REACTION INPUT--#
    
    allreactants = {}
    allproducts = {}
    for reactant in reactants:
        (n,mol) = stoichiometric_coefficient(reactant)
        atoms = count_atoms(mol)
        for atom in atoms.keys():
            if atom in allreactants.keys():
                allreactants[atom] += n * atoms[atom]
            else:
                allreactants[atom] = n * atoms[atom]
    for product in products:
        (n,mol) = stoichiometric_coefficient(product)
        atoms = count_atoms(mol)
        for atom in atoms.keys():
            if atom in allproducts.keys():
                allproducts[atom] += n * atoms[atom]
            else:
                allproducts[atom] = n * atoms[atom]
    atoms = sorted(allreactants.keys())
    nrows = len(atoms)
    print ('Number of rows:',nrows)
    ncolumns = len(reactants) + len(products) - 1
    print ('Number of columns:',ncolumns)
    mat = np.zeros((nrows,ncolumns))
    column = 0
    for compound in reactants[1:]:
        (n,mol) = stoichiometric_coefficient(compound)
        atoms_here = count_atoms(mol)
        for atomtype in range(len(atoms)):
            if atoms[atomtype] in atoms_here.keys():
                mat[atomtype,column] = -atoms_here[atoms[atomtype]]
                #print(mol,atomtype,atoms[atomtype],atoms_here[atoms[atomtype]])
        column += 1
    for compound in products:
        (n,mol) = stoichiometric_coefficient(compound)
        atoms_here = count_atoms(mol)
        for atomtype in range(len(atoms)):
            if atoms[atomtype] in atoms_here.keys():
                mat[atomtype,column] = atoms_here[atoms[atomtype]]
                #print(mol,atomtype,atoms[atomtype],atoms_here[atoms[atomtype]])
        column += 1
    print(mat)
    rhs = np.zeros(nrows)
    (n,mol) = stoichiometric_coefficient(reactants[0])
    atoms_here = count_atoms(mol)
    for atomtype in range(len(atoms)):
        if atoms[atomtype] in atoms_here.keys():
            rhs[atomtype] = atoms_here[atoms[atomtype]]
    #print(rhs)
    res = scipy.linalg.lstsq(mat, rhs)
    print(res)
            charge = int(charge)
        molecule = molecule.split('+')[0]
    if negcharge != -1:
        charge = molecule.split('-')[1]
        if charge == '':
            charge = -1
        else:
            charge = -int(charge)
        molecule = molecule.split('-')[0]
    position = 0
    atoms = {}
    if charge != 0:
        atoms['e'] = charge
    current_atom = ''
    current_number = ''
    for x in molecule:
        if x.isupper():
            if current_atom != '':
                if current_number == '': 
                    if current_atom in atoms.keys():
                        atoms[current_atom] += 1
                    else:
                        atoms[current_atom] = 1
                else:
                    if current_atom in atoms.keys():
                        atoms[current_atom] += int(current_number)
                    else:
                        atoms[current_atom] = int(current_number)
                current_atom = ''
                current_number = ''
            current_atom += x
        elif x.islower():
            current_atom += x
        else:
            current_number += x
    if current_atom != '':
        if current_number == '':
            if current_atom in atoms.keys():
                atoms[current_atom] += 1
            else:
                atoms[current_atom] = 1
        else:
            if current_atom in atoms.keys():
                atoms[current_atom] += int(current_number)
            else:
                atoms[current_atom] = int(current_number)
    return atoms
    
