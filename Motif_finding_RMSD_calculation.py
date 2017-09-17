import Bio.PDB as PDB
import os
import numpy as np
import matplotlib.pyplot as plt

## FILES PARSING AND MATCH FINDING ##

# Create parser object that does not take in account possible warnings
parser = PDB.PDBParser(QUIET = True)

# Create polypeptide builder object that takes in account the inter-CA distance
ppb = PDB.CaPPBuilder(radius = 4)

# Dictionary with IDs as keys and sequences as values
seq_dict = {}

def CA_coords(pentapeptide):
    """This function receives a pentapeptide object and returns the CA coordinates of its residues"""
    for residue in pentapeptide:
        return residue['CA'].get_coord()

# Work over all pdb. files of the top100H directory
for pdb_file in os.listdir('top100H'):

    # Try to perform this operation with all the pdb. files, leaving the damaged ones out (only 1 detected)
    try:
        # Get structure from the parsed file
        structure = parser.get_structure(pdb_file, 'top100H/'+str(pdb_file))

        # Build polypeptide from structure checking if residues are clean (not water or missing CA)
        for poly in ppb.build_peptides((structure), aa_only = 1):
            for i in range(len(poly) - 4):

                # For each pentapeptide: if its sequence is among the dictionary keys
                if str(poly.get_sequence()[i:i + 5]) in seq_dict:

                    # add its CA coordinates as value (by calling CA_coords function)
                    seq_dict[str(poly.get_sequence()[i:i + 5])].append([CA_coords(poly[i:i + 5])])

                else:
                    # If it is not, add it as a key with its CA coordinates
                    seq_dict[str(poly.get_sequence()[i:i + 5])] = [CA_coords(poly[i:i + 5])]

    except:
        pass

## RMSD CALCULATION ##

RMSD_list = []

## 398 matching pentapeptides found ##

for val_pair in seq_dict.values():
    # Discard the pentapeptides that appear more than twice
    if len(val_pair) == 2:

        #Set numpy format for each pair, define X and Y, get transpose of X and correlation matrix (R)
        X = np.asarray(val_pair[0])
        Y = np.asarray(val_pair[1])
        X_trans = X.T
        R = Y * X_trans

        # Perform Single Value Decomposition and return V, S and the transpose of W
        SVD = np.linalg.svd(R)
        V = np.asarray(SVD[1])
        W_trans = np.asarray(SVD[2])
        Z = np.array([1, 1, -1])
        U = W_trans.T * V.T # W * V_trans

        # Apply RMSD
        RMSD = np.sqrt(0.5 * np.sum(np.abs(np.power((X - U*Y), 2))))
        RMSD_list.append(RMSD)

## HISTOGRAM CONSTRUCTION ##

x = np.array(RMSD_list)
plt.hist(x, 50, facecolor = 'green', align = 'mid')
plt.axis([0, 200000, 0, 220])
plt.title('Histogram of RMSD values for equal-sequence pentapeptides')
plt.xlabel('RMSD values')
plt.ylabel('Probability')
plt.grid(True)

plt.show()