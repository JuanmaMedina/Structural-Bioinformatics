import Bio.PDB as PDB

## 1. Download the structure with PDB identifier 2PTC from the PDB database. ##
pdbl = PDB.PDBList()
pdbl.retrieve_pdb_file('2PTC')

## 2. Print out the 3D coordinates of the C-alpha atom of the amino acid with PDB residue number 20 in chain E. ##
# Create parser object
p = PDB.PDBParser()
# Get structure
s = p.get_structure('2PTC', '/home/juanma/PycharmProjects/Bioinformatics/pt/pdb2ptc.ent')

# Extract CA atom from structure and print its 3D coordinates
CA_20_E = s[0]['E'][20]['CA']
print "Coordinates of the CA of the aminoacid number 20 in chain E:"
print CA_20_E.get_coord()

## 3. Calculate the geometric center (centroid) of the above amino acid and print it out. ##

# Extract residue number 20
residue_20 = s[0]['E'][20]

# Bio.PDB.Vector to calculate centroid
atom_sum_20 = PDB.Vector(0., 0., 0.)

# Atoms of the residue number 20 counter
atoms_20 = 0

for atom in residue_20:
    # In-place element-wise addition
    atom_sum_20 += atom.get_vector().get_array()
    atoms_20 += 1

# Divide by the number of atoms to get centroid and print its 3D coordinates
centroid_20 = atom_sum_20.get_array() / atoms_20
print "Coordinates of the centroid of the aminoacid number 20 in chain E:"
print centroid_20

## 4. Create a function "centroid" that takes a residue as a parameter and returns the coordinates of its centroid. ##

def centroid(residue):
    """ Returns 3D coordinates of the centroid of a residue"""
    atom_sum = PDB.Vector(0., 0., 0.,)
    atoms = 0

    for atom in residue:
        atom_sum += atom.get_vector().get_array()
        atoms += 1

    centroid_coordinates = atom_sum.get_array() / atoms

    return centroid_coordinates

# Test of the centroid function
# print centroid(s[0]['E'][20])

## 5. Print out the 3D coordinates of the C-alpha atom of the amino acid with PDB residue number 13 in chain I. ##

CA_13_I = s[0]['I'][13]['CA']
print "Coordinates of the CA of the aminoacid number 13 in chain I:"
print CA_13_I.get_coord()

## 6. Calculate the centroid of the above amino acid and print it out. ##

print "Coordinates of the centroid of the aminoacid number 13 in chain I:"

# Call to the centroid function defined above
print centroid(s[0]['I'][13])