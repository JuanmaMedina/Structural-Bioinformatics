import Bio.PDB as PDB
import numpy as np
import math

## MOVE THE I CHAIN ##

## 0. Download the structure with PDB identifier 2PTC from the PDB database. ##
pdbl = PDB.PDBList()
pdbl.retrieve_pdb_file('2PTC')

## 1. Make two numpy arrays with all the atoms of the residues (if they have CA) belonging chains E and I. ##
# Create parser object
p = PDB.PDBParser()
# Get structure
s = p.get_structure('2PTC', '/home/juanma/PycharmProjects/Bioinformatics/pt/pdb2ptc.ent')

# Extract both chains
chain_E = s[0]['E']
chain_I = s[0]['I']

# Atoms of the CA of both chains
atoms_E = []
atoms_I = []

for res in chain_E:
    # Check if res is sane (not water or missing CA)
    if res.has_id("CA"):
        for atom in res:
            # In-place element-wise addition
            atoms_E.append(atom.get_coord())

for res in chain_I:
    if res.has_id("CA"):
        for atom in res:
            atoms_I.append(atom.get_coord())

# Turn lists into numpy arrays with all the 3D atomic coordinates
coord_E = np.array(atoms_E)
coord_I = np.array(atoms_I)

## 2. Use these arrays to calculate the center of the chains and the direction from chain E to chain I. ##
# Use the .mean method along the COLUMNS (axis = 0) of the arrays
center_E = coord_E.mean(0)
center_I = coord_I.mean(0)

print "Coordinates of the chain E center: %s" % center_E
print "Coordinates of the chain I center: %s" % center_I

# Get the direction vector from enzyme to inhibitor
diff = center_I - center_E

# Normalize vector (make it unit length) Each coordinate = coord / sqrt(coord1^2 + coord2^2 + coord3^2))
direction = diff / np.linalg.norm(diff)

print "Direction of chain I from chain E: %s" % direction

## 3. Move the I chain 20 A in that direction. ##

for atom in chain_I.get_atoms():
    # Determine current position of its atoms
    position = atom.get_coord()
    # Set their new coordinates
    atom.set_coord(position + 20 * direction)

## 4. Write the structure in a PDB file, investigate the result in Pymol and save an image
# Create PDBIO object
io = PDB.PDBIO()

# Set the structure and save it
io.set_structure(s)
io.save('/home/juanma/Descargas/2PTC.split.out.pdb')

## ROTATE THE I CHAIN ##

## 5. Use the rotation matrix to rotate the inhibitor 90 degrees counterCW around its own center
##    (around the axis between enzyme and inhibitor)

def rotation_matrix(axis, theta):
    """ Return the rotation matrix associated with counterclockwise
    rotation about the given axis by theta radians. """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

# Direction of chain I from chain E
axis = [0.11879864,  0.90471053, -0.40912792]

# Rotation matrix object as a numpy array
matrix = rotation_matrix(axis, 0.5 * math.pi)
R = np.array(matrix)

# New coordinates of the inhibitor
atoms_I = []

for res in chain_I:
    if res.has_id("CA"):
        for atom in res:
            atoms_I.append(atom.get_coord())

coord_I = np.array(atoms_I)

# Find center of chain I again
center_I = coord_I.mean(0)

print "New coordinates of the chain I center: %s" % center_I

for atom in chain_I.get_atoms():
    # Determine current position of chain I atoms
    position = atom.get_coord()
    # Set their new coordinates: moving center to the origin, multiplying by R (rotating), and moving back
    rotate_position = R.dot(position - center_I) + center_I
    atom.set_coord(rotate_position)

## 6. Save the structure and check the rotation with Pymol

io.save('/home/juanma/Descargas/2PTC.split_rot.out.pdb')

## MODIFY THE B-FACTORS OF THE I CHAIN ##

## 7. Change the b-factor of the inhibitor and save the structure

for i, atom in enumerate(chain_I.get_atoms()):
    # Alteration of the b-factor following the formula, where i is each atom index of chain I
    b = np.cos(i / 100. * 2 * np.pi) * 20 + 25
    atom.set_bfactor(b)

io.save('/home/juanma/Descargas/2PTC_altered.out.pdb')