## EXERCISE 1: Write a script to sort the base pairs and their probabilities from a .eps (parsed) file ##

import numpy as np

# Parsed file with the following structured lines: 1st base pair element, 2nd base pair element, probability(base pair)
file = open("/home/juanma/PycharmProjects/Bioinformatics/sequence1_parsed.txt").readlines()

# List of lists with every bp and its probability
bp = []

for line in file:
    bp.append(line.split())

# Get the numpy array
bp_np = np.array(bp)

# Sort the bp by their probabilities (reversed)
inverse_bp = bp_np[bp_np[:, 2].argsort()]

# Base-pairs sorted from highest to lowest probability
ordered_bp = inverse_bp[::-1]

# Printed array (view and debugging)
# print ordered_bp

# List to store probabilities
probs = []

# Get only probabilities
for i in bp:
    probs.append(float(i[2]))

counter = 0

for i in probs:
    if i > 0.8:
        counter += 1

print "BP with p > 80%: " + str(counter)


## EXERCISE 3: Write a script to compute the Hamming and base-pair distances between two RNA seqs. of equal length ##

# sequence_A_full = ('..(((((((....)))))))  .....(((((((((((.((..(((...)))..)).)))))))))))....')
# sequence_B_full = ('..(((((((....)))))))  .(((((((..(((((((.....)))))))......)))))))........')
# sequence_consen = ('..(((((((....)))))))  .((..((((((((((((((((((...)))))))..)))))))))))..))')

# As can be seen, the first 20 bases of the two sequences A and B form the same pairs.
# So this part of both sequences can be removed for easier scripting terms.

# Paste here the two sequences to be compared
seq_A = ('.....(((((((((((.((..(((...)))..)).)))))))))))....')
seq_B = ('.(((((((..(((((((.....)))))))......)))))))........')

# Counter for differing letters
ham_dist = 0

# Element-wise comparison: If a letter is different between the two seqs. in a certain position, add 1 to the counter
for i in range(len(seq_A)):
    if seq_A[i] != seq_B[i]:
        ham_dist += 1

print "Hamming distance for the given sequences: " + str(ham_dist)

# Lists for the indexes of the first '(' and second ')' elements of each base-pair in sequence A
seq_A_open = []
seq_A_closed = []

for i in range(len(seq_A)):
    if seq_A[i] == ('('):
        seq_A_open.append(i)
    elif seq_A[i] == (')'):
        seq_A_closed.append(i)

# Reverse the order of sequence A closed list to form the adequate base-pairs
seq_A_closed_rev = seq_A_closed[::-1]

# Form base pairs and present them as list of tuples
seq_A_bp = ([(seq_A_open[i], seq_A_closed_rev[i]) for i in range(len(seq_A_open))])

# Printed list of bp in seq. A (view and debugging)
# print "List of base pairs in sequence A: " + str(seq_A_bp)

# Lists for the indexes of the first '(' and second ')' elements of each base-pair in sequence B
seq_B_open = []
seq_B_closed = []

for i in range(len(seq_B)):
    if seq_B[i] == ('('):
        seq_B_open.append(i)
    elif seq_B[i] == (')'):
        seq_B_closed.append(i)

# Reverse the order of the sequence B closed list to form the adequate base-pairs
seq_B_closed_rev = seq_B_closed[::-1]

# Form base-pairs and present them as list of tuples
seq_B_bp = ([(seq_B_open[i], seq_B_closed_rev[i]) for i in range(len(seq_B_open))])

# Printed list of bp in seq. B (view and debugging)
# print "List of base pairs in the seq_B sequence: " + str(seq_B_bp)

# Differing base pairs IN BOTH sequences counter
basepair_distance = 0

# For each base-pair present in sequence A and not in sequence B, add 1 to the counter
for i in seq_A_bp:
    if i not in seq_B_bp:
        basepair_distance += 1

# For each base-pair present in sequence B and not in sequence A, add 1 to the counter
for i in seq_B_bp:
    if i not in seq_A_bp:
        basepair_distance += 1

print "Base pair distance for the given sequences: " + str(basepair_distance)