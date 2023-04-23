# Importing necessary libraries
import sys
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

# Read the sequences and align
fileName = sys.argv[1]
align = AlignIO.read(fileName,'fasta')

# Calculate the distance matrix
calculator = DistanceCalculator('identity')

# Reading from input matrix
file = open("hardCodedDistance.txt", 'r')
readMatrix = []

for row in file:
    elements = []
    readMatrix.append([x for x in row.split()[1:]])

# Copying to DistanceMatrix object
distMatrix = calculator.get_distance(align)
for i in range(len(readMatrix)):
    for j in range(len(readMatrix)):
        distMatrix[i][j] = readMatrix[i][j]

# Creating a DistanceTreeConstructor object
constructor = DistanceTreeConstructor()

# Construct the phlyogenetic tree using NJ algorithm
NJTree = constructor.nj(distMatrix)

# Draw the phlyogenetic tree
Phylo.draw(NJTree)

# Printing the phlyogenetic tree using terminal
# Phylo.draw_ascii(NJTree)