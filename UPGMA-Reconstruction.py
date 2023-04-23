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
file = open("distanceMatrix.txt", 'r')
readMatrix = []

for row in file:
    readMatrix.append([x for x in row.split()])

# Copying to DistanceMatrix object
distMatrix = calculator.get_distance(align)
for i in range(len(readMatrix)):
    for j in range(len(readMatrix)):
        distMatrix[i][j] = readMatrix[i][j]

# Creating a DistanceTreeConstructor object
constructor = DistanceTreeConstructor()

# Construct the phlyogenetic tree using UPGMA algorithm
UPGMATree = constructor.upgma(distMatrix)

# Draw the phlyogenetic tree
Phylo.draw(UPGMATree)

# Printing the phlyogenetic tree using terminal
# Phylo.draw_ascii(UPGMATree)