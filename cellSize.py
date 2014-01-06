# Data files
abundanceFile = open('511145-E.coli_whole_organism-integrated_dataset.txt', 'r')
UniprotFile = open('2013121951YKSHD5B0.tab', 'r')
FASTAFile = open('ecoli_sequences.fasta', 'r')
PSORTFile = open('psort.txt', 'r')

# Output files
resultFile = open('results.txt', 'w')

#Global variables
cellLife = 3000 # seconds
abundanceError = 0.0 # parts per million

#Data dictionary
proteinDict = {} # proteinDict: UniprotID -> [abundance, localisation, length, cost]
StringToUniprot = {}

for line in UniprotFile:
    data = line.split('\t')
    StringToUniprot[data[0]] = data[1][:-1]

for line in abundanceFile:
    if line[0] != '#':
        data = line.split('\t')
        try:
            proteinDict[StringToUniprot[data[1]]] = float(data[2][:-1])
        except KeyError:
            abundanceError += float(data[2])
        

abundanceFile.close()
UniprotFile.close()
FASTAFile.close()
PSORTFile.close()
resultFile.close()
