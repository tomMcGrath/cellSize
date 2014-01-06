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
StringToUniprot = {} # StringToUniprot: StringID -> UniprotID (helps with abundance data, UniprotID is more widely used)

for line in UniprotFile:
    data = line.split('\t')
    StringToUniprot[data[0]] = data[1][:-1]
    UniprotToString[data[1][:-1]] = data[0]

# Start protein dictionary with abundances. If it doesn't exist on Uniprot, add its abundance to the abundanceError
for line in abundanceFile:
    if line[0] != '#':
        data = line.split('\t')
        try:
            proteinDict[StringToUniprot[data[1]]] = [float(data[2][:-1])]
        except KeyError:
            abundanceError += float(data[2])

# Process FASTA data (currently just gets AA length, degradation deprecated because M-excision & N-end mutually exclusive in prokaryotes)
FASTAFileContents = FASTAFile.read()
FASTAData = FASTAFileContents.split('>sp|')
for datum in FASTAData[1:]:
    header = datum.partition('\n')[0]
    sequence = datum.partition('\n')[2][:-1] # extract only the sequence data
    try:
        proteinDict[header[0:6]].append(header[7:])
        proteinDict[header[0:6]].append(len(sequence))
    except KeyError:
        print header[0:6], 'is not 1-1 with Uniprot ID' # may happen v.rarely
        continue
    
        

abundanceFile.close()
UniprotFile.close()
FASTAFile.close()
PSORTFile.close()
resultFile.close()
