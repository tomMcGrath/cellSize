# Data files
abundanceFile = open('511145-E.coli_whole_organism-integrated_dataset.txt', 'r')
UniprotFile = open('ecoli_conversion.tab', 'r')
FASTAFile = open('ecoli_sequences.fasta', 'r')
PSORTFile = open('psort_ecoli.txt', 'r')

# Output files
resultFile = open('results.txt', 'w')

#Global variables
cellLife = 3000.0 # seconds (3000 for wt e.coli)
ribosomeNumber = 20000 # dimensionless; ribosomes per cell (20000 for wt e. coli)
ribosomeSynthesisRate = 25 # amino acids per ribosome per second (25 for wt e. coli)
abundanceError = 0.0 # parts per million
proteinsPerCell = 5000000 # proteins per cell (3-10e6 for wt e. coli ref BioNumbers)

#Data dictionary
proteinDict = {} # proteinDict: UniprotID -> [abundance, details, length, localisation, localisation certainty, multiple locations?, absolute abundance, energy cost]
StringToUniprot = {} # StringToUniprot: StringID -> UniprotID (helps with abundance data, UniprotID is more widely used)
UniprotToString = {}

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
FASTAData = FASTAFileContents.split('>')
for datum in FASTAData[1:]:
    header = datum.partition('\n')[0]
    sequence = datum.partition('\n')[2][:-1] # extract only the sequence data
    try:
        proteinDict[header[3:9]].append(header[7:])
        proteinDict[header[3:9]].append(len(sequence))
    except KeyError:
        print header[3:9], 'is not 1-1 with Uniprot ID (seq length)' # may happen v.rarely
        continue

# Process localisation data from PSORT-B (needs updating for eukaryotes)
PSORTFileContents = PSORTFile.read()
PSORTData = PSORTFileContents.split('SeqID: ')

for datum in PSORTData[1:]:
    header = datum.partition('\n')[0][3:9]
    headerData = datum.partition('\n')[0][9:]
    analysis = datum.partition('\n')[2]

    result = analysis.partition('Final Prediction:')
    locationData = result[2].split('\n')[1]

    if locationData.find('Unknown') != -1:
        location = 'Unknown'
        probability = 10.0
        possMany = False

    else:
        location = locationData.split(' ')[4]
        probability = locationData.split(' ')[-1][:-1]
        if locationData.split(' ')[5] != '':
            possMany = True
        else:
            possMany = False
    try:
        proteinDict[header].append(location)
        proteinDict[header].append(float(probability))
        proteinDict[header].append(possMany)

    except KeyError:
        print header[0:6], 'is not 1-1 with Uniprot ID (localisation)' # may happen v.rarely
        continue

# Calculate absolute protein abundance
totalAbundance = 0.0
totalCost = 0.0

for key in proteinDict.keys():
    try:
        totalAbundance += proteinDict[key][0]*proteinsPerCell/1000000 # converts from PPM to absolute abundance
        proteinCost = 4*(proteinDict[key][0]*proteinsPerCell/1000000)*proteinDict[key][2]/cellLife # calculates cost of this protein in ATP/s
        totalCost += proteinCost
        
        proteinDict[key].append(proteinDict[key][0]*proteinsPerCell/1000000)
        proteinDict[key].append(proteinCost)

    except IndexError:
        print key

# Calculate % cost
for key in proteinDict.keys():
    proteinDict[key].append(proteinDict[key][7]*100/totalCost)

# Calculate distribution %
cytoSum = 0.0
membraneSum = 0.0
peripSum = 0.0
cytoMemSum = 0.0
unknownSum = 0.0

for key in proteinDict.keys():
    if proteinDict[key][3] == 'Cytoplasmic':
        cytoSum += proteinDict[key][8]
    elif proteinDict[key][3] == 'OuterMembrane':
        membraneSum += proteinDict[key][8]
    elif proteinDict[key][3] == 'Periplasmic':
        peripSum += proteinDict[key][8]
    elif proteinDict[key][3] == 'CytoplasmicMembrane':
        cytoMemSum += proteinDict[key][8]
    elif proteinDict[key][3] == 'Unknown':
        unknownSum += proteinDict[key][8]

print 'Cytoplasmic percentage: ', cytoSum + unknownSum/2
print 'Membrane & periplasm percentage: ', membraneSum + peripSum + cytoMemSum + unknownSum/2
        
abundanceFile.close()
UniprotFile.close()
FASTAFile.close()
PSORTFile.close()
resultFile.close()
