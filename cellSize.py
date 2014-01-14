import os

# Data files
abundanceFile = open('4932-S.cerevisiae_whole_organism-integrated_dataset.txt', 'r')
UniprotFile = open('scerevisiae_conversion.tab', 'r')
FASTAFile = open('scerevisiae_sequences.fasta', 'r')
localisationFile = open('LocTreeData/559292_Saccharomyces_cerevisiae.euka.lc3', 'r')
localisationType = 'LocTree3'

# Output files
resultFile = open('scerevisiae_results.txt', 'w')

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

# Process localisation data from localisation file according to formatting from source tool
if localisationType == 'PSORTb':
    localisationFileContents = localisationFile.read()
    localisationData = localisationFileContents.split('SeqID: ')

    for datum in localisationData[1:]:
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

elif localisationType == 'LocTree3':
    locErrorCount = 0
    for line in localisationFile:
        if line[0] == '#':
            continue # ignore header files
        data = line.split('\t')
        proteinID = data[0][3:9]
        localisation = data[2]
        probability = data[1]
        possMany = False # LocTree3 doesn't seem to account for this

        try:
            proteinDict[proteinID].append(localisation)
            proteinDict[proteinID].append(float(probability)/10)
            proteinDict[proteinID].append(possMany)

        except KeyError:
            locErrorCount += 1
            #print proteinID, 'is not 1-1 with Uniprot ID (localisation)' # may happen v.rarely
            continue

# Calculate absolute protein abundance
totalAbundance = 0.0
totalCost = 0.0
abundErrorCount = 0

for key in proteinDict.keys():
    try:
        totalAbundance += proteinDict[key][0]*proteinsPerCell/1000000 # converts from PPM to absolute abundance
        proteinCost = 4*(proteinDict[key][0]*proteinsPerCell/1000000)*proteinDict[key][2]/cellLife # calculates cost of this protein in ATP/s
        totalCost += proteinCost
        
        proteinDict[key].append(proteinDict[key][0]*proteinsPerCell/1000000)
        proteinDict[key].append(proteinCost)

    except IndexError:
        abundErrorCount += 1
        #print key

# Calculate % cost
costErrorCount = 0
for key in proteinDict.keys():
    try:
        proteinDict[key].append(proteinDict[key][7]*100/totalCost)

    except IndexError:
        costErrorCount += 1
        #print key

# Calculate distribution %
cytoSum = 0.0
membraneSum = 0.0
peripSum = 0.0
cytoMemSum = 0.0
extraCellularSum = 0.0
cellWallSum = 0.0
unknownSum = 0.0

localisations = {}
for key in proteinDict.keys():

    try:
        localisations[proteinDict[key][3]] += proteinDict[key][8]

    except KeyError:
        try:
            localisations[proteinDict[key][3]] = proteinDict[key][8]
        except IndexError:
            print key

    except IndexError:
        print key

for key in localisations.keys():
    print key, localisations[key]


##    if localisationType == 'PSORTb':
##        if proteinDict[key][3] == 'Cytoplasmic':
##            cytoSum += proteinDict[key][8]
##        elif proteinDict[key][3] == 'OuterMembrane':
##            membraneSum += proteinDict[key][8]
##        elif proteinDict[key][3] == 'Periplasmic':
##            peripSum += proteinDict[key][8]
##        elif proteinDict[key][3] == 'CytoplasmicMembrane':
##            cytoMemSum += proteinDict[key][8]
##        elif proteinDict[key][3] == 'Unknown':
##            unknownSum += proteinDict[key][8]
##        elif proteinDict[key][3] == 'Extracellular':
##            extraCellularSum += proteinDict[key][8]
##        elif proteinDict[key][3] == 'Cellwall':
##            cellWallSum += proteinDict[key][8]
##        else:
##            print proteinDict[key][3]
##
##    if localisationType == 'LocTree3':
##        if proteinDict[key][3] == 'nucleus':
##            nucleusSum += proteinDict[key][8]
##        elif proteinDict[key][3] == 'cytoplasm':
##            cytoSum
##
##print 'Cytoplasmic percentage: ', cytoSum
##print 'Outer membrane percentage: ', membraneSum
##print 'Periplasmic percentage: ', peripSum
##print 'Cytoplasmic membrane percentage: ', cytoMemSum
##print 'Cell wall percentage: ', cellWallSum
##print 'Extracellular percentage: ', extraCellularSum
##print 'Unknown percentage: ', unknownSum
##print 'Bulk k: ', cytoSum + unknownSum/2 + extraCellularSum
##print 'Radial k: ', membraneSum + peripSum + cytoMemSum + unknownSum/2 + cellWallSum

# Output data to file
resultFile.write('BEGIN RESULTS HEADER:\n')
for key in localisations.keys():
    resultFile.write(key + str(localisations[key]))

##resultFile.write('Cytoplasmic percentage: ' + str(cytoSum) + '\n')
##resultFile.write('Outer membrane percentage: ' + str(membraneSum))
##resultFile.write('Periplasmic percentage: ' +str(peripSum))
##resultFile.write('Cytoplasmic membrane percentage: ' + str(cytoMemSum))
##resultFile.write('Cell wall percentage: ' + str(cellWallSum) + '\n')
##resultFile.write('Extracellular percentage: ' + str(extraCellularSum) + '\n')
##resultFile.write('Unknown percentage: ' + str(unknownSum) + '\n')
##resultFile.write('Bulk k: ' + str(cytoSum + unknownSum/2 + extraCellularSum) + '\n')
##resultFile.write('Radial k: ' + str(membraneSum + peripSum + cytoMemSum + unknownSum/2 + cellWallSum) + '\n')
##resultFile.write('END OF RESULTS HEADER\n')

for key in proteinDict.keys():
    resultText = key
    for data in proteinDict[key]:
        resultText += (','+ str(data))
    resultText += '\n'

    resultFile.write(resultText)

        
abundanceFile.close()
UniprotFile.close()
FASTAFile.close()
localisationFile.close()
resultFile.close()
