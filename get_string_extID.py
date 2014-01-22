infile = open('10090-integrData_mouse_liver.txt', 'r')
outfile = open('string_IDs.txt', 'w')

for line in infile:
    if line[0] == '#':
        continue
    
    else:
        data = line.split('\t')
        outfile.write(data[1] + '\n')

infile.close()
outfile.close()
        
