import os
import sys

def getlevels(dataset,column):
    #levels = {}
    options = []
    levels = []
    index=int(column)-1
    
    file = open(dataset.get_file_name(), 'r')
    line = file.readline()
    
    cpt=0
    for line in file:
        levels.append(line.split('\t')[index].strip('"'))
        cpt=cpt+1        
        if cpt > 100:
            break

    levels = list(set(levels))
    levels.sort();

    for level in levels:
        #option = {'name':level,'value':level,'options':[],'selected':0}
	options.append((level,level, False))
    file.close()
    
    return options

