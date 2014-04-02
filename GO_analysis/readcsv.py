# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
import csv

def read_csv(filename,delim):
   
    output=[]
    with open(filename, 'rb') as csvfile:
    
        read = csv.reader(csvfile, delimiter=delim)
        for row in read:
            output_temp=[]
            for j in range(len(row)):
                output_temp.append(row[j])
            output.append(output_temp)
    
    return output
# <codecell>


