# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import csv
import pickle
from Bio import Entrez

def run(prefix,datatype,filename):
	
	mapping=pickle.load(open(prefix+datatype+'_name_map.p','rb'))

	item=[]
	not_mapped=[]

	with open(prefix+datatype+'_results/'+datatype+'_'+str(filename)+'_list.dat') as csvfile:
    		reader=csv.reader(csvfile,delimiter='\t')
    		for row in reader:
        		item.append(row)
        out_record=[]
	if len(item)>0:
		for i in item:
			if len(i)>0:

				#print mapping[i[0].lower()]
				if i[0].lower() in mapping:
					out_record.append(mapping[i[0].lower()])
				else:
					not_mapped.append(i[0].lower())

		with open(prefix+datatype+'_results/'+datatype+'_'+str(filename)+'_conv.txt','wb') as csvfile:
		    writer=csv.writer(csvfile,delimiter='\n')
		    
		    #for i in range(len(out_record)):
		    writer.writerow(out_record)
	
	pickle.dump(not_mapped,open(prefix+datatype+'_results/'+datatype+'_not_mapped_id_convert.p','wb'))

