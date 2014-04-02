# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import csv
import pickle


from Bio import Entrez

##Map the RNA list appropriately (to K-12 strain) to use in david. 

# *Always* tell NCBI who you are
Entrez.email = "johnrhouser@utexas.edu"

out_record=[]
not_mapped=[]
mapped=[]
out_common=[]

exit=False

item=[]
with open('../analysis/fit_stepwise/RNA_classify/RNA_all_combined_list.dat') as csvfile:
    reader=csv.reader(csvfile,delimiter='\t')
    for row in reader:
        item.append(row)
        
#item = ['cfa','pntB','gcvR','pflA','fepC','yadR']
animal = 'Escherichia coli' 
for i in item:
    if len(i)>0:
        if 'RF' not in i[0]:
                
            search_string = i[0]+'[Gene] Escherichia coli[Organism]'
            
            handle = Entrez.esearch(db="Gene", term=search_string)
            record = Entrez.read(handle)
            ids = record['IdList']
            
            if len(ids)>0:
                handle = Entrez.efetch(db="Gene", id=ids[0],retmode="xml") 
                records=Entrez.parse(handle)
            
                
                for record in records:
                #print record['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
                    rtemp=record['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
                #out_record.append(rtemp)
                    
            
                handle = Entrez.efetch(db="Gene", id=rtemp,retmode="xml") 
                records=Entrez.parse(handle)
                record in records
                try:
                    name=record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
                    
                except:
                    not_mapped.append(i)
                    
                search_string = str(name)+'[Gene] Escherichia coli[Organism]'
                handle = Entrez.esearch(db="Gene", term=search_string) 
                record = Entrez.read(handle)
                ids = record['IdList']
                
                for k in range(len(ids)):
                    if ids[k][0]=='9' and ids[k][1]=='4':
                        out_record.append(ids[k])
                        exit=True
                        break
                
                
                if exit==False:
                    not_mapped.append(i)
                else:
                    out_common.append(name)
                    mapped.append(i)
                    
                exit=False
                
            else:
                not_mapped.append(i)
            #break
        else:
            not_mapped.append(i)
            

with open('./RNA_analysis/RNA_all_combined_conv_list.txt','wb') as csvfile:
    writer=csv.writer(csvfile,delimiter='\n')
    
    #for i in out_record:
    writer.writerow(out_record)
    
pickle.dump((mapped,out_record,out_common,not_mapped),open('./RNA_analysis/RNA_all_name_mapping.p','wb'))

