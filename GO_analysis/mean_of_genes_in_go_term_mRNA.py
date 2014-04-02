# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from readcsv import read_csv
import numpy as np
import pickle
import pylab as plt

prefix='/media/HD1_/Documents/AG3C_local/PythonClient/RNA_analysis/'

(mapped,out_record,out_common,not_mapped)=pickle.load(open(prefix+'RNA_turn_on_name_mapping.p','rb'))

filename=prefix+'RNA_turn_on_combined_conv_list.txt.chartReport'
delim='\t'
output=read_csv(filename,delim)

#put genes in a given annotation cluster into a matrix (rows are different clusters, columns are different proteins)
genes=[]
annotations=[]
fdr=[]

for i in range(len(output)):
    
    if i!=0:
        if float(output[i][12])<=0.05:
            genes.append(output[i][5].split(', '))
            annotations.append(output[i][1])
        
#get the preped data and average the protein time courses for that go-term
(avout,stdout,av_reference,err,accept)=pickle.load(open('../RNAseq_av_part.p','rb'))
av_ofav_ingo=[]
for i in range(len(genes)):
    a=list(set(map(str.lower,genes[i])))
    temp_avs=[]
    for j in range(len(a)):
        idx=np.where(np.array(out_record)==a[j])[0][0]
        idy=np.where(np.array(av_reference)==mapped[idx][0])[0][0]
        temp_avs.append(avout[idy])
    
    av_ofav_ingo.append(np.average(temp_avs,axis=0))
    
    
#plot results

t=[3,4,5,6,8,12,24,24*7,14*24]

av2=np.ndarray.transpose(np.array(av_ofav_ingo))
plt.figure()
plt.semilogx(t,av2)

#plt.legend(annotations,bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
plt.xlabel('time(hr)')
plt.ylabel('average relative expression')
plt.show()

