# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from readcsv import read_csv
import numpy as np
import pickle   
import csvwriter

filename='./RNA_all_combined_conv_list.txt.chartReport'
delim='\t'
output=read_csv(filename,delim)

output=np.array(output).T
#put genes in a given annotation cluster into a matrix (rows are different clusters, columns are different proteins)
genes=[]
annotations=[]
fdr=[]

for i in range(len(output)):
    genes_temp=[]
    annotation_temp=[]
    #if len(output[i])==0 or i==0: #we are in a new annotation cluster
    #    if i!=0:
    #        genes.append(list(set(genes_temp))) #only keep unique set of genes in an annotation cluster
    #        annotations.append(annotation_temp)
    #        
    #    genes_temp=[]
    #    annotation_temp=[]
        
    #check if we are on a row with annotation output
    #else:
            
    #if 'GO' in output[i][0]:
    if len(output[i])>0:
        genes_temp=genes_temp+output[i][5].split(', ')
        annotation_temp.append(output[i][1])

    genes.append(genes_temp)
    annotations.append(annotation_temp)
     
#find all genes in a given annotation cluster that have a t1< than time t
t=np.linspace(3,14*24,1000)
#load the fit results
(k_out_av,k_out_std,av_reference,err,a)=pickle.load(open('../fit_results_RNA_1.p','rb'))

#get the t1 (initial inflection times)

t1=np.array([0.]*len(k_out_av))
for i in range(len(k_out_av)):
    t1[i]=k_out_av[i]['tau1']
    
    
(mapped,out_record,out_common,not_mapped)=pickle.load(open('/media/HD1_/Documents/AG3C_local/PythonClient/RNA_analysis/RNA_all_name_mapping.p','rb'))

cluster_tau=[]
for i in genes:
    cluster_tau_temp=[]
    for j in i:
        idx=np.where(np.array(out_record)==j)[0]
        if len(idx)>0:
            idy=np.where(np.array(av_reference)==mapped[idx[0]])[0]
            cluster_tau_temp.append(t1[idy[0]])
        
    cluster_tau.append(cluster_tau_temp)
    
T=np.linspace(3,14*24,1000)
n_clust=[]
#find the number of genes in a cluster with a tau1 <t
for i in cluster_tau:
    n_clust_temp=[]
    for t in T:
        n_clust_temp.append(len(np.where(np.array(i)<t)[0]))
    n_clust.append(np.array(n_clust_temp)/float(np.max(n_clust_temp)))
    
#sort the go terms based upon when they hit 50% of max
nclust_combine=[]
annotation_combine=[]

idx1=[0]*len(n_clust)
idy1=[0]*len(n_clust)

for i in range(len(n_clust)):
    idxx=np.where(np.min(np.abs(n_clust[i]-0.5))==np.abs(n_clust[i]-0.5))[0]
    if len(idxx)>0:
        idx1[i]=idxx[0]
        idy1[i]=i

idx_sorted=np.argsort(idx1,axis=None)

annotation_sorted=['']*len(idx1)
time_rise=[]
for i in range(len(idx1)):
    annotation_sorted[i]=[annotations[idx_sorted[i]]]
    time_rise.append([T[idx1[idx_sorted[i]]]])

writeoutput=[]
for i in range(len(time_rise)):
    if len(annotation_sorted[i][0])>0:
        writeoutput.append((annotation_sorted[i][0][0],time_rise[i][0]))

filename='./RNA_sorted_annotation_by_rise_time_rise_chart_all.csv'
csvwriter.write_csv(writeoutput,filename,'\t')

