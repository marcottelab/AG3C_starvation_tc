# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from readcsv import read_csv
import numpy as np
import pickle
import csvwriter
def run(prefix,datatype):

	filename=prefix+datatype+'_combined_conv.txt.chartReport'
	delim='\t'
	output=read_csv(filename,delim)

	#put genes in a given annotation cluster into a matrix (rows are different clusters, columns are different proteins)
	genes=[]
	annotations=[]
	fdr=[]

	for i in range(len(output)):
	    
	    if output[i][0]=='' or i==0: #we are in a new annotation cluster
		if i!=0:
		    genes.append(list(set(genes_temp))) #only keep unique set of genes in an annotation cluster
		    annotations.append(annotation_temp)
		    
		genes_temp=[]
		annotation_temp=[]
		
	    #check if we are on a row with annotation output
	    else:
		    
		if 'GO' in output[i][0]:
		    genes_temp.append(list(set(output[i][6].split(', '))))
		    annotation_temp.append(output[i][2])

	annotation=annotation_temp
	genes=genes_temp

	     
	#find all genes in a given annotation cluster that have a t1< than time t
	t=np.linspace(3,14*24,1000)
	#load the fit results
	(yout_av,yout_std,k_out_av,k_out_std,av_reference,a1)=pickle.load(open(prefix+datatype+'_fits.p','rb'))

	#get the t1 (initial inflection times)

	t1=np.array([0.]*len(k_out_av))
	for i in range(len(k_out_av)):
	    t1[i]=k_out_av[i]['tau1']

	#n is the count of genes (rows) that rise before a given time t (column)
	n=np.array([[0.]*len(t)]*len(genes))

	for i in range(len(t)):
	    #find the terms that have a t1<=t
	    idx=np.where(t1<=t[i])[0]
	    if len(idx)>0:
		ref_temp=[]
		for j in range(len(idx)):
		    ref_temp.append(av_reference[idx[j]])
		
		for j in range(len(genes)):
		    n[j][i]=len(list(set(map(str.lower,ref_temp)) & set(map(str.lower,genes[j]))))/float(len(genes[j]))


	#sort the go terms based upon when they hit 50% of max
	nclust_combine=[]
	annotation_combine=[]
	
	nclust_combine=n
	annotation_combine=annotation
	#for i in range(len(n_clust)):
	#    nclust_combine.append(n_clust[i])
	#    annotation_combine.append(annotation_clust[i])
	#for i in range(len(n_not_clust)):
	#    nclust_combine.append(n_not_clust[i])
	#    annotation_combine.append(annotation_not_clust[i])

	idx1=[0]*len(nclust_combine)
	idy1=[0]*len(nclust_combine)

	for i in range(len(nclust_combine)):
	    idx1[i]=np.where(np.min(np.abs(nclust_combine[i]-0.5))==np.abs(nclust_combine[i]-0.5))[0][0]
	    idy1[i]=i

	idx_sorted=np.argsort(idx1,axis=None)

	annotation_sorted=['']*len(idx1)
	time_rise=[]
	for i in range(len(idx1)):
	    annotation_sorted[i]=[annotation_combine[idx_sorted[i]]]
	    time_rise.append([annotation_sorted[i],t[idx1[idx_sorted[i]]]])

	
	filename=prefix+datatype+'_sorted_annotation_by_rise_time_rise.csv'
	csvwriter.write_csv(time_rise,filename,'\t')

