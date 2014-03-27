# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import scipy
import math
#import Pycluster as pc
import fnmatch
import os
import csv
import pickle
        
def load_mult_csv(files,delim):
    #import data from several different csv files. output to a list
    
    
    output=[]
    
    for i in range(len(files)):
       temp=[]
       with open(files[i], 'rb') as csvfile:
            read = csv.reader(csvfile, delimiter=delim, quotechar='|')
            first=1
            for row in read:
                
                if first==0:
                    temp.append(row)
                else:
                    first=0
                
                if not row:
                    print 'empty row'
       
       output.append(temp)
    
    return output

#def format_data():
output=[]
index=[16,17,18,19,20,21,22,23,24,97,98,99,100,101,102,103,104,105]
index=np.array(index)
#load the data
prefix='/media/HD1_/Documents/AG3C_local/experiments/glucose_time_course/sample'
files=[]
 
for i in index:
    
    for file in os.listdir(prefix+str(i)+'/RNA/'):
        if fnmatch.fnmatch(file, 'MURI_'+str(i)+'*.txt'):
            files.append(prefix+str(i)+'/RNA/'+file)
            
out=load_mult_csv(files,'\t')

ref=np.array(out[0]).T[1]
#average the responses

avindex=[]
for i in range(9):
    avindex.append([np.where(index==16+i)[0],np.where(index==97+i)[0]])


av=[]
std=[]
count=[]
not_observed=[]
file_len=[]
cfout=[]
not_ob_names=[]

for i in range(9):
    
    repeat=[]
    c_repeat=[]
    correction_factor_temp=[]
    not_ob_names_temp=[]
    for j in range(len(avindex[0])):
        
        
        
        
        temp=[float(k[0]) for k in np.array(out)[avindex[i][j]].T[3]]
        temp2=[float(k[0]) for k in np.array(out)[avindex[i][j]].T[2]]
        tempnames=[k[0] for k in np.array(out)[avindex[i][j]].T[1]]
        
        idx=np.where(np.array(temp)==0)[0]
        
        not_ob_names_temp.append(np.array(tempnames)[idx])
        
        
        not_observed=len(np.where(np.array(temp)==0)[0])
        file_len=len(temp2)
        #print not_observed
        
        correction_factor_temp.append(file_len/float(file_len-not_observed))
        
        repeat.append(temp)
        c_repeat.append(temp2)
    
    
    not_ob_names.append(not_ob_names_temp)
    print correction_factor_temp
    correction_factor=np.mean(correction_factor_temp)
    
    cfout.append(correction_factor)
    
    repeat=np.array(repeat).T
    c_repeat=np.array(c_repeat).T
    
    count.append([np.max(k) for k in c_repeat])
    temp_av=[]
    temp_std=[]

    for k in repeat:
        #k=np.array(k)
        if len(k[np.nonzero(k)])>0:
            temp_av.append(np.mean(k[np.nonzero(k)])*correction_factor)
            temp_std.append(np.std(k[np.nonzero(k)])*correction_factor)
            
        elif len(k[np.nonzero(k)])==1:
            temp_av.append(k[np.nonzero(k)]*correction_factor)
            temp_std.append(0.3*k[np.nonzero(k)]*correction_factor)
        else:
            temp_av.append(0.0)
            temp_std.append(0.0)
            
    av.append(temp_av)
    std.append(temp_std)

#print av[0]    
av=np.array(av).T
std=np.array(std).T
count=np.array(count).T

av=list(av)
std=list(std)
count=list(count)

#normalize results
for i in range(len(av)):
    std[i]=std[i]/(np.max(av[i]))
    av[i]=av[i]/(np.max(av[i]))

av_out_temp=[]
std_out_temp=[]
ref_out_temp=[]
#filter out proteins such that the minimum total count (over the full time course) is <7
for i in range(len(count)):
    
    if np.sum(count[i])>3*2*9.0:
        av_out_temp.append(list(av[i]))
        std_out_temp.append(list(std[i]))
        ref_out_temp.append(ref[i])

#filter out RNAs that don't change much (e.g. <1.5x)
av_out=[]
std_out=[]
ref_out=[]
for i in range(len(av_out_temp)):
    
    a=np.array(av_out_temp[i])
    if np.max(a)/np.min(a[np.nonzero(a)])>1.5 and np.max(a)>0.1:
        av_out.append(av_out_temp[i])
        std_out.append(std_out_temp[i])
        ref_out.append(ref_out_temp[i])
        
results=(av_out,std_out,ref_out)
pickle.dump(results,open('./RNAseq_data.p','wb'))

    #return av_out,std_out,ref_out

#av_out,std_out,ref_out=format_data()

