# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pylab as plt
import numpy as np
from scipy import stats
import pickle
import csv

def classify(k,kstd):
    
    err=0.1
    
    #case1 turns on stays on
    cond1=(k['A3']-kstd['A3']-err)>(k['A2']+kstd['A2']) and (k['A2']+kstd['A2']+err)>=(k['A1']-kstd['A1'])
    cond2=(k['A2']-kstd['A2']-err)>(k['A1']+kstd['A1']) and (k['A3']+kstd['A3']+err)>=(k['A2']-kstd['A2'])
    cond2B=(k['A3']-kstd['A3']-err)>(k['A1']+kstd['A1']) and (k['A2']+kstd['A2']+err)>=(k['A1']-kstd['A1']) and (k['A2']-kstd['A2']-err)<=(k['A1']+kstd['A1'])

    #case2 turns off and stays off
    cond3=k['A3']+kstd['A3']+err<k['A2']-kstd['A2'] and k['A2']-kstd['A2']<=k['A1']+kstd['A1']+err
    cond4=k['A2']+kstd['A2']+err<k['A1']-kstd['A1'] and k['A3']-kstd['A3']<=k['A2']+kstd['A2']+err
    cond4B=k['A3']+kstd['A3']+err<k['A1']-kstd['A1'] and k['A2']-kstd['A2']<=k['A1']+kstd['A1']+err and k['A2']+kstd['A2']>=k['A1']-kstd['A1']-err

    #case3 pulsed on
    cond5=k['A2']-kstd['A2']>k['A1']+kstd['A1']+err and k['A3']+kstd['A3']+err<k['A2']-kstd['A2']
    
    
    #case4 pulsed off
    cond6=k['A2']+kstd['A2']+err<k['A1']-kstd['A1'] and k['A3']-kstd['A3']>k['A2']+kstd['A2']+err
    
    classify='other'
    if cond1 or cond2 or cond2B:
        classify='on'
        
        
    elif cond3 or cond4 or cond4B:
        classify='off'
        
    elif cond5:
        
        classify='pulse_up'    
       
    elif cond6:
        
        classify='pulse_down'
         
    else:
        classify='other'
        
    return classify

def run(prefix,datatype):

	(yav,ystd,k_out_av,k_out_std,av_reference,a1)=pickle.load(open(prefix+datatype+'_fits.p','rb'))
	

	I=len(k_out_av)
	classification=['']*I
	for i in range(I):
	    	classification[i]=classify(k_out_av[i],k_out_std[i])
	
	pickle.dump(classification,open(prefix+datatype+'_classify.p','wb'))
	
	classes=['on','off','pulse_up','pulse_down','other']
	
	for i in range(5):
		groupn=[]
		for j in range(I):
			if classification[j]==classes[i]:
				groupn.append(av_reference[j])
	#list_type=[]
		#print groupn
		with open(prefix+datatype+'_'+classes[i]+'_list.dat','wb') as csvfile:
        		writer =csv.writer(csvfile,delimiter='\n')
       			writer.writerow(groupn)

	
	with open(prefix+datatype+'_combined_list.dat','wb') as csvfile:
        	writer =csv.writer(csvfile,delimiter='\n')
       		writer.writerow(av_reference)

if __name__=="__main__":
	run()

