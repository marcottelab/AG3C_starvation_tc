# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

##author: John Houser
##format, average, and normalize AG3C data

import csv
import numpy as np
import math
import scipy
import math
import pickle
import numpy
#from collections import defaultdict
#import dictlist

class Dictlist(dict):
    def __setitem__(self, key, value):
        try:
            self[key]
        except KeyError:
            super(Dictlist, self).__setitem__(key, [])
        self[key].append(value)

class ngram(dict):
    """Based on perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return super(ngram, self).__getitem__(item)
        except KeyError:
            value = self[item] = type(self)()
            return value

class format_data(object):
    def __init__(self,delim='\t'):
        #self.files=files
        self.delim=delim
        self.idata=[]
        self.out=[]
        self.minchange=1.5
        self.cutoff=10
   
    def main(self):
        
        
        #get rid of peptides that don't have a minimum number of observed peptides for the whole time course
        self.cleandata()
	
	#convert raw counts to fractional counts
	self.raw2fraction()
	
        #average the output data
        self.average()
	#print self.avout['YP_003046375.1']
	#keep only those changing by >1.5 through the time cours
	self.filter()
        #normalize to initial value
  	self.normalize()
      
            
    def load_mult_csv(self,files,delim):
        #import data from several different csv files. output to a list
        
        
        output=dict()
	reference_cn=dict()
        with open(files, 'rb') as csvfile:
            read = csv.reader(csvfile, delimiter=delim, quotechar='|')
            first=1
            for row in read:
                    
                if first==0:
			if len(row)==4:
                		output[row[0]]=row[3]
                	else:
				output[row[0]]='0'
			
			reference_cn[row[0]]=row[1]
		else:
                        first=0
                    	#print row

                if not row:
                        print 'empty row'
           
        
	self.common_name=reference_cn

        return output

    def raw2fraction(self):
        
	T=self.time
        #print self.idata
        #outdata=numpy.array(self.idata)
        outdata=self.idata
	'''
	for j in range(3):
		if j==0:
		    ref=self.sample_refa_r
		elif j==1:
		    ref=self.sample_refb_r
		elif j==2:
		    ref=self.sample_refc_r
		#print ref
                sit=[]
		for i in outdata:
		    temps=1.0
		    count=0
		    for t in T:
			#print outdata[i][t],ref
		        k=ref[str(t)]
			#print str(k),outdata[i][t][str(k)]
			if type(outdata[i][t][str(k)])!=type(0.0):
				outdata[i][t][str(k)]=0.0
			if outdata[i][t][str(k)]!=0:
			        temps*=float(outdata[i][t][str(k)])
				count+=1
		    s_itj=[]
		    if count==0:
		        count=1

		    for t in T:
			k=ref[str(t)]	    
		     	s_itj.append(outdata[i][t][str(k)]/(temps**(1/count)))
		    
		    sit.append(s_itj)
		#print sit
		sit=np.array(sit).T
		sit=[s[np.ndarray.nonzero(s)] for s in sit]
		#print sit
	        norm_s=[np.median(s) for s in sit]

		norm_sd=dict()
		for t in range(len(T)):
		 	norm_sd[str(T[t])]=norm_s[t]
			

		for i in outdata:
		    for t in T:
			k=ref[str(t)]
			
			outdata[i][t][str(k)]=float(outdata[i][t][str(k)])/norm_sd[str(t)]
	
	'''
	for t in T:
	    for k in self.sample_ref_r[str(t)]:
		norm=0
        	for i in outdata:
			#print outdata[i][t][k]
			if type(outdata[i][t][k])!=type(0.0):
				outdata[i][t][k]=0.0
			else:
				norm+=outdata[i][t][k]             
			        	
		for i in outdata:
			outdata[i][t][k]=outdata[i][t][k]/norm	
		
        self.output=outdata
    
    def average(self):
	
	data=self.output
	T=self.time
	
	outav=ngram()
	outst=ngram()

	for i in data:
		for j in data[i]:
			dv=np.array(data[i][j].values())
			
			outav[i][j]=np.mean(dv[np.nonzero(dv)])			
			outst[i][j]=np.std(dv[np.nonzero(dv)])		
	
	self.avout=outav
	self.stdout=outst
    
    def filter(self):
	
	data=self.avout
	stdata=self.stdout
	
	nonresponders=[]
	for i in data:
		dv=np.array(data[i].values())
		if np.max(dv[np.nonzero(dv)])/np.min(dv[np.nonzero(dv)])<1.5:
			nonresponders.append(i)
	
	for i in nonresponders:
		data.pop(i,0)
		stdata.pop(i,0)
	
	self.avout=data
	self.stdout=stdata
	
    def normalize(self):
	
	data=self.avout
	stdata=self.stdout
	
	#dataout=ngram()
	for i in data:
		dv=np.array(data[i].values())
		norm=np.max(np.nan_to_num(dv[np.nonzero(dv)]))
		
		for j in data[i]:
			if norm!=0:			
				data[i][j]=np.nan_to_num(data[i][j])/norm
				stdata[i][j]=np.nan_to_num(stdata[i][j])/norm
			else:
				data[i][j]=0.5
				stdata[i][j]=0.5
	self.avout=data
	self.stdout=stdata
	
    def cleandata(self):
    	removelist=[]
        cutoff=self.cutoff
        data=self.idata
	for i in data:
		totals=0
		for j in data[i]:
			for k in data[i][j]:
				totals+=data[i][j][k]
	
		if totals<cutoff:
			#data.pop(i,0)
			removelist.append(i)
	for i in removelist:
		data.pop(i,0)

        self.output=data
            

def run():
	#load the data
	t=[3,4,5,6,8,24,48,168,336]
	#create a dictionary to map the sample number back to the time point.
	sample_ref=dict()
	sample_refa_r=dict()
	sample_refb_r=dict()
	sample_refc_r=dict()
	for i in range(9):
		sample_ref[str(16+i)]=t[i]
		sample_refa_r[str(t[i])]=str(16+i)
		sample_ref[str(25+i)]=t[i]
		sample_refb_r[str(t[i])]=str(25+i)
		sample_ref[str(97+i)]=t[i]
		sample_refc_r[str(t[i])]=str(97+i)

	#reverse the index
	sample_ref_r=Dictlist()
	for i in sample_ref:
		#if str(sample_ref[i]) in sample_ref_r:
		#	sample_ref_r[str(sample_ref[i])].append(i)
		#else:
		sample_ref_r[str(sample_ref[i])]=i
		#sample_ref_r.setdefault(str(sample_ref[i]),[])

	#print sample_ref_r
	D=format_data()
	D.time=t
	D.sample_ref_r=sample_ref_r
	D.sample_refa_r=sample_refa_r
	D.sample_refb_r=sample_refb_r
	D.sample_refc_r=sample_refc_r

	out=ngram()
	#load the raw spectral counts
	#for j in t:
	    #for i in sample_ref_r[str(j)]:
	prefix='/media/HD1_/Documents/AG3C_data/data/experiments/time_course_P0/'
        files=prefix+'/glu_all_reps_raw_count_nor.csv'
	out_temp=D.load_mult_csv(files,',')
	    	
	for k in range(len(out_temp)):
		for j in range(len(out_temp[k])))
			if k!=0: 
				a=''
				end=False
				for text in out_temp[0][j]:				
					if end=True:
						b=text
					if text!='_':
						a=a+text
						end=True
				out[k][a][b]=float(out_temp[k])
	
	#print out
	D.idata=out
	D.main()

	avout=[]
	stdout=[]
	refout=[]
	for i in D.avout:
		
		if i in D.common_name:
			if D.common_name[i]=='':
				refout.append(i)
			else:	
				refout.append(D.common_name[i])
			avout.append(D.avout[i].values())
			stdout.append(D.stdout[i].values())
		
	

	pickle.dump((avout,stdout,refout),open('/media/HD1_/Documents/AG3C/Results/formated_data.p','wb'))



if __name__ == "__main__":

	run()