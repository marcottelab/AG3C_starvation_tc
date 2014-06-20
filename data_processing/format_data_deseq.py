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
        #self.cleandata()
	
	#convert raw counts to fractional counts
	self.raw2fraction()	
        
	#get rid of genes that don't ever change signifcantely throughout the time course
	self.cleandata()
	
	#average the output data
        self.average()
        #normalize to initial value
  	self.normalize()
  

    def cleandata(self):
	
	pval=self.pval
	output=self.output
	#print len(output)
	#print output
	not_significant=[]
	out_clean=ngram()
	for i in pval:
		#print output[i]
		if str(i) in output:
			out_clean[str(i)]=output[str(i)]
	#for i in pval[0:3:1]:
	#	print i
	#for i in output[0:3:1]:
	#	print i

	#print out_clean
	self.output=out_clean

    def load_mult_csv(self,files,delim):
        #import data from several different csv files. output to a list
        
        
        output=ngram()
	reference_cn=dict()
        with open(files, 'rb') as csvfile:
            read = csv.reader(csvfile, delimiter=delim, quotechar='"')
            first=1
	    for row in read:
                if first==1:
			row0=row
	
		#print row
		for j in range(len(row)-1):
			end=False
			a=''
			for k in row0[j+1]:
				if end==True:
					b=k
				if k!='_' and end==False:
					a=a+str(k)
				else:
					end=True
		        if first!=1:
			
				output[row[0]][a][b]=float(row[j+1])
		
		reference_cn[row[0]]=row[0]
		if first==1:
			first=0
           
        
	self.common_name=reference_cn

        return output

    def load_size_factors(self,files,delim):
        #import data from several different csv files. output to a list
        
        T=self.time
        output=ngram()
	reference_cn=dict()
	J=0
	norm=1.0
	for f in files:
		
		with open(f, 'rb') as csvfile:
		    read = csv.reader(csvfile, delimiter=delim, quotechar='"')
		    first=1
		    J+=1
		    for row in read:
			if first==1:
				row0=row
		
			#print row
			for j in range(len(row)-1):
				if first!=1:
					if 't3_' in row[0]:
						if row[0][3:]=='1':
							norm=float(row[1])
						output[row[0]][str(T[J])]=float(row[1])/norm
						
					else:
						output[row[0]]=float(row[1])/norm
			if first==1:
				first=0
           
        

        return output
    
    def load_p_values(self,files,delim):
        #import data from several different csv files. output to a list
        
        T=self.time
        output=[]
	reference_cn=dict()
	J=0
	for f in files:
		
		with open(f, 'rb') as csvfile:
		    read = csv.reader(csvfile, delimiter=delim, quotechar='"')
		    first=1
		    J+=1
		    for row in read:
			if first==1:
				row0=row
		
			#if row[0] not in output:
				
			output.append(row[1])
					
			if first==1:
				first=0
		   
        
	output=list(set(output))

        return output
    
    def raw2fraction(self):
        
	T=self.time
        #print self.idata
        #outdata=numpy.array(self.idata)
        outdata=self.idata
	sF=self.sF

	for t in T:
	#for k in ['1','2','3']:
		norm=0
		for i in sF:
			if 't3_' in i:
				norm=np.mean(sF[i].values())
			else:
				norm=sF[i]

			
			if i[2]=='_':
					
				for j in outdata:
					if type(outdata[j][i[1:2:1]][i[3]])!=type(0):
						outdata[j][i[1:2:1]][i[3]]=0
					else:
						outdata[j][i[1:2:1]][i[3]]=outdata[j][i[1:2:1]][i[3]]/norm
			elif i[3]=='_':
				for j in outdata:
					
					if type(outdata[j][str(i[1:3:1])][str(i[4])])!=type(0):	
						outdata[j][str(i[1:3:1])][str(i[4])]=0	
					else:	
						outdata[j][str(i[1:3:1])][str(i[4])]=outdata[j][str(i[1:3:1])][str(i[4])]/norm


	
			elif i[4]=='_':
				for j in outdata:
					
					if type(outdata[j][str(i[1:4:1])][str(i[5])])!=type(0):

					#print i,i[0:3:1],str(i[4])
						outdata[j][str(i[1:4:1])][str(i[5])]=0
						
					else:
						outdata[j][str(i[1:4:1])][str(i[5])]=outdata[j][str(i[1:4:1])][str(i[5])]/norm
						
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
	
            

def run(prefix,datatype):
	#load the data
	t=[3,4,5,6,8,24,48,168,336]
	#create a dictionary to map the sample number back to the time point.
	sample_ref=dict()
	for i in range(len(t)):
		sample_ref[str(i)]=t[i]

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

	out=ngram()
	#load the raw spectral counts
	#for j in t:
	    #for i in sample_ref_r[str(j)]:
	#prefix='/media/HD1_/Documents/AG3C_data/data/experiments/time_course_P0/'
        files=prefix+datatype+'_data.csv'
	#print files
	out=D.load_mult_csv(files,'\t')
	    	
	
	#print out
	D.idata=out
	
	files=[]
	for j in t:
		if j!=3:
			files.append(prefix+'size_factorst3_t'+str(j))
		
	sF=D.load_size_factors(files,',')

	D.sF=sF

	files=[]
	for j in t:
		if j!=3:
			files.append(prefix+'sig_gene_t3_t'+str(j))
		
	
	pval=D.load_p_values(files,',')

	D.pval=pval

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
			#print D.avout[i]
			temp_av=[]
			temp_st=[]
			for k in t:
				temp_av.append(D.avout[i]['t'+str(k)])
				temp_st.append(D.stdout[i]['t'+str(k)])

			avout.append(temp_av)
			stdout.append(temp_st)
		
	

	pickle.dump((avout,stdout,refout),open(prefix+datatype+'_data_format.p','wb'))



if __name__ == "__main__":

	run()
