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

class format_data(object):
    def __init__(self,delim='\t'):
        #self.files=files
        self.delim=delim
        self.idata=[]
        self.out=[]
        self.minchange=1.5
        self.cutoff=7
   
    def main(self):
        
        output=self.idata
        
        #strip the header row
        for i in range(len(output)):
            output[i][0]=[0,0,0]
        #convert to int
        for i in range(len(output)):
            for j in range(len(output[i])):
                for k in range(3):
                    try:
                        output[i][j][k]=float(output[i][j][k])
                    except:
                        output[i][j][k]=0.0
        
        self.idata=output
        
        #convert from raw data to fractions
        self.raw2fraction()
        
        #get rid of peptides that are not observed 2/3 of the time
        self.cleandata()
        output2=self.output
        reference2=self.refout
        av_out=[]
        std_out=[]
        n_out=[]
        exitflag=0
        #average the output data
        for j in range(len(output2[0])):
            temp=[]
            temp2=[]
            temps=[]
            for i in range(len(output2)):
                N=0.
                o=0.
                s=0.
                for k in range(len(output2[i][j])):
                    if output2[i][j][k]!=0.:
                        o=o+output2[i][j][k]
                        N+=1.
                 
                if N!=0.:
                    o=o/N
                    no=o
                else:
                    o=0.
                    no=0.
                    
                for k in range(len(output2[i][j])):
                    if output2[i][j][k]!=0.:    
                        s=np.add(s,(np.subtract(output2[i][j][k],o))**2)
                
                if N!=0.:
                    s=numpy.sqrt(s/float(N))
                else:
                    s=0.
                    
                temp.append(o)
                temp2.append(no)
                temps.append(s)
                
            av_out.append(temp)
            std_out.append(temps)
            n_out.append(temp2)
            
        av_out2=[]
        std_out2=[]
        n_out2=[]
        ref2=[]
        #make a subtable where the fold change is >=1.5
        for i in range(len(av_out)):
            
            t=numpy.where(numpy.array(av_out[i])>0)
            idx=[]
            for j in range(len(t[0])):
                idx.append(av_out[i][int(t[0][j])])
            
            if not idx:
                mina=1
                maxa=1
            else:
                mina=min(idx)
                maxa=max(av_out[i])
            
            if maxa/mina>=self.minchange:
                av_out2.append(av_out[i])
                std_out2.append(std_out[i])
                n_out2.append(n_out[i])
                ref2.append(reference2[0][i])

        #normalize to initial value
        av=[]
        std=[]
        for i in range(len(av_out2)):
           temp=[]
           temp2=[]
           for j in range(len(av_out2[i])):
                try:
                    ma=max(av_out2[i])
                    if ma==0:
                        ma=1.
        
                    norma=av_out2[i][j]/ma
                    normb=std_out2[i][j]/ma
                    temp.append(norma)
                    temp2.append(normb)
                except:
                    temp.append(0.0)
                    temp2.append(0.0)
           
           av.append(temp)
           std.append(temp2)
         
           self.avout=np.array(av)
           self.std=std
           self.refout=ref2
           self.nout=np.array(n_out2)
                    
    def load_mult_csv(self,files,delim):
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

    def raw2fraction(self):
        
        #print self.idata
        #outdata=numpy.array(self.idata)
        outdata=self.idata
        #data[0][1][2]
        for i in range(len(outdata)):
            N=np.sum(outdata[i],0)
            #print(N)
            for k in [0,1,2]:
                if N[k]!=0:
                    for j in range(len(outdata[i])):
                        outdata[i][j][k]=outdata[i][j][k]/(N[k])
                     
        
        self.output=list(outdata)
    
    def cleandata(self):
    
        cutoff=self.cutoff
        data=self.output
        counts=self.counts
        ref=self.reference
        #require that the minimum number of counts is > a particular cutoff
        output=[]
        refout=[]
        
        for i in range(len(data)):
            temp=[]
            temp2=[]
            for j in range(len(data[i])):
                
                if not min(counts[j]):
                    a=0
                else:
                    a=min(counts[j])
                
                
                if int(a)>cutoff:
                    #try:
                    temp.append(data[i][j])
                        
                    if i==0:
                        temp2.append(ref[j])
                    #except:
                        #'problem with different data lengths'
            output.append(temp)
            if i==0:
                refout.append(temp2)
            
            self.output=output
            self.refout=refout
            
output=[]
if __name__ == "__main__":
	#load the data
	prefix='../data_prep/samples1-33/'
	files=[prefix+'3hr.csv',prefix+'4hr.csv',prefix+'5hr.csv',prefix+'6hr.csv',prefix+'8hr.csv',prefix+'24hr.csv',prefix+'48hr.csv',prefix+'1wk.csv',prefix+'2wk.csv']
	D=format_data()
	out=D.load_mult_csv(files,'\t')
	D.idata=out
	#output=load_mult_csv(files,'\t')

	#print len(output[0])
	#load the reference 
	files=['../data_prep/samples1-33/LabelsMaster.csv']
	reference=D.load_mult_csv(files,' ')
	D.reference=reference[0]

	counts=D.load_mult_csv(['../data_prep/SampleCounts.csv'],'\t')
	D.counts=counts[0]

	D.main()

	pickle.dump((D.avout,D.std,avref),open('./formated_data2b.p','wb'))

