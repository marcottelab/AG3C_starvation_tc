# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

## Author: John Houser

import pickle
from scipy.cluster.vq import kmeans2
import numpy as np
import pylab as plt
import sys


def clust(Nclust=25,loadfile='./formated_data2.p',writefile='./cluster_kmeans.p',plots=False):

    
    #load the formated data
    (avout,stdout,av_reference)=pickle.load(open(loadfile,'rb'))
    
    TFav=[] 
    TFref=[] 
    TFstd=[] 
    
    #do k-means clustering
    (o,n)=kmeans2(np.array(avout),Nclust)
    
    #save the output
    pickle.dump((o,n),open(writefile,'wb'))

    if plots==True:
        #make the plots
        t=[3.,4.,5.,6.,8.,24.,48,7*24.,14*24.]
        plt.figure()
        for i in range(Nclust):
            
            idx=np.where(n==i)[0]
           
            plt.subplot(np.ceil(sqrt(Nclust)),np.ceil(sqrt(Nclust)),i+1)
            plt.hold(True)
            for j in range(len(idx)):
                plt.semilogx(t,avout[idx[j]])
                plt.ylim([0,1])
                
        plt.show()


clust()
