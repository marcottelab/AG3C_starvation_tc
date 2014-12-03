# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import copy
import pickle
    
class find_param(object):
    def __init__(self,data=None,method='DE'):
        
        self.method,self.data=method,data
        self.maxit=10
        self.Nagents=10
        self.t0=0
        self.CR=0.75
        self.F=0.6
        self.T=0.5
        self.k=dict()
        
        for i in range(3):
            self.k['A'+str(i+1)]=0
        for i in range(5):
            self.k['tau'+str(i+1)]=0
        
    def main(self):
        
        if self.method=='DE':
            self.DE_run()
        elif self.method=='MCMC':
            self.MCMC_run()
        else:
            print "Error: I do not have a parameter search method called "+self.method
    
    def Objective(self,Y):
        #unpack the data
        td,data,std=self.data
        ts,sim=Y

        keys=data.keys()

        Err=0.
        
        d=data['x']
	sim=np.array(sim)/np.mean(sim)
	d=np.array(d)/np.mean(d)

        for j in range(len(d)):
            if d[j]==0:
                d[j]=0.5
                std[j]=0.5
	    if std[j]==0:
		std[j]=0.1*d[j]
                
            Err=Err+((d[j]-sim[j])**2)/(std[j]**2)

        
        return Err
    
    def func(self,k,data):
        
        td,data,std=data
        
        T1=k['tau1']+k['tau2']
        T2=k['tau1']+k['tau2']+k['tau3']
        F=np.zeros_like(td)
        F=np.float32(F)
        td=np.float32(td)
        
        for i in range(len(td)):
            
            if td[i]<k['tau1']:
                F[i]=k['A1']
                
            if td[i]>=k['tau1'] and td[i]<(T1):
                F[i]=k['A1']+(td[i]-k['tau1'])*(k['A2']-k['A1'])/(k['tau2'])
            
            if td[i]>=(T1) and td[i]<(T2):
                
                F[i]=k['A2']
                
            if td[i]>=(T2) and td[i]<(T2+k['tau4']):
                F[i]=k['A2']+(td[i]-T2)*(k['A3']-k['A2'])/(k['tau4'])
            
            if td[i]>=(T2+k['tau4']):
                F[i]=k['A3']
            
            #else:
            #    F.append(k['A3'])
        
        Y=td,F
        return Y 
        #np.linspace()
        
        
    def DE_run(self):
        
        accepted=0
        Agents0=[]
        Err0=[]
        k=self.k
        k_keys=k.keys()
        Errout=[]
        #print k_keys
        
        #make some random initial guesses for parameters
        for i in range(self.Nagents):    
            for j in range(len(k_keys)):
                if 'A' in k_keys[j]:
                    k[k_keys[j]]=np.random.uniform(low=0.01,high=1)
                    
                else:
                    k[k_keys[j]]=2+3*np.random.lognormal(mean=0,sigma=0.5)
           
            Agents0.append(copy.deepcopy(k))
            #Solve the equations using the initial parameter guess
            
            #get the output
            Y0=self.func(k,self.data)
            #calculate the Error function
            Err0.append(self.Objective(Y0))
        Agents=copy.deepcopy(Agents0)
        
        for i in range(self.maxit):
            
            for j in range(self.Nagents):
                agent_list=[]
                agent_list=np.array(range(self.Nagents))
                np.delete(agent_list,j)
                

                a=[]
                for l in range(2):  
                    idx=np.random.randint(0,high=len(agent_list)-1)
                    a.append(Agents0[agent_list[idx]])
                    np.delete(agent_list,idx)

                #pick a random parameter that will be changed no matter what
                ridx=np.random.randint(0,len(k_keys))
                tempa=[]
                for l in range(len(k_keys)):
                    
                    
                    randomn=np.random.uniform()
                    if((ridx==l) or (randomn<self.CR)):
                        
                        if 'A' in k_keys[l]:
                            b=0.3
                        else:
                            b=100
                            
                        e=0
                        #np.random.uniform(-b,b)
                        tempa.append(np.abs(Agents0[j][k_keys[l]]+self.F*(a[0][k_keys[l]]-a[1][k_keys[l]])+e))
                    
                        if 'tau1' in k_keys[l]:
                            if tempa[l]<3:
                                tempa[l]=3.5+np.random.uniform(-0.4,0.4)
                        if 'tau' in k_keys[l]:
                            if tempa[l]<0.5:
                                tempa[l]=3.5+np.random.uniform(-0.4,0.4)
                        if 'A' in k_keys[l]:
                            if tempa[l]>1:
                                tempa[l]=0.95+np.random.uniform(-0.4,0.4)
                                
                    else:
                        tempa.append(Agents0[j][k_keys[l]])
                        
                        
                    Agents[j][k_keys[l]]=copy.copy(tempa[l])
                
            #print tempa
            
            #compute the fitness of the evolved agents
            Err=[]
            Errtemp=[]
            tempA=[]
            Yout=[]
            
            for j in range(self.Nagents):
                #print Agents[j]
                k=Agents[j]
                Y=self.func(k,self.data)
                Yout.append(Y)
                #calculate the Error function
                Err.append(self.Objective(Y))
                
                r=np.random.uniform()
                dE=np.exp((Err0[j]-Err[j])/self.T)
                if Err[j]<Err0[j]:# or (r<dE):
                    accepted+=1
                    tempA.append(Agents[j])
                    Errtemp.append(Err[j])
                else:
                    tempA.append(Agents0[j])
                    Errtemp.append(Err0[j])
            
            for j in range(self.Nagents):
                Err0[j]=copy.deepcopy(Errtemp[j])
                
                Agents0[j]=copy.deepcopy(tempA[j])
                
            Errout.append(Errtemp)
            
        self.Yout=Yout    
        self.Errout=Errout
        self.Agents_final=Agents
        self.accept=accepted
        print 'accepted='+str(accepted)
    

import time
def run(prefix,filename,time=None):
	if time==None:
		td=[3,4,5,6,8,24,48,1*7*24,2*7*24]
	
	else:
		td=time

	(avout,stdout,av_reference)=pickle.load(open(prefix+filename+'_results/'+filename+'_data_format.p','rb'))

	k_out_av=['']*len(avout)
	k_out_std=['']*len(avout)
	ftype=[]
	yav_out=[]
	ystd_out=[]	
	data_out=[]
	Errout=[]
	#ta=time.time()

	for i in range(len(avout)):
		
	    d=dict()
	    d['x']=avout[i]
	  
	    #td=[3,4,5,6]
	    data=td,d,stdout[i]
	    
	    data_out.append(data)
	    
	    search=find_param(data)
	    search.Nagents=12
	    search.maxit=150
	    search.main()
	    keys=search.k.keys()
	    idx=np.where(np.array(search.Errout[search.maxit-1])<=np.median(np.array(search.Errout[search.maxit-1])))
	    yc=[]
	    paramsc=[]
	    
	    for j in range(len(idx[0])):
		tout,yout=search.Yout[idx[0][j]]
		yc.append(yout)
	    
	    for j in range(len(keys)):
		params=[]
		for k in range(len(idx[0])):
		    
		    params.append(search.Agents_final[idx[0][k]][keys[j]])
		
		paramsc.append(params)
		   
	    yav=np.average(yc,axis=0)
	    ystd=np.std(yc,axis=0)
	    pav=np.average(paramsc,axis=1)
	    pstd=np.std(paramsc,axis=1)
	    yav_out.append(yav)
	    ystd_out.append(ystd)

	    Errout.append(search.Errout)
	    k=dict()
	    kstd=dict()
	    for j in range(len(keys)):
		k[keys[j]]=pav[j]
		kstd[keys[j]]=pstd[j]
	    k_out_av[i]=k
	    k_out_std[i]=kstd

	#tb=time.time()
	#print tb-ta
	results = (yav_out,ystd_out,k_out_av,k_out_std,av_reference,Errout)
	pickle.dump(results, open(prefix+filename+'_results/'+filename+'_fits.p', "wb" ) )
if __name__=="__main__":
	run()

