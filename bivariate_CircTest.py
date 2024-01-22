"""
Created on Sat Oct 24 10:56:35 2020

@author: mahsaaltafi
"""

def bivariate_CircTest(d1, d2, m1, m2, iter):
    
    # A homogeneity test using likelihood ratio test to compare two bivariate circular distributions
    # H0: The two bivariate distributions are not different
    #
    # Inputs:
    #        d1 :a 2d numpy array containing a bivariate distribution in rad (0, 2pi]
    #        d2 :a 2d numpy array containing another bivariate distribution in rad (0, 2pi]
    #        m1 :maximum number of modes to assume for first variable
    #        m2 :maximum number of modes to assume for 2nd variable
    #       iter  :the maximum number of iteration for estimation (usually set to 1000 or more)
    #
    # Output: 
    #       DataFrame   :containing the loglikelihhod for each assumption of number of modes, Chi^2 p-value, df etc.  
    #
    #Reference:
    #Fernández-Durán, J. J. & Domínguez, M. M. G. A Likelihood Ratio Test for   
    #Homogeneity in Circular Data. J. Biom. Biostat. 1, 1–9 (2010).
    
    from scipy.io import loadmat
    import numpy as np
    import pandas as pd
    from scipy import stats
    import sys

    d= np.concatenate([d1,d2], axis=0)
    
    DataFrame= pd.DataFrame([])
    k=0
    for i in range(m1+1):
        for j in range(m2+1):
            R1=multi_var_estim(d1, [i,j], iter)   #loglikelihood estimation of 1st bivariate distribution
            R2=multi_var_estim(d2, [i,j], iter)  #loglikelihood estimation of 2nd bivariate distribution
            R=multi_var_estim(d, [i,j], iter)   #loglikelihood estimation of the combined distributions
        
            df=2*(i+1)*(j+1)-2
            likHom= -2*(R-(R1+R2))   #likelihood ratio
            pval=1-stats.chi2.cdf(likHom,df)  #p-value of Chi^2 test
        
            DataFrame=DataFrame.append(pd.DataFrame({'M1':i, 'M2':j, 'loglik1':R1, 'loglik2':R2, 'loglik':R,\
                                                 'Homogenity test statistics':likHom, 'df':df, 'chi^2 p value':pval}, index=[k]))
        
            k +=1
    return DataFrame    
        


def multi_var_estim(data, M, iter):
    
    # Computes the loglikelihood estimation of a circular bivariate distribution having maximum (M[0], M[1]) number of modes 
    # Inputs:
    #        data :a 2d numpy array containing the bivariate distribution in rad (0, 2pi]
    #        M    :a 2d list containig the number of modes for each variable 
    #       iter  :the maximum number of iteration for estimation (usually set to 1000 or more)
    #
    #Output: 
    #       loglik   :the loglikelihood estimation 
    
    

    import numpy as np 
    import pandas as pd
    import math
    
    data = np.asarray(data)
    n = np.size(data,0)
    dim = np.size(M)
    
    
    sec=[]
    
    if dim==1:
        sec=list(np.arange(M[0]+1))
    elif dim==2:
        sec=[(x0, x1) for x0 in range(M[0]+1) for x1 in range(M[1]+1)]    #list of all possible number of modes which should be changed for dim>=3
    else: 
        return "Incorrect number of variables (for dimensions above 2, the code needs small adjusments)"
    ind = np.asarray(sec)
    
  
    statmat= np.zeros([np.prod([x+1 for x in M]),n])
    
    statmat= np.exp(-(np.matmul(ind, np.transpose(data)))*complex(0,1))  
    
    
    c0=np.mean(statmat, axis=1)  #initial point can be adjusted
    c0=c0/math.sqrt(np.sum(np.abs(c0)**2))
    c0= np.reshape(c0, (-1,1))
    
    eta= np.zeros([np.prod([x+1 for x in M]), 1])
    
    for k in range(n):
        eta += np.concatenate((1/n)* 1/np.matmul(np.transpose(np.conj(c0)), np.reshape(statmat[:,k], (-1,1)))) * np.reshape(statmat[:,k], (-1,1))
    
    eta -= c0;
    newtonmanifold = c0 + eta
    newtonmanifold = newtonmanifold/math.sqrt(np.sum(np.abs(newtonmanifold)**2))
    newtonmanifold = newtonmanifold * np.exp(-(complex(0,1)) * np.angle(newtonmanifold[0,0]))
    
    newtonmanifoldold = newtonmanifold
    for j in range(iter):
        eta = np.zeros([np.prod([x+1 for x in M]), 1])
        for k in range(n):
            eta += np.concatenate((1/n) * 1/np.matmul(np.transpose(np.conj(newtonmanifold)), np.reshape(statmat[:,k], (-1,1)))) * np.reshape(statmat[:,k], (-1,1))
        eta -= newtonmanifold
        newtonmanifold = newtonmanifold +eta
        newtonmanifold = newtonmanifold/math.sqrt(np.sum(np.abs(newtonmanifold)**2))
        newtonmanifold = newtonmanifold * np.exp(-(complex(0,1)) * np.angle(newtonmanifold[0,0]))
        
        if j==(iter-1):
            err= math.sqrt(np.sum(np.abs(newtonmanifold - newtonmanifoldold)**2))
            newtonmanifoldold = newtonmanifold
            
    newtonmanifold = newtonmanifold / math.sqrt(2*math.pi)**dim 
    
    if np.sum(M)==0:
       return (-n * dim * math.log(2*math.pi))
    statmats= np.exp(np.matmul(ind, np.transpose(data))*complex(0,1))  
    aux= np.matmul(np.transpose(np.asarray(newtonmanifold)), statmats)
    res = np.real(aux*np.conj(aux))
    res = np.transpose(res)
    loglik = np.sum(np.log(res.flatten()))
    return loglik












