# -*- coding: utf-8 -*-
### Packages
import numpy as np
import tools 
import oapackage as oa
import time 
import pandas as pd
from itertools import chain 
from HadamardTools import selectIsomorphismClasses

### Basic infos 
N = 16
maxfac = 15
verbose = True

### Initiate algorithm
def CatInit(N,maxfac):
    global r
    r = int(np.log2(N))
    bf = [2**i for i in range(r)]
    firstGenCandi = [i-1 for i in [2**j for j in range(1,r+1)] if i>2]
    global startcols
    startcols = [bf+[i] for i in firstGenCandi]
    global B
    B = tools.Bmat(r)
    global G 
    G = tools.Gmat(r)
    global cat
    cat = pd.DataFrame()
    for cols in startcols:
        cat = cat.append({'n':r+1, 'cols' : cols},ignore_index=True)


#%% rLmin procedure
CatInit(N,maxfac)
print('---rLmin procedure\t','-'*50)
print('%i runs - up to %i factors\n'%(N,maxfac))
firstTic = time.time()
for k in range(r+1,maxfac):
    tic = time.time()
    # Get candidates
    nbrCandi = 0
    isocols = []
    for cols in startcols:
        candicols = tools.STselect(N,cols)
        nbrCandi += len(candicols)
        
        # rLmin form checking
        for cols in candicols:
            candiReducedMat = G[:,[i-1 for i in cols]]
            if tools.rLmin(candiReducedMat):
                isocols.append(cols)
                cat = cat.append({'n':k+1,'cols':cols},ignore_index=True)
    toc = time.time()
    # Initiate next run
    if verbose:
        print('%i factors: %i candidates - %i representatives\t\t (%.2f sec.)'%(k+1,nbrCandi,len(isocols),tic-toc))
    startcols = isocols
    
totalTime = time.time() - firstTic
print('\nrLmin procedure: Catalog generated in %.2f seconds'%totalTime)
print('-'*70,'\n')
       


#%% NAUTY procedure
def NAUTYiso(al,candicols):
    # Define the isomorphism classes
    ind,isoClass = selectIsomorphismClasses(al, verbose=0)
    # Select one rep. per class
    vals,zz = np.unique(ind, return_index=True)
    zz.sort()
    return [candicols[x] for x in list(zz)]

CatInit(N,maxfac)
print('---NAUTY procedure\t','-'*50)
print('%i runs - up to %i factors\n'%(N,maxfac))
firstTic = time.time()
for k in range(r+1,maxfac):
    tic = time.time()
    # Get candidates
    candicols = []
    for cols in startcols:
        candicols.append(tools.DOPselect(N,cols))
    candicols = list(chain(*candicols))
       
    # Isomorphism reduction
    al = []
    for cols in candicols:
        candiD = B[:,[i-1 for i in sorted(cols)]]
        al.append(oa.array_link(candiD))   
    isocols = NAUTYiso(al,candicols)
    for cols in isocols:    
        cat = cat.append({'n':k+1,'cols':cols},ignore_index=True)
    toc = time.time()
    # Initiate next run
    if verbose:
        print('%i factors: %i candidates - %i representatives\t\t (%.2f sec.)'%(k+1,nbrCandi,len(isocols),tic-toc))
    startcols = isocols
    
totalTime = time.time() - firstTic
print('\nNAUTY procedure: Catalog generated in %.2f seconds'%totalTime)
print('-'*70,'\n')
