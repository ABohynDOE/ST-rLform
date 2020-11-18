#%% Packages
import tools
import pandas as pd
import numpy as np

# #%% ST-rLform Catalog generation
# ### Input parameters
# N = 16;
# verbose = True;

# ### Initiate catalog
# df = pd.DataFrame()
# r = int(np.log2(N))
# bf = [int(2**i) for i in range(r)]
# startcols = bf

# ### Start iterating
# for n in range(r+1,15):
#     if verbose:
#         print('%i factors\n'%n)
#     k = n-r
#     outcols = []
#     # First iteration case
#     if n == r+1:
#         candicols = [int(2**i)-1 for i in range(r+1) if i >1 ]
#         outcols = [startcols + [i] for i in candicols]
        
#     else:
#         # Following iteration
#         for cols in startcols:
#             candicols = tools.STselect(N,cols)
#             if verbose:
#                     print('\t%i candidates\n'%len(candicols))
#             for d in candicols:
#                 if tools.__rLform(N,d):
#                     outcols.append(d)
                  
#     if verbose:
#         print('\t%i unique designs\n'%len(outcols))
#     # Input in df 
#     for d in outcols:
#         df = df.append({'n':n,'k':k,'cols':d},ignore_index=True)

#     startcols = outcols

#%% NAUTY cattalog generation
### Input parameters
N = 32;
verbose = True;

### Initiate catalog
df = pd.DataFrame()
r = int(np.log2(N))
bf = [int(2**i) for i in range(r)]
startcols = bf
candiMat = pd.DataFrame()

### Start iterating
for n in range(r+1,16):
    if verbose:
        print('%i factors\n'%n)
    k = n-r
    outcols = []
    # First iteration case
    if n == r+1:
        candicols = [int(2**i)-1 for i in range(r+1) if i >1 ]
        outcols = [startcols + [i] for i in candicols]
    
    else:
        candicols = []
        for cols in startcols:
            tempcandi = tools.DOPselect(N,cols)
            for d in tempcandi:
                candicols.append(d)
        if verbose:
                print('\t%i candidates\n'%len(candicols))
        outcols = tools.NAUTYiso(candicols,r)
    
    if verbose:
        print('\t%i unique designs\n'%len(outcols))
    candiMat = candiMat.append({'n' : n, 'candidates' : len(candicols), 'unique' : len(outcols)}, ignore_index=True)
    # Input in df 
    ind = 0
    for d in outcols:
        ind += 1
        df = df.append({'n':n,'k':k,'i' : ind, 'cols':d},ignore_index=True)

    startcols = outcols
