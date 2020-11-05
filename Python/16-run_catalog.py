#%% Packages
import tools
import pandas as pd
import numpy as np

#%% Catalog generation
### Input parameters
N = 16;
verbose = True;

### Initiate catalog
df = pd.DataFrame()
r = int(np.log2(N))
bf = [int(2**i) for i in range(r)]
startcols = bf

### Start iterating
for n in range(r+1,14):
    k = n-r
    outcols = []
    # First iteration case
    if n == r+1:
        candicols = [int(2**i)-1 for i in range(r+1) if i >1 ]
        outcols = [startcols + [i] for i in candicols]
        
    else:
        # Following iteration
        for cols in startcols:
            candicols = tools.STselect(N,cols)
            for d in candicols:
                if tools.__rLform(N,d):
                    outcols.append(d)
                    
    # Input in df 
    for d in outcols:
        df = df.append({'n':n,'k':k,'cols':d},ignore_index=True)

    startcols = outcols
