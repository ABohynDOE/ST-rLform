#%% Packages
import numpy as np
import oapackage as oa 
from itertools import combinations,permutations
from HadamardTools import selectIsomorphismClasses
import sys

#%% Matrix
### Basic factor matrix
def Rmat(r):
    """ Creates the basic factors matrix (N x r), for r basic factors."""
    b = Gmat(4)
    return np.flip(b.T,axis=1)

### Reduced generalized interaction matrix
def Gmat(r):
    """ Creates the reduced generalized interaction matrix G (r x N-1), for r basic factors."""
    a = np.array(range(1,int(2**(r))),dtype=np.uint8)
    a = a[...,None]
    b = np.unpackbits(a.T,axis=0,bitorder='little',count=r)
    return np.array(b,dtype=int)

### Full generalized interaction matrix
def Bmat(r):
    """ Creates the generalized interaction matrix B (N x N-1), for r basic factors."""
    return (np.matmul(Rmat(r),Gmat(r))%2).astype(int)

#%% Selection
### DOP selection
def __isMA(wlplst):
    """ Compare the WLP of a list"""
    for ii in range(len(wlplst)-1):
        if wlplst[ii] < wlplst[-1]:
            return False
    return True
    
def DOPselect(N,cols):
    """
    Generates columns of candidate designs according to the DOP selection criterion.

    Parameters
    ----------
    cols : list
        List of columns numbers.

    Returns
    -------
    out : list
        List of lists of columns numbers, in 1-based indexing, of the selected designs.

    """
    r = int(np.log2(N))
    B = Bmat(r)
    candicols = [i for i in range(1,N) if i not in cols]
    out = []
    for col in candicols:
        tempcols = cols + [col]
        Dc = B[:,[i-1 for i in tempcols]]
        wlplst = []
        for j in range(Dc.shape[1]):
            DOP = np.delete(Dc,j,axis=1)
            wlplst.append(oa.array_link(DOP).GWLP()[3:])
        if __isMA(wlplst):
            out.append(tempcols)    
    return out

### ST selection
def __pow2Fac(x,char=False):
    """ Decompose a number into the smallest amount of power of 2"""
    powers = []
    i = 1
    while i <= x:
        if i & x:
            powers.append(i)
        i <<= 1
    return powers

def STselect(N,cols,order="rL"):
    """
     Generates columns of candidate designs according to the search-table selection method, using the order specified.

    Parameters
    ----------
    D : Design
        Regular two-level design.
    order : {'rL','cL'}, optional
        Type of order to use in the search-table. The default is "rL". 
        'rL' is for reverserd lexicographic and 'cL' is for conditional lexicographic ordering.

    Returns
    -------
    list
        List of lists of columns numbers, in 1-based indexing, of the selected designs.

    """
    if order == "rL":
        ordint = list(range(1,N))
    elif order == "cL":
        ordint = [sum(j) for j in sorted([__pow2Fac(i) for i in range(1,N)],key=lambda x: (len(x),x))]
    else:
        ValueError("Unknown order")
    candicols = [i for i in range(1,N) if i not in cols]
    return [cols+[col] for col in candicols if ordint.index(col)>ordint.index(cols[-1])]



#%% Isomorphism
### NAUTY isomorphism
def NAUTYiso(C,r):
    """
    Isomorphism class representative selection using NAUTY algorithm.

    Parameters
    ----------
    C : list
        List of lists of column numbers, in 1-based indexing, of the conditate designs
    r : int
        Number of basic factors in the candidate designs.

    Returns
    -------
    list
        List of lists of columns numbers, in 1-based indexing, of the design selected as
        representatives of their isomorphism classes.

    """
    B = Bmat(r)
    # Converts columns to designs in OA
    al = [oa.array_link(B[:,[i-1 for i in x]]) for x in C]
    # Define the isomorphism classes
    ind,isoClass = selectIsomorphismClasses(al, verbose=0)
    # Select one rep. per class
    vals,zz = np.unique(ind, return_index=True)
    zz.sort()
    return [C[x] for x in list(zz)]

### LMC isomorphism
def bi2de(binary):
    # Same basic function as matlab bi2de function
    # Convert a binary array to a decimal number
    # Read from right to left in binary array
    bin_temp = 0
    bin_res = np.zeros(len(binary), dtype=int)
    for i in range(len(binary)):
        for j in range(len(binary[i])):
            bin_temp = bin_temp + binary[i][j] * (2 ** j)
        bin_res[i] = bin_temp
        bin_temp = 0
    return bin_res


def rLsmaller(L,Ls):
    a = [i for i in bi2de(L.T) if i not in bi2de(Ls.T)]
    if a == []:
        return False
    b = [i for i in bi2de(Ls.T) if i not in bi2de(L.T)]
    return a[0] < b[0]    

def __rLform(N,cols):
    """
    Checks if a matrix is in rL-form

    Parameters
    ----------
    G : array
        Design matrix in reduced design matrix form (r x n), with columns in first-to-last-added order.
    lastfac : int
        Index of the last added factor, in the input array.

    Returns
    -------
    bool
        Returns True if the array is in rL-form.

    """
    # Reduced design matrix of the columns
    r = int(np.log2(N))
    a = np.array(sorted(cols),dtype=np.uint8)
    a = a[...,None]
    b = np.unpackbits(a.T,axis=0,bitorder='little',count=r)
    S = np.array(b,dtype=int)
    lastfac  = cols[-1]
    
    # Reference L matrix
    Lref =  S[:,S.sum(0)> 1]
    # All B.F. set possibilities
    rposs = list(combinations(range(S.shape[1]),r))
    # All row permutations 
    perm = list(permutations(range(r)))
    for r in rposs:
        # # Reject set where last added factor is not included
        # if lastfac-1 not in r:
        #     continue
        # Define new set of B.F.
        R = S[:,r]
        # Check if R is singular
        if np.linalg.cond(R) >= 1/sys.float_info.epsilon:
            continue 
        # Compute new interactions matrix - L
        K = S[:,[i for i in range(S.shape[1]) if i not in r]]
        L = np.linalg.solve(R,K).astype(int)%2
        # Check that L is the same rank as K
        if np.linalg.matrix_rank(L) != np.linalg.matrix_rank(K):
            continue
        # Permute the rows
        for p in perm:
            Lstar = L[p,:]
            # Re-arrange L in rL order
            Lstar = Lstar[:,np.argsort(bi2de(Lstar.T))]
            # Test if L is rL-smaller than Lref
            if rLsmaller(Lstar,Lref):
                print(R,L)
                return False
    return True

def rLiso(N,C):
    return [col for col in C if __rLform(N, col)]
