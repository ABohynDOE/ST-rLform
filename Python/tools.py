#%% Packages
import numpy as np
import oapackage as oa 
from itertools import combinations,permutations
from HadamardTools import selectIsomorphismClasses

#%% Matrix
### Basic factor matrix
def Rmat(r):
    """ Creates the basic factors matrix (N x r), for r basic factors."""
    R = np.zeros((int(2**r),r))
    # Alternates 0 and 1 in columns until (0,1) is repeated N/2 times
    for i in range(r):
        R[:,i] = ([0]*(2**(r-(i+1)))+[1]*(2**(r-(i+1))))*(2**i)
    return R.astype(int)

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
    al = [oa.array_link(B[:,[i-1 for i in sorted(x)]]) for x in C]
    # Define the isomorphism classes
    ind,isoClass = selectIsomorphismClasses(al, verbose=0)
    # Select one rep. per class
    vals,zz = np.unique(ind, return_index=True)
    zz.sort()
    return [C[x] for x in list(zz)]


def col2num(L):
    r,n = L.shape
    rvec = np.array([2**i for i in range(r)])
    return rvec*L.T

def rLsmaller(L,Ls):
    C = L!=Ls
    colind = np.argmax(C.any(0))
    rowind = (C.shape[0]-1)-np.argmax(C[:,colind][::-1])
    if (colind,rowind) ==(0,0):
        return False
    return L[rowind,colind] < Ls[rowind,colind]

def rLmin(S):
    r = S.shape[0]
    n = S.shape[1]
    Lref = S[:,S.sum(0)>1];
    
    # All B.F. set possibilities
    rposs = list(combinations(range(n),r))
    # All row permutations 
    perm = list(permutations(range(r)))
    for r in rposs:
        # Define new set of B.F.
        R = S[:,r]
        # Check if R is singular
        if np.linalg.cond(R) >= 1/np.finfo(float).eps or np.linalg.det(R)==0:
            continue 
        # Compute new interactions matrix - L
        K = S[:,[i for i in range(n) if i not in r]]
        L = np.linalg.solve(R,K)%2
        # Check that L is binary
        if not np.array_equal(L, L.astype(bool)):
            continue
        # Permute the rows
        for p in perm:
            Lstar = L[p,:].astype(int)
            # Re-arrange L in rL order
            Lstar = Lstar[:,np.argsort(col2num(Lstar))]
            # Test if L is rL-smaller than Lref
            if rLsmaller(Lstar,Lref):
                return False
    return True

