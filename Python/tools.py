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
    G = np.zeros((r,2**r-1))
    for i in range(int(2**r-1)):
        # Creates the B.F. columns
        if np.log2(i+1).is_integer():
            G[int(np.log2(i+1)),i] = 1
        # Creates the interactions columns as sum of the B.F. modulo 2
        else:
            a = (i+1)-int(2**int(np.log2(i+1)))
            b = (i+1)-a
            G[:,i] = G[:,a-1] + G[:,b-1] 
    return G.astype(int)

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
    
def DOPselect(D):
    """
    Generates columns of candidate designs according to the DOP selection criterion.

    Parameters
    ----------
    D : Design() 
        Regular two-level design.

    Returns
    -------
    out : list
        List of lists of columns numbers, in 1-based indexing, of the selected designs.

    """
    B = Bmat(D.r)
    candicols = [i for i in range(1,D.N) if i not in D.cols]
    out = []
    for col in candicols:
        tempcols = D.cols + [col]
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

def STselect(D,order="rL"):
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
        ordint = list(range(1,D.N))
    elif order == "cL":
        ordint = [sum(j) for j in sorted([__pow2Fac(i) for i in range(1,16)],key=lambda x: (len(x),x))]
    else:
        ValueError("Unknown order")
    candicols = [i for i in range(1,D.N) if i not in D.cols]
    return [D.cols+[col] for col in candicols if ordint.index(col)>ordint.index(D.cols[-1])]

    
    

#%% Partition


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
def __in_rLform(G):
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
    # Find the rL order of G
    Gorder = [int(''.join([str(i) for i in G[:,i]][::-1]),2) for i in range(G.shape[1])]
    # re-order G in rL order
    G =  G[:,np.argsort(Gorder)]
    # Find index of the B.F. in G
    LrefAFindex = [i for i in range(G.shape[1]) if sum(G[:,i]) >1]
    # Reference L matrix
    Lref = G[:,LrefAFindex]
    # All B.F. set possibilities
    rposs = list(combinations(range(G.shape[1]),G.shape[0]))
    # All row permutations 
    perm = list(permutations(range(G.shape[0])))
    for r in rposs:
        # # Reject set where last added factor is not included
        # if lastfac-1 not in r:
        #     continue
        # Define new set of B.F.
        R = G[:,r]
        # Check if R is singular
        if np.linalg.det(R) == 0:
            continue 
        # Compute new interactions matrix - L
        K = G[:,[i for i in range(G.shape[1]) if i not in r]]
        L = np.linalg.solve(R,K).astype(int)%2
        # Permute the rows
        for p in perm:
            Lstar = L[p,:]
            # Re-arrange L in rL order
            Lindex = [int(''.join([str(i) for i in Lstar[:,i]][::-1]),2) for i in range(Lstar.shape[1])]
            Lstar = Lstar[:,np.argsort(Lindex)]
            
            # Test if L is rL-smaller than Lref
            if list(Lstar.flatten('F')[::-1]) < list(Lref.flatten('F')[::-1]):
                print(L,r,p)
                return False
    return True

def rLiso(C,r):
    G = Gmat(r)
    al = [G[:,[i-1 for i in cols]] for cols in C]
    return [C[i] for i in range(len(al)) if __in_rLform(al[i])]

#%% Design class
class Design:
    def __init__(self,N,cols):
        self.N = N
        self.cols = cols
        self.n = len(cols)
        self.r = int(np.log2(self.N))
        self.k = self.n - self.r
    
    def __repr__(self):
        return 'Design: %i runs - cols = %s'%(self.N,str(self.cols))