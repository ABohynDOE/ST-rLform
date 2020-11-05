#%% Packages
import tools
from itertools import chain

#%% Catalog generation
### Starting designs
Dlist = [tools.Design(16,[1,2,4,8,i]) for i in [3,7,15]]

### ST selection of the candidates
candilist = [tools.STselect(x) for x in Dlist]
candilist = list(chain(*candilist))

### Isomorphism reduction through rL-form
tools.rLiso()
