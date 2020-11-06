# MATLAB files documentation

This files details each MATLAB function and file in the folder in terms of input, output and purpose.

## Functions
The main functions in the folder are:

  - **rLmin**: determines if a reduced design matrix is rL-minimal. 
     - *Input*: set of columns of the design tested.
     - *Output*: True if the design is rL-minimal.
  - **STselect**: selects candidate designs using the search-table method.
     - *Input*: set of columns of the parent design, order used in the search-table.
     - *Output*: List of columns' set of all the candidate designs.
  - **DOPselect**: selects candidate designs using the delete-one-factor projection method.
     - *Input*: set of columns of the parent design.
     - *Output*:  List of columns' set of all the candidate designs.

The rest of the functions are utility functions and are detailed in the MATLAB code.

## Files
There is a single MATLAB file in the folder: **Catalog_Test.m**.
It produces a catalog for N-run two-level regular designs up to n-factors.
There are 3 variables to change at the begining of the file:

  - **N**: the number of runs.
  - **n**: the maximal number of factors.
  - **verbose**: (boolean) to have info on the computation, and a summary table at the end.
