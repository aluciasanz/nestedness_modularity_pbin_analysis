# Analysis of modularity and nestedness using bimat
Scripts **_NODF_MOD_dynamics.m_** and **_LPBRIM_analysis.m_** were used to analyze nestedness and modularity of $\phi 21$ phage and _E. coli_ bacteria coevolutionary infection dynamics for the publication [nameofpublication]
# (REQUIRED) Install lastest version of <a href="https://bimat.github.io/">Bimat</a> from source
Follow carefuly the installation guidelines https://bimat.github.io/inst/getting_started.html#1
Note that bimat is required to run the next analysis.
# Run *.m files
### 1. If you haven't add bimat source folder to matlab path sources, add the source in line 5  
```matlab
PATH='~/path_to_bimat_folder'; 
GEN_PATH = genpath(PATH);
addpath(GEN_PATH);
mkdir(strcat(PATH,'/myfolder')) % create a folder to save data   
cd(strcat(PATH,'/myfolder')); % use the new folder as the working directory
```
### 2. Set the directory to the matrix you want to analyze in line 12.
```matlab
matrix_raw=readtable(strcat(cd,'/Data/pbin_to_analyze.csv')); % PBIN located in Data
```
Here the matrix data is contained in the folder *Data*
