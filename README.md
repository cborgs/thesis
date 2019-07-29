# Thesis code


All code needed to reproduce the thesis results for
Borgs, C. (2019): Optimal parameter choice for Bloom filter-based Privacy-preserving Record Linkage. Thesis. University of Duisburg-Essen.


Steps:
1. Unpack Diss_CB_Source.tar.gz OR use the R command packrat::unbundle(file = "Diss_CB_Source.tar.gz") (requires packrat to be loaded)
2. Open 00_Main.R using Rstudio
3. Install packrat package from CRAN if not already installed and run packrat::packify()
NOTE: Packrat will now use the repository in the tar.gz. with all packages as used in this thesis 
4. Run the source files in the order given by the source file
NOTE: If only the plots and tables need reproduction, run only the last source file. The required files are already present

