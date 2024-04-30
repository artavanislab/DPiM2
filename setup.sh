conda create -n DPiM
conda activate DPiM
conda install conda-forge::perl
conda install conda-forge::gcc
conda install conda-forge::gsl
conda install anaconda::boost
conda install bioconda::perl-bioperl
conda install bioconda::perl-datetime
conda install bioconda::perl-statistics-descriptive
conda install bioconda::perl-math-cdf
conda install bioconda::bioconductor-org.dm.eg.db
conda install bioconda::bioconductor-org.hs.eg.db
conda install conda-forge::r-biocmanager
R -q -e 'if(!require(data.table)) install.packages("data.table",repos = "http://cran.us.r-project.org")'
R -q -e 'if(!require(GO.db)) BiocManager::install("GO.db")'
R -q -e 'if(!require(reactome.db)) BiocManager::install("reactome.db")'

path=$(pwd)

DPSCR="${path}/src/DPiM/scripts/"; export DPSCR;
DPCPP="${path}/src/DPiM/cpp/"; export DPCPP;
UTILSCR="${path}/src/utilScripts/"; export UTILSCR;
PERL5LIB="${path}/src/DPiM/scripts/lib:${path}/src/utilScripts/HomeBrew/lib${PERL5LIB+:}${PERL5LIB}"; export PERL5LIB;

