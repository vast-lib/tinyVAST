

# Search "cmd" to open Powershell

# Install Docker from CRAN folks
# Only necessary once, and then shows in Docker Desktop
docker pull rocker/r-devel-ubsan-clang


#
docker run --cap-add SYS_PTRACE --rm -ti --security-opt seccomp=unconfined -v ~c/Users/james/OneDrive/Desktop/Git/tinyVAST:/pkg rocker/r-devel-ubsan-clang
# then use RD which is linked to Rdevel

cd pkg
apt update
apt install cmake
apt install libssl-dev
apt-get install gdal-bin
apt-get install libgdal-dev
apt install libudunits2-dev
apt install libgeos-dev geos-bin
apt install proj-bin
apt install libproj-dev
RD -e "install.packages('remotes')"
RD -e "install.packages('sf')"
RD -e "install.packages('BiocManager'); BiocManager::install('graph')"
RD -e "remotes::install_deps(dependencies = TRUE)"
RD CMD INSTALL --preclean .
