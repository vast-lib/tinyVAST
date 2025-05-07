
#############
# Option-1 ... working
#############

# search "cmd" to open Powershell

# Update stuff by typing in powershell
wsl
sudo apt update
sudo apt install r-base build-essential clang

# Use development version if needed
sudo apt install r-base-dev

# Set compiler flags for Clang
mkdir -p ~/.R
nano ~/.R/Makevars
# and copy-paste the following into nano
CFLAGS += -fsanitize=undefined -fno-omit-frame-pointer -fno-sanitize-recover=all
CXXFLAGS += -fsanitize=undefined -fno-omit-frame-pointer -fno-sanitize-recover=all
# and then hit "ctrl+O" and then "Enter" to save, and "ctrl+X" to exit Nano

# Make a tar.gz of the current build

# navigate to that directory
cd [tar.gz directory]

# Run R CMD check
R CMD check [tar.gz name, e.g., tinyVAST_1.0.0.tar.gz]




#############
# Option-2 not working
#############

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
