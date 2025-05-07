
# Get dependencies
sudo apt update
sudo apt install -y \
  build-essential \
  gfortran \
  libreadline-dev \
  libx11-dev \
  libxt-dev \
  libpng-dev \
  libjpeg-dev \
  libcairo2-dev \
  libbz2-dev \
  libzstd-dev \
  liblzma-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  texinfo \
  texlive \
  texlive-fonts-extra \
  texlive-latex-extra \
  wget \
  git \
  subversion \
  gdebi-core

# Get R-devel
svn checkout https://svn.r-project.org/R/trunk R-devel
cd R-devel

