

# Instructions here: https://herbhuang.com/en/posts/old/using-r-on-wsl/

# Open Powershell using search "cmd"
# Enter WSL using
wsl
# check BLAS and LAPACK using
sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu
sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu
# and check R versions using
R
sessionInfo()
q()
# Should show original libraries
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0

# Install openblas using
sudo apt install libopenblas-dev
# update libraries priority
sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu
sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu
# and check R versions using
R
sessionInfo()
q()

