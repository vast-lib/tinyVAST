
# Start -> Search "cmd"
# To install R in command-line:
#  1.  "wsl --install"
#  2.  "wsl sudo apt update"
#  3.  "wsl sudo apt install r-base"
#  4.  "wsl sudo apt install cmake"        # CPP compiles
#  5.  "wsl sudo apt install libssl-dev"      # Unknown
#  5B.  "wsl sudo apt install libudunits2-dev"         # units dependency
#  5C.  "wsl sudo apt install proj-bin libproj-dev"    # PROJ
#  6.  "wsl sudo apt-get install gdal-bin"     # GDAL for spatial stuff
#  7.  "wsl sudo apt-get install libgdal-dev"  # ditto
#  8.  "wsl sudo apt install valgrind"   # valgrind
#  9.  "wsl git clone https://github.com/vast-lib/tinyVAST"
# To open R in command-line:  "wsl R"
# Then use R as normal

# To RUN valgrind
# Start -> Search "cmd"
# cd [copy windows explorer, e.g., C:\Users\james\OneDrive\Desktop\Git\tinyVAST]
# wsl R -d "valgrind" -f scratch/test.R      #  remove  --leak-check=full to match CRAN

#