# Getting started

## Build/install

Make sure you have NGSolve installed. Building the package should be straight forward:
```
git clone https://github.com/NGSolve/ngs-special-functions.git
cd ngs-special-functions
mkdir build && cd build
cmake ..
# automatically installs during build to your NGSolve installation
make -j4 
cd ..
```

## Run examples
```
cd demo
netgen gammaln.py
```
