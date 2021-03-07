# Bringing more special functions to NGSolve 

This is an add-on package to NGSolve, which shows
how to bring special functions implemented elsewhere into 
NGSolve's python interface. We bring in the Bessel and Hankel
function implementations of the [Slatec package](http://www.netlib.org/slatec/) 
and the error function implementation from the 
[Faddeeva package](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package).


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

## Extend

To extend by adding further functions from Slatec,  follow these steps: 

1. In `src/fetch_and_convert_slatec.cmake` add an additional line: `fetch_and_convert_slatec_sources(<what_you_want_to_add>)`
2. Modify `src/specialcf.cpp` to recognize the additional function.
3. Rebuild package (using `cmake` and `make` as above).


