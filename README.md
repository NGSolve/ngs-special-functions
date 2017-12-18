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

## Extend

You can extend by adding further functions from slatec. Follow these steps: 

1. In `src/fetch_and_convert_slatec.cmake` add an additional line: `fetch_and_convert_slatec_sources(<what_you_want_to_add>)`
2. Modify `src/specialcf.cpp` to recognize the additional function.
3. Rebuild package (using `cmake` and `make` as above).


