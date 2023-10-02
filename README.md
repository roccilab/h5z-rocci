# H5Z-ROCCI
ROCCI: Integrated Cyberinfrastructure for In Situ Lossy Compression Optimization Based on Post Hoc Analysis Requirements


# Quickstart with Spack
First we'll set up a spack environment:
```bash
git clone https://github.com/spack/spack
git clone https://github.com/robertu94/spack_packages robertu94_packages
source ./spack/share/spack/setup-env.sh
spack repo add robertu94_packages

spack env create -d rocci_env
spack env activate rocci_env/
spack -e rocci_env config add modules:prefix_inspections:lib64:[LD_LIBRARY_PATH]
spack -e rocci_env config add modules:prefix_inspections:lib:[LD_LIBRARY_PATH]
spack add libpressio-tools ^ libpressio+sz+zfp+sz3 # add compressor plugins as needed here

spack concretize -f
spack install
```

Then we can build H5Z-ROCCI:
```bash
git clone https://github.com/roccilab/h5z-rocci.git
cd h5z-rocci
mkdir build && cd build/
cmake -DBUILD_TESTS=ON ..
make
```

Now you can test the filter using the examples `rocciToHDF5` and `drocciFromHDF5` in `build/test`. To alter the configuration you can alter `rocci.config` in `test/`, this file should be in the working directory when calling the H5Z filter to ensure it is found.

```bash
# build/test/
./rocciToHDF5 -f testfloat_8_8_128.dat 8 8 128
./drocciFromHDF5 testfloat_8_8_128.dat.rocci.h5
```

*Note: while the C API usage in the examples above has been confirmed to work with hdf5 1.14.2, using utilities such as h5repack and h5dump may require a different version of hdf5. Currently, we have confirmed that hdf5 1.10.1 works for both the C API and the hdf5 utility scripts*