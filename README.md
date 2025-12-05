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

Additionally, we must allow hdf5 to find our custom filter by setting `HDF5_PLUGIN_PATH` to the directory containing `libhdf5rocci.so`:

```bash
export HDF5_PLUGIN_PATH=/home/projects/h5z-rocci/build
```

The filter id is 32088, with some parameters to control compression mode and compression parameters. 

Before using the filter, please use the following command to update the plugin path and library path. 

```bash
export HDF5_PLUGIN_PATH=${HDF5_PLUGIN_PATH}:/your/path/to/filter
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/your/path/to/this/project/h5z-rocci/build
```

For the compression, for example, use the following command to do compression and decompression. 

```bash
# compress with fixed error bound 1e-4 using best compressor
h5repack -f UD=32088,0,3,3944497965,1058682594,256 -i original_file.h5 -o compressed_file.h5
# decompress
h5repack -f NONE -i compressed_file.h5 -o decompressed_file.h5
```

Additionally, you can test the filter using the examples `rocciToHDF5` and `drocciFromHDF5` in `build/test`. To alter the configuration you can alter `rocci.config` in `test/`, this file should be in the working directory when calling the H5Z filter to ensure it is found.

```bash
# build/test/
./rocciToHDF5 -f testfloat_8_8_128.dat 8 8 128
./drocciFromHDF5 testfloat_8_8_128.dat.rocci.h5
```

*Note: while the C API usage in the examples above has been confirmed to work with hdf5 1.14.2, using utilities such as h5repack and h5dump may require a different version of hdf5. Currently, we have confirmed that hdf5 1.10.1 works for both the C API and the hdf5 utility scripts*