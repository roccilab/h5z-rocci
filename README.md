# H5Z-ROCCI
ROCCI: Integrated Cyberinfrastructure for In Situ Lossy Compression Optimization Based on Post Hoc Analysis Requirements


# Quickstart with Spack and Libpressio

## 1. Prerequisites (Compiler & Spack)
> **Important Notes:**
> * **Use Pre-built Binaries:** We use pre-built binaries from the Spack cache to speed up the building process. Note that the Spack cache only supports **Linux**.
> * **Compiler Version:** GCC 14/15 treats some warnings as errors, which interrupts the LibPressio build process. Please use **GCC 11, 12, or 13** instead.
> * **Spack Version:** Spack (since v1.0) has known issues with pre-built binary packages. Please use **Spack v0.23.1** strictly.
### Check/Install GCC on Ubuntu
```bash
sudo apt update
sudo apt install -y gcc-13 g++-13 gfortran-13
```
If you are on an HPC cluster, use `module load` instead.
### Check/Install Spack v0.23.1
If you do not have Spack, clone the specific version `v0.23.1`.
```bash
git clone --depth 2 --branch v0.23.1 https://github.com/spack/spack.git
. spack/share/spack/setup-env.sh
```
### Ensure Spack Detects Your System GCC
```bash
spack compiler find
```
## 2. Configure Mirrors & Trust Keys
Add the official binaries and E4S cache, then trust the signing keys to enable binary usage.
```bash
spack mirror add v0.23.1-e4s https://binaries.spack.io/v0.23.1/e4s
spack mirror add v0.23.1 https://binaries.spack.io/v0.23.1
spack buildcache keys --install --trust
```
## 3. Install LibPressio Packages
Run the following command to install LibPressio with standard compressors (SZ, ZFP, SZ3, SZX):
```bash
spack install libpressio-tools+sperr ^libpressio+sz+zfp+sz3+szx ^szx@main
```
> **Tip:** If your Spack compiler list also includes GCC 14 or 15, append `%gcc@13` to the command to strictly enforce the use of GCC 13:
> ```bash
> spack install libpressio-tools+sperr ^libpressio+sz+zfp+sz3+szx ^szx@main %gcc@13
> ```



## 4. Build H5Z-ROCCI:
```bash
git clone https://github.com/roccilab/h5z-rocci.git
cd h5z-rocci
mkdir build && cd build/
cmake -DBUILD_TESTS=ON ..
make
```

Additionally, we must allow hdf5 to find our custom filter by setting `HDF5_PLUGIN_PATH` to the directory containing `libh5z_multi.so`

## 5. Test H5Z-ROCCI

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

The filter parameters can be generated using `conf2int`. 