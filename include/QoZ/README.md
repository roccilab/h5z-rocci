For users: This is the developing branch which is unstable and may have bugs. Do not use it for test.

QoZ: Dynamic Quality Metric Oriented Error Bounded Lossy Compression for Scientific Datasets

## Introduction

This is the source code of QoZ: Dynamic Quality Metric Oriented Error Bounded Lossy Compression for Scientific Datasets

## Dependencies

Please Installing the following dependencies before running the artiact evaluation experiments:

* Python >= 3.6
* numpy 
* pandas 
* qcat (from https://github.com/Meso272/qcat, check its readme for installation guides. Make sure the following executables are successfully installed: calculateSSIM and computeErrAutoCorrelation)

## 3rd party libraries/tools

* Zstandard (https://facebook.github.io/zstd/). Not mandatory to be mannually installed as Zstandard v1.4.5 is included and will be used if libzstd can not be found by
  pkg-config.

## Installation

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include




