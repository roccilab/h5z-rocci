project(H5Z_ROCCI CXX C)

find_package(HDF5 REQUIRED)
find_package(std_compat REQUIRED)
find_package(LibPressio REQUIRED)
find_package(LibPressioTools REQUIRED)
find_package(LibPressioOpt REQUIRED)
find_package(LibDistributed REQUIRED)
find_library(UUID_LIBRARY NAMES uuid libuuid REQUIRED)

option(USE_SZX "Build SLYCE SZx Surrogate" ON)
option(USE_SZ3 "Build SLYCE SZ3 Surrogate" ON)
option(USE_ZFP "Build SLYCE ZFP Surrogate" ON)
option(USE_SZP "Build SLYCE SZP Surrogate" ON)
# option(USE_CUSZP "Build SLYCE CUSZP Surrogate" ON)
option(USE_SPERR "Build SLYCE ZFP Surrogate" ON)

if(USE_SZX)
  list(APPEND slyce_src 
    src/surrogates/SZx/estimate_metric.cc
    src/surrogates/SZx/sampling.cc
    src/surrogates/SZx/szx_float.cc
    src/surrogates/SZx/SZxSurrogate.cc
  )
endif()

if(USE_SZ3)
  list(APPEND slyce_src 
    src/surrogates/SZ3/estimate_cr.cc
    src/surrogates/SZ3/estimate_psnr.cc
    src/surrogates/SZ3/estimate_ssim.cc
    src/surrogates/SZ3/CustomHuffmanEncoder.cc
    src/surrogates/SZ3/SZ3Surrogate.cc
  )
endif()

if(USE_ZFP)
  list(APPEND slyce_src 
    src/surrogates/ZFP/estimate_metric.cc
    src/surrogates/ZFP/ZFPSurrogate.cc
  )
endif()

if(USE_SZP)
  find_package(OpenMP REQUIRED)
  list(APPEND slyce_src 
    src/surrogates/SZp/estimate_metric.cc
    src/surrogates/SZp/SZp_float.cc
    src/surrogates/SZp/SZpSurrogate.cc
  )
endif()

if(USE_CUSZP)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  list(APPEND slyce_src 
    src/surrogates/cuSZp/cuSZp_entry.cu
    src/surrogates/cuSZp/estimatef32.cu
    src/surrogates/cuSZp/cuSZpSurrogate.cc
  )
endif()

if(USE_SPERR)
  list(APPEND slyce_src 
    src/surrogates/SPERR/estimate_metric.cc
    src/surrogates/SPERR/SPERRSurrogate.cc
  )
endif()

add_library(SLYCE SHARED
  src/surrogates/qcat_ssim.cc
  ${slyce_src}
)
target_include_directories(
  SLYCE
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/surrogates>
  $<INSTALL_INTERFACE:SLYCE>
  ${CMAKE_CURRENT_SOURCE_DIR}/include/QoZ/include
)
target_link_libraries(
  SLYCE
  PUBLIC LibPressio::libpressio
)

if(USE_SZP)
  target_link_libraries(
    SLYCE
    PUBLIC OpenMP::OpenMP_CXX
  )
endif()

if(USE_CUSZP)
  set_property(TARGET SLYCE PROPERTY CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(
    SLYCE
    PRIVATE CUDA::cuda_driver
  )
endif()

add_library(ROCCI SHARED
  src/rocci.cc
  src/rocci_ByteToolkit.cc
  src/rocci_FileUtil.cc
  src/rocci_utils.cc
  src/inih/ini.c
)
set_target_properties(ROCCI PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(
  ROCCI
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:ROCCI>
)

target_link_libraries(
  ROCCI
  PUBLIC LibPressio::libpressio LibPressioOpt::libpressio_opt LibDistributed::libdistributed ${UUID_LIBRARY} SLYCE
)

add_library(
  hdf5rocci SHARED
  src/H5Z_ROCCI.cc
  )
target_link_libraries(
  hdf5rocci
  PUBLIC ROCCI ${HDF5_LIBRARIES}
  )
target_include_directories(
  hdf5rocci
  PUBLIC
  ${HDF5_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:hdf5_rocci>
  )

if(BUILD_TESTS)
  add_subdirectory(test)
endif()

add_subdirectory(tools)

install(TARGETS hdf5rocci EXPORT HDF5ROCCI
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_rocci)
export(TARGETS hdf5rocci ROCCI SLYCE FILE HDF5ROCCI.cmake)
