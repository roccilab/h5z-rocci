project(H5Z_ROCCI CXX C)

find_package(HDF5 REQUIRED)
find_package(std_compat REQUIRED)
find_package(LibPressio REQUIRED)
find_package(LibPressioOpt REQUIRED)
find_package(LibDistributed REQUIRED)
find_library(UUID_LIBRARY NAMES uuid libuuid REQUIRED)

# include(CMakePrintHelpers)
# message(HDF5 DIRS:)
# cmake_print_variables(HDF5_LIBRARIES)
# cmake_print_variables(HDF5_INCLUDE_DIRS)
add_library(SECRE SHARED
  src/surrogates/qcat_ssim.cc
  src/surrogates/SZx/estimate_metric.cc
  src/surrogates/SZx/sampling.cc
  src/surrogates/SZx/szx_float.cc
  src/surrogates/SZx/SZxSurrogate.cc

  src/surrogates/SZ3/estimate_cr.cc
  src/surrogates/SZ3/estimate_psnr.cc
  src/surrogates/SZ3/estimate_ssim.cc
  src/surrogates/SZ3/CustomHuffmanEncoder.cc
  src/surrogates/SZ3/SZ3Surrogate.cc

  src/surrogates/ZFP/estimate_metric.cc
  src/surrogates/ZFP/ZFPSurrogate.cc

)
target_include_directories(
  SECRE
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/surrogates>
  $<INSTALL_INTERFACE:SECRE>
)
target_link_libraries(
  SECRE
  PUBLIC LibPressio::libpressio
)

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
  PUBLIC LibPressio::libpressio LibPressioOpt::libpressio_opt LibDistributed::libdistributed ${UUID_LIBRARY} SECRE
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

install(TARGETS hdf5rocci EXPORT HDF5ROCCI
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hdf5_rocci)
export(TARGETS hdf5rocci ROCCI SECRE FILE HDF5ROCCI.cmake)
