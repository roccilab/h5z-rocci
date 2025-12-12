find_package(LibPressio REQUIRED)
find_package(LibPressioTools REQUIRED)

find_package(HDF5 REQUIRED COMPONENTS C)

function(build_hdf5_test)
  get_filename_component(test_name ${ARGV0} NAME_WE)
  add_executable(${test_name} ${ARGV})
  target_link_libraries(${test_name} PUBLIC ROCCI hdf5rocci PRIVATE LibPressioSperr::libpressio_sperr LibPressio::libpressio LibPressioTools::libpressio_meta)
  # target_link_libraries(${test_name} PUBLIC ROCCI hdf5rocci PRIVATE LibPressio::libpressio LibPressioTools::libpressio_meta)
endfunction(build_hdf5_test)

build_hdf5_test(select.cpp)
# build_hdf5_test(tem.cpp)
build_hdf5_test(conf2int.cpp)

# add_executable( tem tem.cpp )
# target_link_libraries( tem PUBLIC ROCCI LibPressioTools::libpressio_meta )

set(CMAKE_CXX_STANDARD 17)

add_library(h5z_multi SHARED h5z_multi.cpp)

target_include_directories(h5z_multi
  PUBLIC  ${HDF5_INCLUDE_DIRS}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include/surrogates
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include/QoZ/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src/surrogates
)

target_link_options(h5z_multi PRIVATE -static-libstdc++ -static-libgcc)

target_link_libraries(h5z_multi
  PUBLIC  ${HDF5_LIBRARIES}
  PRIVATE SLYCE LibPressio::libpressio LibPressioTools::libpressio_meta
)

set_target_properties(h5z_multi PROPERTIES
  C_VISIBILITY_PRESET   default
  CXX_VISIBILITY_PRESET default
  VISIBILITY_INLINES_HIDDEN OFF
  OUTPUT_NAME "h5z_multi"
  PREFIX "lib"
  SUFFIX ".so"
)

# target_compile_options(SLYCE PRIVATE -H)