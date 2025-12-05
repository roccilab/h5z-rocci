function(build_hdf5_test)
  get_filename_component(test_name ${ARGV0} NAME_WE)
  add_executable(${test_name} ${ARGV})
  target_link_libraries(${test_name} PUBLIC ROCCI hdf5rocci PRIVATE LibPressio::libpressio LibPressioTools::libpressio_meta)
endfunction(build_hdf5_test)

find_package(LibPressio REQUIRED)
find_package(LibPressioTools REQUIRED)

build_hdf5_test(rocciToHDF5.cc) 
build_hdf5_test(drocciFromHDF5.cc)
build_hdf5_test(convertBinToHDF5.cc)
build_hdf5_test(print_h5repack_args.cc)
build_hdf5_test(tem.cc)