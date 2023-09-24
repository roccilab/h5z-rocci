find_package(HDF5 REQUIRED)
find_package(std_compat REQUIRED)
find_package(LibPressio REQUIRED)

include(CMakePrintHelpers)
message(HDF5 STUFF:)
cmake_print_variables(HDF5_LIBRARIES)
cmake_print_variables(HDF5_INCLUDE_DIRS)

add_library(ROCCI SHARED
  src/rocci.c
  src/rocci_ByteToolkit.c
  src/rocci_FileUtil.c
)
target_include_directories(
  ROCCI
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:ROCCI>
)

target_link_libraries(
  ROCCI
  PUBLIC LibPressio::libpressio
)

add_library(
  hdf5rocci SHARED
  src/H5Z_ROCCI.c
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
export(TARGETS hdf5rocci ROCCI FILE HDF5ROCCI.cmake)
