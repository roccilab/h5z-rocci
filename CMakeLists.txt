find_package(HDF5 REQUIRED)
add_library(
  hdf5rocci
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
export(TARGETS hdf5rocci FILE HDF5ROCCI.cmake)
