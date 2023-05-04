cmake_minimum_required(VERSION 3.16)

# some CMAKE parameters
set(CMAKE_BUILD_TYPE Debug)
set (CMAKE_CXX_STANDARD 17)

project ( Hzlnt VERSION 0.0.0 )

# some options
option ( Hzlnt_HDF5 "enable HDF5 I/O" OFF )
option ( Hzlnt_MPI "enable MPI parallelism" ON )

# add Kokkos CMake files that will set compiler options
add_subdirectory( external/kokkos )

add_executable( Hzlnt main.cpp )

if(Hzlnt_MPI)
  add_compile_definitions("WITH_MPI")
  find_package(MPI REQUIRED)
  target_link_libraries( Hzlnt MPI::MPI_CXX)
endif()

target_link_libraries( Hzlnt Kokkos::kokkos )
