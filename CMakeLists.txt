cmake_minimum_required(VERSION 3.16)

# some CMAKE parameters
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project ( Hzlnt VERSION 0.0.0 )

# some options
option ( Hzlnt_HDF5 "enable HDF5 I/O" OFF )
option ( Hzlnt_HERCULE "enable HERCULE I/O" ON )
option ( Hzlnt_MPI "enable MPI parallelism" ON )

# Hercule config
# Hercule package
if(${Hzlnt_HERCULE})
    set(Hercule_DIR "$ENV{HERCULE_ROOT_DIR}/lib/Hercule/cmake")
    message(STATUS "[Hercule_DIR] : ${Hercule_DIR}")
    set( MachineTypes_DIR "$ENV{HERCULE_ROOT_DIR}/lib/MachineTypes/cmake" )
    set( MachineIo64_DIR "$ENV{HERCULE_ROOT_DIR}/lib/MachineIo64/cmake" )
    set( Lm_DIR "$ENV{HERCULE_ROOT_DIR}/lib/Lm/cmake" )
    set( Gme_DIR "$ENV{HERCULE_ROOT_DIR}/lib/Gme/cmake" )
endif()

# add Kokkos CMake files that will set compiler options
add_subdirectory( external/kokkos )

add_subdirectory( src/io )

add_executable( Hzlnt main.cpp )
target_link_libraries( Hzlnt Hzlnt_io )

if(Hzlnt_MPI)
  add_compile_definitions("WITH_MPI")
  find_package(MPI REQUIRED)
  target_link_libraries( Hzlnt MPI::MPI_CXX)
endif()

target_link_libraries( Hzlnt Kokkos::kokkos )
