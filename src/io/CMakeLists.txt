
# list of source file for directory 'src/io'
set(IO_SRC_FILES "")
set(IO_LIB_LINK "")

if(${Hzlnt_HERCULE})

    if( NOT EXISTS ${Hercule_DIR} )
        message(FATAL_ERROR "Hercule_DIR does not exists : " ${Hercule_DIR})
    endif()

    find_package(Hercule REQUIRED)
    if( ${Hercule_FOUND} )
        message(STATUS "Found Hercule library")
    endif()

    message(STATUS "Hercule headers: ${Hercule_DIR}/include")

    file(GLOB IO_SRC_HERCULE CONFIGURE_DEPENDS "./Reader_Hercule*")
    
    # list of files regarding Hercule
    list(APPEND IO_SRC_FILES ${IO_SRC_HERCULE})

    # Hercule link
    list(APPEND IO_LIB_LINK Hercule::hercule_hic)
endif()

message(STATUS "[src/io] Sources files: ${IO_SRC_FILES}")

add_library(Hzlnt_io SHARED ${IO_SRC_FILES})
set_target_properties(Hzlnt_io PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(Hzlnt_io PUBLIC ${IO_LIB_LINK} Kokkos::kokkos)