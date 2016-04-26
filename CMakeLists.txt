cmake_minimum_required (VERSION 2.8)
project (insitu_error)


FIND_PACKAGE(VTK REQUIRED)
include(${VTK_USE_FILE})
add_definitions(-DVTK)

# settings

set(lib_dir "${CMAKE_SOURCE_DIR}/../lib" CACHE PATH "path to all libs")
set(osuflow_dir "${CMAKE_SOURCE_DIR}/osuflow"   CACHE PATH "path to osuflow library")
set(flowvis_dir "${CMAKE_SOURCE_DIR}/../flowvis" CACHE PATH "path to flowvis project")
option(warnings_on "Show Warnings" OFF)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wno-unused-result")

if (${warnings_on})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif()



# defaults for cmake_build_type
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "DEBUG;RELEASE")

# settings for large data
add_definitions(-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE=1) # large data

# platform settings
IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  add_definitions(-U__STRICT_ANSI__ -D__STDC_FORMAT_MACROS)
  add_definitions(-Doff64_t=off_t -Dlseek64=lseek) # large data
ENDIF()



# automatic settings
SET(libOSUFlow  optimized OSUFlow     debug OSUFlow_d)
SET(libJCLib    optimized jclib     debug jclib)


# directories
INCLUDE_DIRECTORIES( 
    ${osuflow_dir}/include/osuflow       # for make install
    ${osuflow_dir}/include/osuflow/renderer
    ${osuflow_dir}/include/osuflow/renderer/gcb
    ${lib_dir}
    ${flowvis_dir}/common
    )

LINK_DIRECTORIES( 
    ${LINK_DIRECTORIES}
    ${osuflow_dir}/src        
    ${CMAKE_SOURCE_DIR}/osuflow/renderer
    ${CMAKE_SOURCE_DIR}/osuflow/renderer/libgcb
    ${osuflow_dir}/lib/osuflow # for make install
    ${osuflow_dir}/lib/osuflow/renderer
    ${osuflow_dir}/lib/osuflow/renderer/gcb
    ${lib_dir}
    )

# profiling

# openmap
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
message("OPENMP NOT FOUND")
endif()

##################################


add_executable(flow_error flow_error.cpp)
TARGET_LINK_LIBRARIES(flow_error  ${libOSUFlow} ${libJCLib} m)


add_executable(merge_error merge_error.cpp)
TARGET_LINK_LIBRARIES(merge_error ${libOSUFlow} ${libJCLib} m ${VTK_LIBRARIES})

#add_executable(insitu_pathline insitu_pathline.cpp)
#TARGET_LINK_LIBRARIES(insitu_pathline  ${libOSUFlow} ${libJCLib} m)

#add_executable(cmp_pathline cmp_pathline.cpp)
#TARGET_LINK_LIBRARIES(cmp_pathline ${libOSUFlow} ${libJCLib} m ${VTK_LIBRARIES})

