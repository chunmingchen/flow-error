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


# directories
INCLUDE_DIRECTORIES( 
    ${osuflow_dir}/include       # for make install
    ${osuflow_dir}/include/renderer
    ${osuflow_dir}/include/renderer/gcb
    ${lib_dir}/jclib
    ${flowvis_dir}/common
    )

LINK_DIRECTORIES( 
    ${LINK_DIRECTORIES}
    ${osuflow_dir}/src        
    ${CMAKE_SOURCE_DIR}/osuflow/renderer
    ${CMAKE_SOURCE_DIR}/osuflow/renderer/libgcb
    ${osuflow_dir}/lib # for make install
    ${osuflow_dir}/lib/renderer
    ${osuflow_dir}/lib/renderer/gcb
    ${lib_dir}
    )

SET(libOSUFlow  optimized OSUFlow_r     debug OSUFlow_d)
SET(libJCLib    optimized jclib     debug jclibd)

# profiling

##################################


add_executable(flow_error flow_error.cpp)
TARGET_LINK_LIBRARIES(flow_error  ${libOSUFlow} ${libJCLib} m)

add_executable(merge_error merge_error.cpp)
TARGET_LINK_LIBRARIES(merge_error ${libOSUFlow} ${libJCLib} m ${VTK_LIBRARIES})

#add_executable(cmp_pathline cmp_pathline.cpp)
#TARGET_LINK_LIBRARIES(cmp_pathline ${libOSUFlow} ${libJCLib} m ${VTK_LIBRARIES})

