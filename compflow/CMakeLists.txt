cmake_minimum_required(VERSION 3.3.0)

#-----------------------------------------------------------------------------
# set the project name
project(CompFlow)

#-----------------------------------------------------------------------------
# Specify source directory
set (PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

#-----------------------------------------------------------------------------
# Necessary directories
set(SV_DIR "/usr/local/sv/simvascular/2022-07-20")
set(VTK_DIR "FULL_PATH/vtk/install")
#For include dir for files like "<vtkSmartPointer.h>" (not included in SV thirdparty)
#(local version headers used) 
set(VTK_MAJOR_VERSION "8")  
set(VTK_MINOR_VERSION "2")
set(SV_LIB_PREFIX "_simvascular_") #Libraries start with lib_sim..._
list(APPEND CMAKE_PREFIX_PATH ${SV_DIR}/lib )
list(APPEND CMAKE_PREFIX_PATH ${SV_DIR}/svExternals/lib)

#-----------------------------------------------------------------------------
# Print variables
message("prefix: ${CMAKE_FIND_LIBRARY_PREFIXES}")
message("suffix: ${CMAKE_FIND_LIBRARY_SUFFIXES}")
message("sv dir: ${SV_DIR}")
message("VTK dir: ${VTK_DIR}")
message("")

#Find lib_simvascular_module_simulation.so (for sv4gui_SimulationUtils)
find_library(SIM_LIB "${SV_LIB_PREFIX}module_simulation") 
find_library(VTK_LIB "vtkCommonCore-8.1") #Not sure why, but linking this VTK library was enough

#-----------------------------------------------------------------------------
# Add the executable
add_executable(compflow "${PROJECT_SOURCE_DIR}/compute_flow.cxx")

#-----------------------------------------------------------------------------
#compflow needs C++17 to use std::filesystem
set_property(TARGET compflow PROPERTY CXX_STANDARD 17)

#-----------------------------------------------------------------------------
#Add necessary libraries
target_link_libraries(compflow PUBLIC "${SIM_LIB}" "${VTK_LIB}")

#Add the search path for include files
target_include_directories(compflow PUBLIC "${SV_DIR}/include/core")
target_include_directories(compflow PUBLIC "${VTK_DIR}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")


#-----------------------------------------------------------------------------
# Install by copying executable to src/../bin
install(TARGETS compflow DESTINATION "${PROJECT_SOURCE_DIR}/../bin")
