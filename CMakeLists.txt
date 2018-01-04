project(opengm-converter)
cmake_minimum_required(VERSION 2.8.12)

# C++11
add_compile_options(-std=c++11)

add_subdirectory(external/opengm)
include_directories(external/opengm/include)

# for opengm
find_package(HDF5 1.8.15 REQUIRED)
add_library(HDF5 INTERFACE)
target_include_directories(HDF5 INTERFACE ${HDF5_INCLUDE_DIR})
target_link_libraries(HDF5 INTERFACE ${HDF5_LIBRARIES})
target_compile_definitions(HDF5 INTERFACE ${HDF5_DEFINITIONS})
message(STATUS ${HDF5_LIBRARIES})
message(STATUS ${HDF5_INCLUDE_DIR})

add_executable(convert_mrf_to_uai convert_mrf_to_uai.cpp)
target_link_libraries( convert_mrf_to_uai HDF5 ) 

add_executable(convert_multicut_to_text convert_multicut_to_text.cpp)
target_link_libraries( convert_multicut_to_text HDF5 ) 
