cmake_minimum_required (VERSION 3.4)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project (cuba CXX)
 
find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
 
add_executable (cuba cuba.cpp)
 
target_link_libraries (cuba seqan3::seqan3)
