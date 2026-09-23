cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

project(MyThesisAnalysis)

find_package(O2 REQUIRED)
find_package(O2Physics REQUIRED)

add_executable(my-analysis-task myTableCombinations.cxx)

# Explicitly inject the O2 and O2Physics include folders
target_include_directories(my-analysis-task PUBLIC 
    $ENV{O2_ROOT}/include 
    $ENV{O2PHYSICS_ROOT}/include
)

target_link_libraries(my-analysis-task PUBLIC O2::Framework)