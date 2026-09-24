cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

project(MyThesisAnalysis)

find_package(O2 REQUIRED)
find_package(O2Physics REQUIRED)

# Point to the official tutorial file
add_executable(tutorial-task collisionTracksIteration.cxx)

target_include_directories(tutorial-task PUBLIC 
    $ENV{O2_ROOT}/include 
    $ENV{O2PHYSICS_ROOT}/include
)

target_link_libraries(tutorial-task PUBLIC O2::Framework)