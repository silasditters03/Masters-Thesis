cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

project(MyThesisAnalysis)

find_package(O2 REQUIRED)
find_package(O2Physics REQUIRED)

# Point to the official tutorial files: <executable name> <source file>
add_executable(tutorial-task collisionTracksIteration.cxx)
add_executable(tutorial-histograms histograms.cxx)
add_executable(tutorial-newcollections newCollections.cxx)
add_executable(tutorial-associated-example associatedExample.cxx)

set(TUTORIAL_TARGETS
    tutorial-task
    tutorial-histograms
    tutorial-newcollections
    tutorial-associated-example
)

foreach(target IN LISTS TUTORIAL_TARGETS)
  # Include directories
  target_include_directories(${target} PUBLIC
      $ENV{O2_ROOT}/include
      $ENV{O2PHYSICS_ROOT}/include
  )

  # Link libraries
  target_link_libraries(${target} PUBLIC O2::Framework)
endforeach()