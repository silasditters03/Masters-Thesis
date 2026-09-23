cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

project(MyThesisAnalysis)

find_package(O2 REQUIRED)
find_package(O2Physics REQUIRED)

add_executable(my-analysis-task myTableCombinations.cxx)

target_link_libraries(my-analysis-task PUBLIC O2::Framework O2Physics::Core)