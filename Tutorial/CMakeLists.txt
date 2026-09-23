cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(MyThesisAnalysis)

# Force CMake to look inside the professor's custom O2 module paths
list(APPEND CMAKE_MODULE_PATH $ENV{O2_ROOT}/lib/cmake/O2)

# We only strictly need O2Physics for this task
find_package(O2Physics REQUIRED)

add_executable(my-analysis-task myTableCombinations.cxx)
target_link_libraries(my-analysis-task PUBLIC O2Physics::O2Physics)
