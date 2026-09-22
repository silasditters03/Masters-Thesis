cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(MyThesisAnalysis)

# Find the O2 and O2Physics libraries loaded in the environment
find_package(O2 REQUIRED)
find_package(O2Physics REQUIRED)

# Create the executable from your C++ file
add_executable(my-analysis-task myTableCombinations.cxx)

# Link your executable against the ALICE O2 framework
target_link_libraries(my-analysis-task PUBLIC 
    O2Physics::O2Physics 
    O2::O2
)
