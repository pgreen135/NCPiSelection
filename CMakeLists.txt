cmake_minimum_required(VERSION 3.10)

set(BUILD_SHARED_LIBS OFF)
find_package(ROOT REQUIRED)

project(AnalyzeNCPi LANGUAGES CXX)

add_executable(AnalyzeNCPi src/AnalyseNCPi.cxx src/CreateTrainingTree.cxx src/EventContainer.cxx src/Selection.cxx src/SelectionDriver.cxx src/StackedHistTool.cxx src/Utility.cxx src/Hypfit.C src/PhysdEdx.C)
target_link_libraries(AnalyzeNCPi PUBLIC ROOT::Core ROOT::Physics ROOT::Hist ROOT::Tree ROOT::Gpad ROOT::MathCore ROOT::MathMore)
target_include_directories(AnalyzeNCPi PRIVATE include)

ROOT_GENERATE_DICTIONARY(map_analyzenue_dict MODULE AnalyzeNCPi LINKDEF include/LinkDef.h)

project(ProcessNTuple LANGUAGES CXX)

add_executable(ProcessNTuple src/ProcessNTuple.cxx src/EventContainer.cxx src/Selection.cxx src/Utility.cxx src/Hypfit.C src/PhysdEdx.C)
target_link_libraries(ProcessNTuple PUBLIC ROOT::Core ROOT::Physics ROOT::Hist ROOT::Tree ROOT::Gpad ROOT::MathCore ROOT::MathMore)
target_include_directories(ProcessNTuple PRIVATE include)

ROOT_GENERATE_DICTIONARY(map_processntuple_dict MODULE ProcessNTuple LINKDEF include/LinkDef.h)
