cmake_minimum_required(VERSION 3.10)

set(BUILD_SHARED_LIBS OFF)
find_package(ROOT REQUIRED)

find_package(xgboost REQUIRED)

project(AnalyzeNCPi LANGUAGES CXX)

add_executable(AnalyzeNCPi src/AnalyseNCPi.cxx src/BDTTool.cxx src/CreateTrainingTree.cxx src/EventContainer.cxx src/Selection.cxx src/SelectionDriver.cxx src/StackedHistTool.cxx src/Utility.cxx)
target_link_libraries(AnalyzeNCPi PUBLIC xgboost::xgboost ROOT::Core ROOT::Physics ROOT::Hist ROOT::Tree ROOT::Gpad)
target_include_directories(AnalyzeNCPi PRIVATE include)

ROOT_GENERATE_DICTIONARY(map_analyzenue_dict MODULE AnalyzeNCPi LINKDEF include/LinkDef.h)

project(ProcessNTuple LANGUAGES CXX)

add_executable(ProcessNTuple src/ProcessNTuple.cxx src/BDTTool.cxx src/EventContainer.cxx src/Selection.cxx src/Utility.cxx)
target_link_libraries(ProcessNTuple PUBLIC xgboost::xgboost ROOT::Core ROOT::Physics ROOT::Hist ROOT::Tree ROOT::Gpad)
target_include_directories(ProcessNTuple PRIVATE include)

ROOT_GENERATE_DICTIONARY(map_processntuple_dict MODULE ProcessNTuple LINKDEF include/LinkDef.h)
