# build executable
add_executable(BunchShuffling BunchShuffling.cc)
target_link_libraries(BunchShuffling AsymmetryCalc ${ROOT_LIBRARIES})
