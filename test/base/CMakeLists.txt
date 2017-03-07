add_ant_test(StdExt)
add_ant_test(Intervals)
add_ant_test(Tree)
add_ant_test(WrapTFile)
add_ant_test(Detector_t)
add_ant_test(OptionsList)
add_ant_test(ParticleType)
add_ant_test(Vec)
add_ant_test(Interpolator)
add_ant_test(Printable)
add_ant_test(FloodFillAverages)
add_ant_test(SavitzkyGolay)

add_ant_test(WrapTTree)
# fixes strange bug with nasty test in Release mode
set_source_files_properties(TestWrapTTree.cc
                            PROPERTIES COMPILE_FLAGS "-O0")
