include_directories(.)

set(ANT_PATH_DATABASE "${CMAKE_SOURCE_DIR}/database")
configure_file(Paths.h.in Paths.h)

set(SRCS
  types.cc
  Detector_t.cc
  ParticleType.cc
  interval.h
  piecewise_interval.h
  printable.cc
  Logger.cc
  tmpfile_t.cc
  std_ext.h
  Format.h
  Tree.h
  root_printable.h
  WrapTFile.cc
  ReadTFiles.cc
  mapped_vectors.h
  CmdLine.h
  Paths.h
  detail/easylogging++.h
  detail/format.cc
  detail/CmdLine.h
  detail/tclap/ArgException.h
  detail/tclap/Arg.h
  detail/tclap/ArgTraits.h
  detail/tclap/CmdLineInterface.h
  detail/tclap/CmdLineOutput.h
  detail/tclap/Constraint.h
  detail/tclap/DocBookOutput.h
  detail/tclap/HelpVisitor.h
  detail/tclap/IgnoreRestVisitor.h
  detail/tclap/MultiArg.h
  detail/tclap/MultiSwitchArg.h
  detail/tclap/OptionalUnlabeledTracker.h
  detail/tclap/StandardTraits.h
  detail/tclap/StdOutput.h
  detail/tclap/SwitchArg.h
  detail/tclap/UnlabeledMultiArg.h
  detail/tclap/UnlabeledValueArg.h
  detail/tclap/ValueArg.h
  detail/tclap/ValuesConstraint.h
  detail/tclap/VersionVisitor.h
  detail/tclap/Visitor.h
  detail/tclap/XorHandler.h
  detail/tclap/ZshCompletionOutput.h
)

set_source_files_properties(detail/format.cc
        PROPERTIES COMPILE_FLAGS -Wno-effc++)

add_library(base ${SRCS})
target_link_libraries(base ${ROOT_LIBRARIES})
