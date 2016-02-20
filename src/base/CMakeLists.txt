set(ANT_PATH_DATABASE "${CMAKE_SOURCE_DIR}/database")
set(ANT_PATH_GITREPO  "${CMAKE_SOURCE_DIR}")

configure_file(Paths.h.in Paths.h)

set(SRCS_BASE
  types.h
  Detector_t.cc
  ParticleType.cc
  ParticleTypeTree.cc
  interval.h
  piecewise_interval.h
  printable.h
  Logger.cc
  tmpfile_t.cc
  Format.h
  Tree.h
  root_printable.h
  WrapTFile.cc
  mapped_vectors.h
  CmdLine.h
  Paths.h
  matrixstack.cc
  iterators.h
  enumfield.h
  GitInfo.cc
  OptionsList.cc
  ProgressCounter.cc
  TF1Ext.h
  TH1Ext.h
  WrapTTree.h
)

set(SRCS_CEREAL
  cereal/cereal.hpp
  cereal/external/base64.hpp
  cereal/external/rapidjson/document.h
  cereal/external/rapidjson/license.txt
  cereal/external/rapidjson/stringbuffer.h
  cereal/external/rapidjson/internal/stack.h
  cereal/external/rapidjson/internal/pow10.h
  cereal/external/rapidjson/internal/strfunc.h
  cereal/external/rapidjson/rapidjson.h
  cereal/external/rapidjson/prettywriter.h
  cereal/external/rapidjson/writer.h
  cereal/external/rapidjson/genericstream.h
  cereal/external/rapidjson/filestream.h
  cereal/external/rapidjson/reader.h
  cereal/external/rapidxml/license.txt
  cereal/external/rapidxml/rapidxml.hpp
  cereal/external/rapidxml/rapidxml_print.hpp
  cereal/external/rapidxml/rapidxml_iterators.hpp
  cereal/external/rapidxml/rapidxml_utils.hpp
  cereal/external/rapidxml/manual.html
  cereal/macros.hpp
  cereal/details/helpers.hpp
  cereal/details/static_object.hpp
  cereal/details/util.hpp
  cereal/details/polymorphic_impl.hpp
  cereal/details/traits.hpp
  cereal/archives/xml.hpp
  cereal/archives/json.hpp
  cereal/archives/binary.hpp
  cereal/archives/portable_binary.hpp
  cereal/archives/adapters.hpp
  cereal/access.hpp
  cereal/types/map.hpp
  cereal/types/queue.hpp
  cereal/types/string.hpp
  cereal/types/array.hpp
  cereal/types/boost_variant.hpp
  cereal/types/unordered_map.hpp
  cereal/types/forward_list.hpp
  cereal/types/deque.hpp
  cereal/types/common.hpp
  cereal/types/utility.hpp
  cereal/types/memory.hpp
  cereal/types/chrono.hpp
  cereal/types/set.hpp
  cereal/types/tuple.hpp
  cereal/types/base_class.hpp
  cereal/types/vector.hpp
  cereal/types/polymorphic.hpp
  cereal/types/valarray.hpp
  cereal/types/unordered_set.hpp
  cereal/types/complex.hpp
  cereal/types/stack.hpp
  cereal/types/bitset.hpp
  cereal/types/list.hpp
)

set(SRCS_DETAIL
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
  detail/tclap/ValuesConstraintExtra.h
)

set(SRCS_MATH_FUNCTIONS
  math_functions/Faddeeva.cc
  math_functions/Faddeeva.h
  math_functions/voigtian.h
  math_functions/voigtian.cc
)

set(SRCS_STD_EXT
  std_ext/math.h
  std_ext/memory.h
  std_ext/misc.h
  std_ext/string.h
  std_ext/time.h
  std_ext/vector.h
  std_ext/system.cc
  std_ext/convert.h
)

set(SRCS
  ${SRCS_BASE}
  ${SRCS_CEREAL}
  ${SRCS_DETAIL}
  ${SRCS_MATH_FUNCTIONS}
  ${SRCS_STD_EXT}
)

append_header_files(SRCS)

set_source_files_properties(detail/format.cc
        PROPERTIES COMPILE_FLAGS -Wno-effc++)

add_library(base ${SRCS})
target_link_libraries(base ${ROOT_LIBRARIES})
