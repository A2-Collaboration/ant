set(SRCS
  RootDraw.cc
  HistogramFactory.cc
  PromptRandomHist.cc
  CutTree.h
  HistStyle.cc
)

add_library(analysis_plot ${SRCS})
target_link_libraries(analysis_plot third_party_interface)
