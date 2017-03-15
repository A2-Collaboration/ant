set(SRCS
  root_draw.cc
  HistogramFactory.cc
  PromptRandomHist.cc
  CutTree.h
  HistStyle.cc
)

add_library(analysis_plot ${SRCS})