set(ANT_BASE
    base/types.cc
    base/types.h
    ParticleType.h
    ParticleType.cc
    Track.h
    Track.cc
    Particle.h
    Particle.cc
    Event.h
    Event.cc
    EventManager.h
    EventManager.cc
    AntPhysics.h
    AntPhysics.cc
    base/printable.h
    base/printable.cc
    utils/combinatorics.h
    plot/Histogram.h
    plot/Histogram.cc
    base/interval.h
    TaggerHit.h
    TaggerHit.cc
    plot/root_draw.h
    plot/root_draw.cc
    Detector.h
    Detector.cc
    utils/matcher.h
    TriggerInfo.h
    TriggerInfo.cc
    A2GeoAcceptance.h
    A2GeoAcceptance.cc
    plot/SmartHist.h
    plot/SmartHist.cc
    plot/HistogramFactories.h
    plot/HistogramFactories.cc
    )

set(ANT_PHYSICS
    analysis/TestPhysics.h
    analysis/TestPhysics.cc
    analysis/DeltaPlusPhysics.h
    analysis/DeltaPlusPhysics.cc
    analysis/omega.h
    analysis/omega.cc
    analysis/omega_bottomup.h
    analysis/omega_bottomup.cc
    analysis/GeoAcceptance.h
    analysis/GeoAcceptance.cc
    analysis/TestAPLCON.h
    analysis/TestAPLCON.cc
    analysis/DataOverview.h
    analysis/DataOverview.cc
    analysis/mctrue_acceptance.h
    analysis/mctrue_acceptance.cc
    )



set(GOAT_EXTRA
    inc/GDataChecks.h
    src/GDataChecks.cc
    inc/GParticleReconstruction.h
    src/GParticleReconstruction.cc
    inc/GMesonReconstruction.h
    src/GMesonReconstruction.cc
    inc/GSort.h
    src/GSort.cc
    inc/GoAT.h
    src/GoAT.cc
    )

# look for APLCON fitter package
# don't forget to build the library  
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH})
find_package(APLCON REQUIRED)
  
foreach(f ${GOAT_BASE} ${GOAT_EXTRA})
    list(APPEND ANT_GOAT_BASE ../${f})
endforeach()

add_library(antlib SHARED ${ANT_BASE} ${ANT_PHYSICS})

include_directories(.)
include_directories(${APLCON_INCLUDE_DIR})

aux_source_directory("progs" ANT_PROGS)
foreach(PROG_SRC ${ANT_PROGS})

    get_filename_component(name ${PROG_SRC} NAME_WE)
    message(STATUS "Adding ant executable: ${name}")
    add_executable(${name}
        ${PROG_SRC}
        )

    target_link_libraries(${name}
    	antlib
    	${APLCON_LIBRARIES}
        goatbase
        ${LIBS}
    	${ROOT_LIBRARIES}
    )
endforeach()
