# check some dependencies in more detail to ensure the build process
# does not fail because of the problems below

find_package(GSL REQUIRED)
# gsl is required for ant to build
if(NOT GSL_FOUND)
    message(FATAL_ERROR "gsl is needed! Please make sure to install libgsl-dev or gsl, \
        depending on your Linux distribution.")
endif()

# check that the GSL include path is not empty --> gsl headers missing, required for ant/base
if("${GSL_INCLUDE_DIR}" STREQUAL "")
    message(FATAL_ERROR "gsl include directory empty, please make sure the \
        gsl-dev package of your distribution is installed.")
endif()

# if gsl was missing during ROOT compilation, the mathmore package might be missing
# check if it is available (needed for Interpolator)
execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --features
    OUTPUT_VARIABLE ROOT_FEATURES)
if(NOT ROOT_FEATURES MATCHES ".*mathmore.*")
    message(FATAL_ERROR "ROOT was not build with the mathmore feature.\
        Please make sure gsl or libgsl-dev is installed while compiling ROOT.")
endif()
