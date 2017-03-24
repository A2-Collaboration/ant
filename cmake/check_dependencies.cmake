# check some dependencies in more detail to ensure the build process
# does not fail because of the problems below

# if gsl was missing during ROOT compilation, the mathmore package might be missing
# check if it is available (needed for Interpolator)
execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --features
    OUTPUT_VARIABLE ROOT_FEATURES)
if(NOT ROOT_FEATURES MATCHES ".*mathmore.*")
    message(FATAL_ERROR "ROOT was not build with the mathmore feature.\
        Please make sure gsl or libgsl-dev is installed while compiling ROOT.")
endif()
