# Define the executable and adding the source files required to build the application.
# When we are defining fortran source files if we name the source files with the extension
# .F90 (using a capital *F* and not a lowercase one) CMake will automatically run the *C*
# preprocessor on the file for us, saving us from specifying that through other means.
# This is only necessary if the file requires the the *C* preprocessor.
add_executable(finite_elasticity_uniaxial finite_elasticity_uniaxial.F90)

# Set the required libraries for this executable.  Here we let CMake know which libraries
# are required to link the application.  Now because the *opencmisslibs* is known to CMake 
# as a target, there is nothing else we need to do.  The *opencmisslibs* target contains
# enough information for CMake to properly build and link the application with the OpenCMISS
# libraries.
#
target_link_libraries(finite_elasticity_uniaxial PUBLIC opencmisslibs)
