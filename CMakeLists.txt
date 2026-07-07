# Specify the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name and specify that it's a C project
# KDevelop defaults to CXX (C++), so changing this is crucial.
project(OFDMsim C)

# Find the FFTW3 library on your system.
# The 'REQUIRED' keyword will make CMake stop with an error if it can't be found.
find_package(FFTW3 REQUIRED)

# Add the directory containing your header files (.h) to the include path.
# "${CMAKE_CURRENT_SOURCE_DIR}" refers to the project's root directory.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/header)
# Tell CMake to create an executable named 'OFDMsim'.
# List all of your .c source files here.
# For now, we only have main.c, but you will add more later.
add_executable(OFDMsim
    main.c
    # Add other .c files here as you create them, like so:
    source/complex_math.c
    source/transmitter.c
    source/channel.c
    source/receiver.c
    source/simulation.c
)

# Link the FFTW3 library to your executable.
# This tells the compiler to connect your code with the functions from FFTW.
target_link_libraries(OFDMsim ${FFTW3_LIBRARIES} m)
