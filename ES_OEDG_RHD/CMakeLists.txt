cmake_minimum_required(VERSION 3.13.4)
project(rhd_dg)

# Force Release build before deal.II configuration
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

find_package(deal.II 9.7.1 REQUIRED HINTS /home/biswarup-biswas/dealii-install)
deal_ii_initialize_cached_variables()

add_executable(rhd_dg rhd_dg.cc)
deal_ii_setup_target(rhd_dg)

# Add optimization flags for modern CPUs
target_compile_options(rhd_dg PRIVATE -O3 -march=native -ffast-math)

# Custom target to clean EVERYTHING (build artifacts + VTK + PNG + CMake debris)
add_custom_target(cleanall
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND rm -f *.vtk *.png
    COMMAND rm -rf CMakeFiles CMakeCache.txt* cmake_install.cmake Makefile
    COMMAND find . -name "CMakeCache.txt*" -delete
    COMMAND find . -name "*.tmp*" -delete
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Cleaning all build files, simulation outputs, and ALL CMake/temp debris..."
)
