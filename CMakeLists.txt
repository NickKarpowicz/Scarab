if(UNIX AND NOT APPLE)
    cmake_minimum_required(VERSION 3.23.0)
    project(scarab LANGUAGES C CXX)
    execute_process(COMMAND git clone https://github.com/kthohr/gcem)
    find_package(PkgConfig REQUIRED)
    find_package(fmt REQUIRED)
    find_package(MKL REQUIRED)
    find_package(TBB REQUIRED)
    pkg_check_modules(GTK4 gtk4)

    include_directories(${MKL_ROOT}/include/fftw)
    include_directories(${GTK4_INCLUDE_DIRS})
    include_directories(build/gcem/include)
    include_directories(/usr/local/OceanOptics/OceanDirect/include)
    include_directories(../LightwaveExplorer/Source)
    link_directories(${GTK4_LIBRARY_DIRS})
    link_directories(/usr/local/OceanOptics/OceanDirect/lib)
    add_compile_options(-std=c++20 -O3 -fopenmp -w -Wall)
    add_executable(scarab 
        Scarab/Scarab.cpp)
    target_link_libraries(scarab ${GTK4_LIBRARIES})
    target_link_libraries(scarab fmt::fmt -lm -loceandirect)
    target_link_libraries(scarab
        -Wl,--start-group 
        ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64.a 
        ${MKL_ROOT}/lib/intel64/libmkl_tbb_thread.a 
        ${MKL_ROOT}/lib/intel64/libmkl_core.a 
        -Wl,--end-group
        ${MKL_ROOT}/../../compiler/latest/lib/libiomp5.a)
    target_link_libraries(scarab TBB::tbb)
    install(TARGETS scarab)
endif()