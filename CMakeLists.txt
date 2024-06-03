if(UNIX AND NOT APPLE)
    cmake_minimum_required(VERSION 3.23.0)
    project(scarab LANGUAGES C CXX)

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
        message("Using existing gcem clone")
    else()
        execute_process(COMMAND git clone --depth=1 https://github.com/kthohr/gcem)
    endif()

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/LightwaveExplorer)
    message("Using existing Lightwave Explorer clone")
    else()
        execute_process(COMMAND git clone --depth=1 https://github.com/NickKarpowicz/LightwaveExplorer)
    endif()
    
    find_package(PkgConfig REQUIRED)
    find_package(fmt REQUIRED)
    find_package(MKL REQUIRED)
    pkg_check_modules(GTK4 gtk4)

    include_directories(${MKL_ROOT}/include/fftw)
    include_directories(${GTK4_INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)
    include_directories(/usr/local/OceanOptics/OceanDirect/include)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/LightwaveExplorer/Source)
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
        ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.a 
        ${MKL_ROOT}/lib/intel64/libmkl_core.a 
        -Wl,--end-group -lgomp -lpthread -ldl)
    install(TARGETS scarab)
endif()