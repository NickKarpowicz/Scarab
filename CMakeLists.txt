if(UNIX AND NOT APPLE)
    cmake_minimum_required(VERSION 3.23.0)
    project(scarab LANGUAGES C CXX)

    if(NOT FLATPAK)
        if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gcem)
            message("Using existing gcem clone")
        else()
            execute_process(COMMAND git clone --depth=1 https://github.com/kthohr/gcem)
        endif()
    endif()

    find_package(PkgConfig REQUIRED)
    find_package(fmt REQUIRED)
    if(FLATPAK)
        find_package(MKL PATHS /app/opt/intel/oneapi/mkl/latest/lib/cmake)
        include_directories(/app/include)
        link_directories(/app/lib)
    else()
        find_package(MKL REQUIRED)
        include_directories(/usr/local/OceanOptics/OceanDirect/include)
        link_directories(/usr/local/OceanOptics/OceanDirect/lib)
    endif()

    pkg_check_modules(GTK4 gtk4)

    include_directories(${MKL_ROOT}/include/fftw)
    include_directories(${MKL_ROOT}/include)
    include_directories(${GTK4_INCLUDE_DIRS})
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/gcem/include)

    link_directories(${GTK4_LIBRARY_DIRS})

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
