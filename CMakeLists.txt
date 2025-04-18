cmake_minimum_required (VERSION 3.16)

project (adaptive_map_optimization
        VERSION 0.0.1
        LANGUAGES CXX
)

set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
#set(EXTERNAL_DIR "/Users/FloriGod/Desktop/Masterarbeit/AdaptiveMapOptimization/external/")
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")

# find_package(OpenMesh REQUIRED)

if(NOT TARGET OpenMeshCore)
    FetchContent_Declare(openmesh
            GIT_REPOSITORY  https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
            GIT_TAG OpenMesh-11.0
            SOURCE_DIR "${EXTERNAL_DIR}/OpenMesh"
            )
    FetchContent_MakeAvailable(openmesh)
    include_directories("${EXTERNAL_DIR}/OpenMesh/src/")
endif()


# find_package(Eigen3 QUIET)
# if(NOT TARGET Eigen3::Eigen)
#     FetchContent_Declare(eigen
#             #GIT_REPOSITORY https://gitlab.com/libeigen/eigen
#             #GIT_TAG 3.4.0
#             GIT_REPOSITORY https://gitlab.com/mheistermann/eigen
#             GIT_TAG fix/spqr-index-vs-storageindex
#             SOURCE_DIR "${EXTERNAL_DIR}/eigen"
#     )
#     FetchContent_Populate(eigen)
#     message("Downloaded Eigen3 to ${eigen_SOURCE_DIR}")
#     add_library(Eigen3::Eigen INTERFACE IMPORTED)
#     target_include_directories(Eigen3::Eigen INTERFACE "$<BUILD_INTERFACE:${eigen_SOURCE_DIR}>")

#     #target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
# endif()

find_package(Eigen3 QUIET)
if (NOT EIGEN3_FOUND)
    message("Eigen3 not found, downloading...")
    include(FetchContent)
    FetchContent_Declare(
            eigen3
            URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
    )
    FetchContent_Populate(eigen3)

    add_library(eigen INTERFACE)
    add_library(Eigen3::Eigen ALIAS eigen)
    target_include_directories(eigen INTERFACE
            $<BUILD_INTERFACE:${eigen3_SOURCE_DIR}>
            )
endif ()


if(NOT TARGET TinyAD)
    FetchContent_Declare(tinyad
            GIT_REPOSITORY  https://github.com/patr-schm/TinyAD.git
            GIT_TAG main
            SOURCE_DIR "${EXTERNAL_DIR}/TinyAD"
    )
    FetchContent_MakeAvailable(tinyad)
    set(TINYAD_INCLUDE_DIR "${EXTERNAL_DIR}/TinyAD/include")
endif()


if(NOT TARGET CoMISo)
    FetchContent_Declare(comiso
            GIT_REPOSITORY https://graphics.rwth-aachen.de:9000/CoMISo/CoMISo.git
            GIT_TAG master
            SOURCE_DIR "${EXTERNAL_DIR}/CoMISo"
    )
    option(COMISO_NO_INSTALL "avoid cmake install/export (hack)" ON)
    FetchContent_MakeAvailable(comiso)
endif()


add_executable(adaptive_map_optimization
        src/adaptive_map.cc
        src/adaptive_map.hh

        src/helpers.cc
        src/helpers.hh

        src/includes.hh

        src/OptimizationElements2D.hh

        src/RemesherBase.hh
        src/RemesherValentin.hh
        src/RemesherValentin.cc
        src/RemesherNaive.hh
        src/RemesherNaive.cc
        src/RemeshingAssisted.hh
        src/RemeshingAssisted.cc





    )

target_include_directories(adaptive_map_optimization
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<BUILD_INTERFACE:${eigen3_SOURCE_DIR}>
        ${EXTERNAL_DIR}/OpenMesh/src/
)
target_link_libraries(adaptive_map_optimization PUBLIC
        OpenMeshCore
        OpenMeshTools
        CoMISo
        PRIVATE Eigen3::Eigen)


#make library for openflipper
add_library(Adopt STATIC
    src/adaptive_map.cc
    src/helpers.cc
    src/adaptive_map.hh
    src/helpers.hh
    src/includes.hh
    src/OptimizationElements2D.hh
    src/RemesherBase.hh
    src/RemesherValentin.hh
    src/RemesherValentin.cc
    src/RemesherNaive.hh
    src/RemesherNaive.cc
)



target_include_directories(Adopt
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
       $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
       $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
       $<BUILD_INTERFACE:${eigen3_SOURCE_DIR}>
       ${EXTERNAL_DIR}/OpenMesh/src/
        ${EXTERNAL_DIR}/CoMISo/src/
        ${EXTERNAL_DIR}/TinyAD/src/

    )

target_link_libraries (Adopt PRIVATE
    OpenMeshCore
    OpenMeshTools
    CoMISo
    Eigen3::Eigen)

if(NOT "${CMAKE_PROJECT_NAME}" MATCHES "Adopt")
    set (ADOPT_FOUND true PARENT_SCOPE)
    set (ADOPT_LIBRARIES Adopt PARENT_SCOPE)
    set (ADOPT_LIBRARY Adopt PARENT_SCOPE)
    set (ADOPT_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src/" PARENT_SCOPE)
    set (ADOPT_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE PATH "The directory where the Adopt libraries can be found.")
endif()

