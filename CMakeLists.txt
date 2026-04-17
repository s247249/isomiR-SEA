# =============================================================================
# CMakeLists.txt
#
# Top-level CMake file setting up project configuration steps.
# -----------------------------------------------------------------------------
#
# author : Marco Capettini <s252912@studenti.polito.it>
#          Emanuele Parisi <emanuele.parisi@polito.it>
#          Gianvito Urgese <gianvito.urgese@polito.it>
# date   : October, 2019
# =============================================================================

cmake_minimum_required(VERSION
                       3.5)

# -----------------------------------------------------------------------------
# Initialize project metadata
# -----------------------------------------------------------------------------
project(isomiR-SEA
        VERSION
        2.0.0.0
        DESCRIPTION
        ""
        HOMEPAGE_URL
        ""
        LANGUAGES
        CXX)

# -----------------------------------------------------------------------------
# Make sure the project is configured for release, if not differently 
# specified
# -----------------------------------------------------------------------------
set(${CMAKE_PROJECT_NAME}_DEFAULT_BUILD_TYPE
    "Release")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        ${${CMAKE_PROJECT_NAME}_DEFAULT_BUILD_TYPE}
        CACHE
        STRING
        "Define build type"
        FORCE)
endif()
message(STATUS
        "${CMAKE_PROJECT_NAME} build type is ${CMAKE_BUILD_TYPE}")

# -----------------------------------------------------------------------------
# Configure isomiR-SEA library and applications
# -----------------------------------------------------------------------------
message(STATUS
        "Configure ${CMAKE_PROJECT_NAME} application library")
add_subdirectory(${CMAKE_SOURCE_DIR}/include)

message(STATUS
        "Configure ${CMAKE_PROJECT_NAME} applications")
add_subdirectory(${CMAKE_SOURCE_DIR}/apps)

# -----------------------------------------------------------------------------
# Configure isomiR-SEA tests
# -----------------------------------------------------------------------------
option(${CMAKE_PROJECT_NAME}_CONFIGURE_TESTS
       "Configure ${CMAKE_PROJECT_NAME} tests"
       NO)
if(${${CMAKE_PROJECT_NAME}_CONFIGURE_TESTS})
    message(STATUS
            "Configure ${CMAKE_PROJECT_NAME} tests")
    add_subdirectory(${CMAKE_SOURCE_DIR}/tests)
endif()
