# =============================================================================
# CMakeLists.txt
#
# Configure application library.
# -----------------------------------------------------------------------------
#
# author : Marco Capettini <s252912@studenti.polito.it>
#          Emanuele Parisi <emanuele.parisi@polito.it>
#          Gianvito Urgese <gianvito.urgese@polito.it>
# date   : October, 2019
# =============================================================================

include(FetchContent)

# -----------------------------------------------------------------------------
# Retrieve SeqAn3 library
# -----------------------------------------------------------------------------
#FetchContent_Declare(
#  seqan3
#  GIT_REPOSITORY https://github.com/seqan/seqan3.git
#  GIT_TAG 3.0.0
#)

#FetchContent_MakeAvailable(seqan3)
#find_package(SeqAn3
#             REQUIRED
#             HINTS
#             ${seqan3_SOURCE_DIR}/build_system
set(SEQAN3_NO_SDSL 1)
find_package(seqan3 REQUIRED)

# -----------------------------------------------------------------------------
# Check the presence of libseqan2-dev
# -----------------------------------------------------------------------------
find_package(SeqAn)
IF (${SeqAn_FOUND})
   MESSAGE(STATUS "Found libseqan2-dev.")
ELSE (${SeqAn_FOUND})
   MESSAGE(STATUS "WARNING: Could not locate libseqan2-dev, please refer to https://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html.")
ENDIF (${SeqAn_FOUND})

# -----------------------------------------------------------------------------
# Check the thread library used by the system
# -----------------------------------------------------------------------------
find_package(Threads
             REQUIRED)

# -----------------------------------------------------------------------------
# Check the presence of OpenMP
# -----------------------------------------------------------------------------
find_package(OpenMP)

# -----------------------------------------------------------------------------
# Create the isomiR-SEA interface library
# -----------------------------------------------------------------------------
add_library(isomiR_SEA_isomiR_SEA
            INTERFACE)
target_include_directories(isomiR_SEA_isomiR_SEA
                           INTERFACE
                           ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(isomiR_SEA_isomiR_SEA
                      INTERFACE
                      seqan3::seqan3
                      Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(isomiR_SEA_isomiR_SEA INTERFACE OpenMP::OpenMP_CXX)
endif()
add_library(isomiR_SEA::isomiR_SEA
            ALIAS
            isomiR_SEA_isomiR_SEA)
