# =============================================================================
# CMakeLists.txt
#
# Configure application suite.
# -----------------------------------------------------------------------------
#
# author : Marco Capettini <s252912@studenti.polito.it>
#          Emanuele Parisi <emanuele.parisi@polito.it>
#          Gianvito Urgese <gianvito.urgese@polito.it>
# date   : October, 2019
# =============================================================================

# -----------------------------------------------------------------------------
# Configure isomiR-SEA targets
# -----------------------------------------------------------------------------
add_executable(isomir_sea
               ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(isomir_sea
                      PUBLIC
                      isomiR_SEA::isomiR_SEA stdc++fs)

install(TARGETS isomir_sea DESTINATION bin)
