cmake_minimum_required(VERSION 3.16)
project(exubercore VERSION 0.3.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Armadillo REQUIRED)

add_library(exubercore src/radf.cpp)
add_library(exubercore::exubercore ALIAS exubercore)

target_include_directories(exubercore PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${ARMADILLO_INCLUDE_DIRS}
)
target_link_libraries(exubercore PUBLIC ${ARMADILLO_LIBRARIES})
target_compile_features(exubercore PUBLIC cxx_std_17)

option(EXUBERCORE_BUILD_TESTS "Build exubercore's test suite" ON)
if(EXUBERCORE_BUILD_TESTS)
  enable_testing()
  add_executable(test_radf tests/test_radf.cpp)
  target_link_libraries(test_radf PRIVATE exubercore::exubercore)
  add_test(NAME test_radf
    COMMAND test_radf ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/golden)
endif()

include(GNUInstallDirs)
install(TARGETS exubercore EXPORT exubercoreTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/exubercoreConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/exubercoreConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/exubercoreConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/exubercore
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/exubercoreConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/exubercoreConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/exubercore
)
install(EXPORT exubercoreTargets
  FILE exubercoreTargets.cmake
  NAMESPACE exubercore::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/exubercore
)
