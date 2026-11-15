cmake_minimum_required(VERSION 3.16)

project(
  Intergenic
  VERSION 0.0.2
  DESCRIPTION "Format-aware extraction of intergenic regions"
  LANGUAGES CXX
)

include(GNUInstallDirs)

option(INTERGENIC_USE_ZLIB "Enable gzip-compressed input" OFF)

add_executable(
  intergenic
  src/Intergenic.cpp
  src/Classes.cpp
  src/Motif.cpp
  src/Utils.cpp
)

target_compile_features(intergenic PRIVATE cxx_std_17)
set_target_properties(
  intergenic
  PROPERTIES
    CXX_EXTENSIONS OFF
    OUTPUT_NAME intergenic
)

target_include_directories(
  intergenic
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

if(MSVC)
  target_compile_options(intergenic PRIVATE /W4 /permissive-)
else()
  target_compile_options(intergenic PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(INTERGENIC_USE_ZLIB)
  find_package(ZLIB REQUIRED)
  target_compile_definitions(intergenic PRIVATE ZLIB)
  target_link_libraries(intergenic PRIVATE ZLIB::ZLIB)
endif()

install(TARGETS intergenic RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(
  FILES README.md LICENSE COPYING
  DESTINATION "${CMAKE_INSTALL_DOCDIR}"
)
install(
  DIRECTORY docs/
  DESTINATION "${CMAKE_INSTALL_DOCDIR}/docs"
  FILES_MATCHING PATTERN "*.md"
)

set(CPACK_PACKAGE_NAME "intergenic")
set(CPACK_PACKAGE_VENDOR "TUBIC, Tianjin University")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(
  CPACK_PACKAGE_FILE_NAME
  "intergenic-${PROJECT_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}"
)
set(CPACK_GENERATOR "TGZ;ZIP")
include(CPack)
