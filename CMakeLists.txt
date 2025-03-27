cmake_minimum_required(VERSION 3.15)
project(Angles)
file(READ ${CMAKE_SOURCE_DIR}/VERSION CMAKE_PROJECT_VERSION) # set project version
add_compile_definitions(EXEC_VERSION="${CMAKE_PROJECT_VERSION}") # set version string for executables

set(CMAKE_CXX_STANDARD 17) # C++ std filesystem

include(FetchContent)
FetchContent_Declare(cmake_common GIT_REPOSITORY https://gitlab.au.dk/ausa/cmake_common.git)
FetchContent_MakeAvailable(cmake_common) # custom CMake scripts commonly used across AUSA projects


find_package(ROOT REQUIRED)
find_package(GSL REQUIRED)
find_package(AUSALIB REQUIRED)
find_package(cxxopts REQUIRED)

add_executable(Angles Angles.cpp)
target_link_libraries(Angles PRIVATE AUSA::AUSA)
target_include_directories(Angles PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Make the install targets
install(TARGETS Angles
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT bin)
# Install headers and resources as well, preserving directory hierarchy
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ # use https://cmake.org/cmake/help/latest/command/target_sources.html#file-sets when CMake version >= 3.23 can be expected
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Angles
    COMPONENT include)

# Report components to be installed to user
list(APPEND INSTALL_COMPONENTS bin include)
list(APPEND INSTALL_PATHS ${CMAKE_INSTALL_BINDIR} ${CMAKE_INSTALL_INCLUDEDIR}/Angles)
