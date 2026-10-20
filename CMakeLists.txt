cmake_minimum_required(VERSION 3.10)
project(NASTS LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB_RECURSE SOURCES "src/*.cpp")
add_library(NastsCore STATIC ${SOURCES})
target_include_directories(NastsCore PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(nasts main.cpp)
target_link_libraries(nasts PRIVATE NastsCore)

enable_testing()
add_subdirectory(tests)