cmake_minimum_required(VERSION 3.15)

project(
  Edit-Distance
  VERSION 1.0
  DESCRIPTION "Edit Distance"
  HOMEPAGE_URL "https://github.com/ucrparlay/Edit-Distance"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_SANITIZER "Enable AddressSanitizer" OFF)

#Til OMP 
#set(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
#set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")

#find_package(OpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib")

# Sanitizer configuration
if(ENABLE_SANITIZER)
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=address)
endif()

option(CILKPLUS "" OFF)
option(OPENCILK "" OFF)
option(SERIAL "" OFF)
option(STDALLOC "" OFF)

if(CILKPLUS)
  add_compile_options(-DPARLAY_CILKPLUS -DCILK -fcilkplus)
elseif(OPENCILK)
  add_compile_options(-DPARLAY_OPENCILK -DCILK -fopencilk)
elseif(SERIAL)
  add_compile_options(-DPARLAY_SEQUENTIAL)
else()
  add_compile_options(-pthread)
endif()

add_compile_options(-O3 -mcx16 -march=native)


if(STDALLOC)
  add_compile_options(-DPARLAY_USE_STD_ALLOC)
endif()

add_subdirectory(parlaylib)

add_executable(test_framework test_framework.cpp)
target_link_libraries(test_framework PUBLIC edit_distance_parallel parlay)

add_executable(test_framework_real test_framework_real.cpp)
target_link_libraries(test_framework_real PUBLIC edit_distance_parallel parlay)

add_executable(test_framework_ratio test_framework_ratio.cpp)
target_link_libraries(test_framework_ratio PUBLIC edit_distance_parallel parlay)

add_library(edit_distance_parallel edit_distance_parallel.cpp)
target_link_libraries(edit_distance_parallel PUBLIC parlay)

add_executable(suffix_array_test suffix_array_test.cpp)
target_link_libraries(suffix_array_test PUBLIC suffix_array_sequential parlay)

add_library(suffix_array_sequential suffix_array_sequential.cpp)
