cmake_minimum_required(VERSION 3.24)
project(microBUDE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(microBUDE main.cpp)

find_package(OpenMP REQUIRED)

#target_compile_options(microBUDE PRIVATE -fsanitize=address)
#target_link_options(microBUDE PRIVATE -fsanitize=address)

target_link_libraries(microBUDE PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C)

target_compile_options(microBUDE PRIVATE "$<$<CONFIG:RelWithDebInfo>:-Ofast;-march=native>")
target_compile_options(microBUDE PRIVATE "$<$<CONFIG:Release>:-Ofast;-march=native>")
target_compile_options(microBUDE PRIVATE "$<$<CONFIG:Debug>:-O2;-fno-omit-frame-pointer>")
