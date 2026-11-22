cmake_minimum_required(VERSION 3.10)
project(StressSolver C)

set(CMAKE_C_STANDARD 11)
add_compile_options(-Wall -Wextra -Wpedantic)

# Create shared library
add_library(stress_lib SHARED stress.c)
target_include_directories(stress_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(stress_lib m) # Link math library
set_target_properties(stress_lib PROPERTIES OUTPUT_NAME stress)

# Create executable
add_executable(stress main.c)
target_link_libraries(stress stress_lib m)

# Installation options (optional)
install(TARGETS stress_lib stress
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(FILES stress.h DESTINATION include)

enable_testing()
add_executable(test_stress tests/test_stress.c)
target_link_libraries(test_stress stress_lib m)
add_test(NAME test_stress COMMAND test_stress)
