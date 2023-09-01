# Defining the executable target
add_executable(05 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(05 PUBLIC my_md_lib)

# Copying "input" file simulation_test_input.xyz to the same directory as executable
#add_input_file(05 simulation_test_input.txt)

add_compile_options(-Wall)