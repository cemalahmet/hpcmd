# Defining the executable target
add_executable(06 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(06 PUBLIC my_md_lib)

#add_input_file(06 lj54.xyz)

add_compile_options(-Wall)