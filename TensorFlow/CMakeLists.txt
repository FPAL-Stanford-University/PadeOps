set(tf_load_src_files 
    #tfloadscript_lstm.cc
    tfloadscript_unet.cc
)

add_library(tf_load ${tf_load_src_files})
#set_target_properties(tf_load PROPERTIES LINKER_LANGUAGE Fortran) 

target_include_directories(tf_load PUBLIC
    ${TensorFlow_INCLUDE_DIRS}
    ${TensorFlow_INCLUDE_DIRS}/src
)
message(STATUS "TensorFlow_INCLUDE_DIRS=${TensorFlow_INCLUDE_DIRS}")

target_link_libraries(tf_load PUBLIC
    ${TensorFlow_LIBRARIES}
    ${TensorFlow_CC_LIB}
    ${TensorFlow_FRAMEWORK_LIB}
)
