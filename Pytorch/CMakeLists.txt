set(pytorch_load_src_files
     pytorchscript_interface.cc
     pytorchscript_interface_CUDA.cc
)

add_library(pytorch_load OBJECT ${pytorch_load_src_files})
#set_target_properties(pytorch_load PROPERTIES LINKER_LANGUAGE Fortran)


#find_package(Torch REQUIRED)
#message(STATUS "CMAKE Torch FOUND")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

target_include_directories(pytorch_load PUBLIC
  ${TORCH_INCLUDE_DIRS}
)

message(STATUS "Using Torch_INCLUDE_DIRS: ${TORCH_INCLUDE_DIRS}")
#set(GLIBC_DIR "/ocean/projects/phy230038p/awu2/glibc_install/glibc-2.29")
#include_directories(${GLIBC_DIR}/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -Wl,--rpath=${GLIBC_DIR}/lib")
#set(CMAKE_EXE_LINKER_FLAGS "-L${GLIBC_DIR}/lib -Wl,--dynamic-linker=${GLIBC_DIR}/lib/ld-linux-x86-64.so.2")


target_link_libraries(pytorch_load PUBLIC ${TORCH_LIBRARIES})
set_property(TARGET pytorch_load PROPERTY CXX_STANDARD 17)
message(STATUS "Using Torch_LIBRARIES: ${TORCH_LIBRARIES}")
#execute_process(
#        COMMAND ${Python_EXECUTABLE} -c "import torch; print(torch.utils.cmake_prefix_path)"
#        OUTTPUT_VARIABLE TORCH_CMAKE_PREFIX_PATH
#        OUTPUT_STRIP_TRAILING_WHITESPACE
#)
#set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${TORCH_CMAKE_PREFIX_PATH})
#find_package(torch REQUIRED)

#target_include_directories(pytorch_load PUBLIC ${Python3_INCLUDE_DIRS})
#target_link_libraries(pytorch_load PUBLIC ${Python3_LIBRARIES} pybind11::embed)
message("CMAKE Torch LINKED")

