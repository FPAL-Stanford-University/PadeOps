# This cmake module configures find_package(Pytorch)
# we note this is NOT needed for torch c++ and pytorch geometric stuff for now!

if (DEFINED $ENV{Pytorch_ROOT})
    set(Pytorch_ROOT $ENV{Pytorch_ROOT})
endif()

list(APPEND CMAKE_PREFIX_PATH "${Pytorch_ROOT}/share/cmake/Torch")
#list(APPEND CMAKE_PREFIX_PATH "${Pytorch_ROOT}")

find_package(Torch REQUIRED)
message(STATUS "CMAKE Torch FOUND")

# init comment
find_package(PkgConfig)

#if (PKG_CONFIG_FOUND AND NOT Pytorch_ROOT)
#    pkg_check_modules(PKG_TensorFlow QUIET "tensorflow_cc" "tensorflow_framework")
#endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES})

#find_path( TensorFlow_INCLUDE_DIRS
#    NAMES "farmhash.h"
#    PATHS ${TensorFlow_ROOT} ${PKG_TensorFlow_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
#    PATH_SUFFIXES "include" "inc"
#) 

#find_library( TensorFlow_CC_LIB
#    NAMES "tensorflow_cc" "libtensorflow_cc.so.2"
#    PATHS ${TensorFlow_ROOT} ${PKG_TensorFlow_LIBRARY_DIRS} ${INCLUDE_INSTALL_DIR}
#)

#find_library( TensorFlow_FRAMEWORK_LIB
#    NAMES "tensorflow_framework" "libtensorflow_framework.so.2"
#    PATHS ${TensorFlow_ROOT} ${PKG_TensorFlow_LIBRARY_DIRS} ${INCLUDE_INSTALL_DIR}
#)

#if (TensorFlow_CC_LIB)
#   set(TensorFlow_CC_LIB_FOUND TRUE)
#    set(TensorFlow_LIBRARIES ${TensorFlow_LIBRARIES} ${TensorFlow_CC_LIB})
#    add_library(TensorFlow::cc INTERFACE IMPORTED)
#    set_target_properties(TensorFlow::cc PROPERTIES
#        INTERFACE_INCLUDE_DIRECTORIES "${TensorFlow_INCLUDE_DIRS}"
#        INTERFACE_LINK_LIBRARIES "${TensorFlow_CC_LIB}")
#    message(STATUS "Using TensorFlow_CC_LIB: ${TensorFlow_CC_LIB}")
#else()
#    set(TensorFlow_CC_LIB_FOUND FALSE)
#    message(STATUS "${TensorFlow_CC_LIB}")
#endif()


#if (TensorFlow_FRAMEWORK_LIB)
#    set(TensorFlow_FRAMEWORK_LIB_FOUND TRUE)
#    set(TensorFlow_LIBRARIES ${TensorFlow_LIBRARIES} ${TensorFlow_FRAMEWORK_LIB})
#    add_library(TensorFlow::framework INTERFACE IMPORTED)
#    set_target_properties(TensorFlow::framework PROPERTIES
#        INTERFACE_INCLUDE_DIRECTORIES "${TensorFlow_INCLUDE_DIRS}"
#        INTERFACE_LINK_LIBRARIES "${TensorFlow_FRAMEWORK_LIB}")
#    message(STATUS "Using TensorFlow_FRAMEWORK_LIB: ${TensorFlow_FRAMEWORK_LIB}")
#else()
#    set(TensorFlow_FRAMEWORK_LIB_FOUND FALSE)
#    message(STATUS "${TensorFlow_FRAMEWORK_LIB}")
#endif()

# the two below commented
set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV})

include(FindPackageHandleStandardArgs)

#find_package_handle_standard_args(TensorFlow
#    REQUIRED_VARS
#    TensorFlow_INCLUDE_DIRS
#    HANDLE_COMPONENTS
#)

#mark_as_advanced(
#    TensorFlow_INCLUDE_DIRS
#    TensorFlow_LIBRARIES
#    TensorFlow_CC_LIB
#    TensorFlow_FRAMEWORK_LIB
#)
