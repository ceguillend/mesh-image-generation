#http://www.cmake.org/cmake/help/cmake2.6docs.html
cmake_minimum_required (VERSION 2.6)
project (Method)

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")

FIND_PACKAGE( OpenCV REQUIRED )


ADD_EXECUTABLE(Method main.cpp delaunay.cpp basic_geo.cpp curve.cpp thinning.cpp color_utils.cpp) 
TARGET_LINK_LIBRARIES(Method ${OpenCV_LIBS} )
