cmake_minimum_required(VERSION 3.13)

project(gradient LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_CURRENT_LIST_DIR}/lib)

set(GUI_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

file(GLOB_RECURSE SOURCE_FILES
    ${CMAKE_CURRENT_LIST_DIR}/*.cpp
    ${CMAKE_CURRENT_LIST_DIR}/*.c
    ${GUI_DIR}/*.cpp
    ${GUI_DIR}/*.c
    ${GUI_DIR}/*.inl)

file(GLOB_RECURSE HEADER_FILES
    ${CMAKE_CURRENT_LIST_DIR}/*.hpp
    ${CMAKE_CURRENT_LIST_DIR}/*.h
    ${GUI_DIR}/*.hpp
    ${GUI_DIR}/*.h
    ${GUI_DIR}/*.inl)

include_directories(${CMAKE_CURRENT_LIST_DIR}/lib)
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

find_package(SFML REQUIRED COMPONENTS graphics audio window system)
find_library(SFML-WIDGETS_LIB libsfml-widgets.a)
LINK_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/lib)

if(UNIX AND NOT APPLE)
    #GNU/Linux specific code
    target_link_libraries(gradient PUBLIC sfml-graphics sfml-audio sfml-window sfml-system GL ${SFML-WIDGETS_LIB})
elseif(APPLE)
    #OSX specific code
    target_link_libraries(gradient PUBLIC sfml-graphics sfml-audio sfml-window sfml-system "-framework OpenGL" ${SFML-WIDGETS_LIB})
endif()
