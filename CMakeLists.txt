project(UnscentedKF)

cmake_minimum_required(VERSION 3.5)

#FILE(GLOB Eigen_plugins_SRCS "*.h")
#
#INSTALL(FILES
#  ${Eigen_plugins_SRCS}
#  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/plugins COMPONENT Devel
#  )

add_definitions(-std=c++0x)

set(sources
        src/ukf.cpp
        src/main.cpp
        src/tools.cpp)

add_executable(UnscentedKF ${sources})