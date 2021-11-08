


#--- Targets ---

set(BINARY ${PROJECT_NAME})

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp *.hpp)
set(SOURCES ${SOURCES})

#add_executable(${BINARY} ${SOURCES})


#target_link_libraries(${BINARY} ${LIBS_BPP})
#target_link_libraries(${BINARY} glog::glog)

#add_executable(${BINARY}_run ${SOURCES})

add_library(${BINARY}_lib STATIC ${SOURCES})


MESSAGE(STATUS "SRC Targets OK")
