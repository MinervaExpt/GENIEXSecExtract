#Find the source and header files.  This isn't great practice because you have to
#rerun cmake to detect new files to link.
add_library(${PROJECT_NAME} SHARED TreeWrapper.cxx ChainWrapper.cxx XSec.cxx XSecLooper.cxx)
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} MAT MAT-MINERvA)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets DESTINATION lib)
install(EXPORT ${PROJECT_NAME}Targets DESTINATION lib/cmake/${PROJECT_NAME})

export(TARGETS ${PROJECT_NAME} FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)
