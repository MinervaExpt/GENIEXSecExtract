add_executable(runXSecLooper_CCQENuInclusive runXSecLooper_CCQENuInclusive.cpp)
target_link_libraries(runXSecLooper_CCQENuInclusive ${ROOT_LIBRARIES} MAT MAT-MINERvA GENIEXSecExtract)
install(TARGETS runXSecLooper_CCQENuInclusive DESTINATION bin)

add_executable(runCCIncForMultiNeutron runCCIncForMultiNeutron.cpp)
target_link_libraries(runCCIncForMultiNeutron ${ROOT_LIBRARIES} MAT MAT-MINERvA GENIEXSecExtract GenVector)
install(TARGETS runCCIncForMultiNeutron DESTINATION bin)

add_executable(runCCIncForNuEMEC runCCIncForNuEMEC.cpp)
target_link_libraries(runCCIncForNuEMEC ${ROOT_LIBRARIES} MAT MAT-MINERvA GENIEXSecExtract)
install(TARGETS runCCIncForNuEMEC DESTINATION bin)
