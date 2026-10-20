set(TEST_FILES
    testGrid.cpp
    testFvMesh.cpp
    testMatrices.cpp
    testScalarFields.cpp
    testLaplacian.cpp
    testGrad.cpp
    )

foreach(test_source IN LISTS TEST_FILES)

    get_filename_component(test_name ${test_source} NAME_WE)

    add_executable(${test_name} ${test_source})

    target_link_libraries(${test_name} PRIVATE NastsCore)

    add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    
endforeach()