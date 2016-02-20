include_directories( ${INCLUDES} )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
add_executable( simulation_loop_model_single_shape
    simulation1.cpp
    basalganglia.cpp
    csntc.cpp
    layer.cpp
    md5.cpp
    reservoir.cpp
    simulator.cpp
    utils.cpp
    )

set_target_properties(simulation_loop_model_single_shape 
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/simulation_loop_model_single_shape )

target_link_libraries ( simulation_loop_model_single_shape ${LIBS} )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
add_executable( simulation_loop_model_multiple_shapes
    simulation2.cpp
    basalganglia.cpp
    csntc.cpp
    layer.cpp
    md5.cpp
    reservoir.cpp
    simulator.cpp
    utils.cpp
    )


set_target_properties( simulation_loop_model_multiple_shapes 
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/simulation_loop_model_multiple_shapes )

target_link_libraries ( simulation_loop_model_multiple_shapes ${LIBS} )
            

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
add_executable( bg_test 
    basalganglia.cpp
    basalganglia_test.cpp
    layer.cpp
    utils.cpp
    )

set_target_properties( bg_test 
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/basalganglia )

target_link_libraries ( bg_test ${LIBS} )
            

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
add_executable( res_test 
    reservoir.cpp
    reservoir_test.cpp
    md5.cpp
    utils.cpp
    )

set_target_properties( res_test 
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/reservoir )

target_link_libraries ( res_test ${LIBS} )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
if(WITH_CENSLIB)
    add_executable( simulation_loop_model_hand
        hand.cpp
        hand_engine.cpp
        hand_controller.cpp
        basalganglia.cpp
        csntc.cpp
        layer.cpp
        md5.cpp
        reservoir.cpp
        simulator.cpp
        utils.cpp
        )


    set_target_properties( simulation_loop_model_hand 
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/simulation_loop_model_hand )

    target_link_libraries ( simulation_loop_model_hand ${LIBS} )
endif()
