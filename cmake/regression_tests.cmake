# Python required for application testing
find_package(Python3 3.10 COMPONENTS Interpreter Development)
if (Python3_FOUND)
    message("-- Python 3 interpreter version: " ${Python3_VERSION})
else()
    message("-- Python 3 interpreter not found")
endif()

# Enable ctest
enable_testing()

# copy the test the test directory to build directory
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/test
            ${PROJECT_BINARY_DIR}/test)

# adding a test to ctest
add_test(
    NAME IR_Berry
    COMMAND pytest test_IR_Berry.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/IR_Berry
)
add_test(
    NAME IR_Wannier_Ph
    COMMAND pytest test_IR_Wannier_Ph.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/IR_Wannier_Ph
)
add_test(
    NAME IR_Wannier_whole
    COMMAND pytest test_IR_Wannier_whole.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/IR_Wannier_whole
)
add_test(
    NAME Raman_Berry
    COMMAND pytest test_Raman_Berry.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Raman_Berry
)
add_test(
    NAME Raman_DFPT
    COMMAND pytest test_Raman_DFPT.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Raman_DFPT 
)
add_test(
    NAME Power_pos_mw
    COMMAND pytest test_Power_pos_mv.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Power/pos_mw
)
add_test(
    NAME Power_pos_nomw
    COMMAND pytest test_Power_pos_nomw.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Power/pos_nomw
)
add_test(
    NAME Power_vel_mw
    COMMAND pytest test_Power_vel_mw.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Power/vel_mw
)
add_test(
    NAME Power_vel_nomw
    COMMAND pytest test_Power_vel_nomw.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Power/vel_nomw
)
add_test(
    NAME Absorption
    COMMAND pytest test_Absorption.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Absorption/normal
)
add_test(
    NAME Absorption_Pade
    COMMAND pytest test_Absorption_Pade.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Absorption/Pade
)
add_test(
    NAME Normal_Modes
    COMMAND pytest test_Normal_Modes.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Normal_Modes
)
add_test(
    NAME IR_Static_Ext
    COMMAND pytest test_IR_Static_Ext.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/IR_Static_Ext
)
add_test(
    NAME IR_Static_NMA
    COMMAND pytest test_IR_Static_NMA.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/IR_Static_NMA
)
add_test(
    NAME Raman_Static_DFPT_Ext
    COMMAND pytest test_Raman_Static_DFPT_Ext.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Raman_Static_DFPT_Ext
)
add_test(
    NAME Raman_Static_DFPT_NMA
    COMMAND pytest test_Raman_Static_DFPT_NMA.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Raman_Static_DFPT_NMA
)
add_test(
    NAME RR_Static_normal_Ext
    COMMAND pytest test_RR_Static_normal_Ext.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Resonance_Raman_Static_Ext/normal
)
add_test(
    NAME RR_Static_Pade_Ext
    COMMAND pytest test_RR_Static_Pade_Ext.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Resonance_Raman_Static_Ext/Pade
)
add_test(
    NAME RR_Static_normal_NMA
    COMMAND pytest test_RR_Static_normal_NMA.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Resonance_Raman_Static_NMA/normal
)
add_test(
    NAME RR_Static_Pade_NMA
    COMMAND pytest test_RR_Static_Pade_NMA.py
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/Resonance_Raman_Static_NMA/Pade
)
