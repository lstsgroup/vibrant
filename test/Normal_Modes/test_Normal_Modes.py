import numpy as np 
import os
import subprocess
import pytest 

@pytest.fixture
def vibrant_binary_path():
    return "../../vibrant"

def run_vibrant(binary_path:str, input_file:str, omp_threads:int=1):
    """ run vibrant using the specified input file """
    # read input file
    with open(input_file) as f:
        in_data = ""
        for line in f.readlines():
            in_data +=line
    # run vibrant calculation
    my_env = {**os.environ, "OMP_NUM_THREADS": str(omp_threads)}
    calc = subprocess.run([binary_path], 
                          env = my_env,
                          input=in_data, 
                          text=True, 
                          capture_output=True
    )
    return calc.returncode, calc.stdout

    
def parse_normal_displ_from_file(fname:str):
    """ parse the generated normal mode displacements """
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    return x, y, z

def test_Normal_Modes(vibrant_binary_path):
    """ run test case normal modes """
    # run vibrant calculation and test for successfull exit
    returncode, stdout = run_vibrant(vibrant_binary_path, "input.txt")
    assert returncode == 0
    
    # file names
    
    reference = "output/normal_mode_displ.txt" 
    test = "normal_mode_displ.txt"

    # parse filenames

    x_ref, y_ref, z_ref = parse_normal_displ_from_file(reference)
    x_test, y_test, z_test = parse_normal_displ_from_file(test)
    
    # compare test against reference
    assert np.allclose(np.abs(x_ref), np.abs(x_test), atol=1e-8)
    assert np.allclose(np.abs(y_ref), np.abs(y_test), atol=1e-8)
    assert np.allclose(np.abs(z_ref), np.abs(z_test), atol=1e-8)
