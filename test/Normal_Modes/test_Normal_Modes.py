import numpy as np 
import os
import subprocess
import pytest 

@pytest.fixture
def vibrant_binary_path():
    return "../../vibrant"

def run_vibrant(binary_path: str, input_file: str, omp_threads: int = 1):
    """Run vibrant using the specified input file as a command-line argument"""
    my_env = {**os.environ, "OMP_NUM_THREADS": str(omp_threads)}
    calc = subprocess.run([binary_path, input_file],
                          env=my_env,
                          text=True,
                          capture_output=True)
    return calc.returncode, calc.stdout

def parse_normal_freq_from_file(fname:str):
    """ parse the generated normal mode frequencies """
    data = np.loadtxt(fname)
    x = data[:]
    return x
    
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
    reference1 = "output/normal_mode_freq.txt"
    test1 = "normal_mode_freq.txt"
    
    reference2 = "output/normal_mode_displ.txt" 
    test2 = "normal_mode_displ.txt"

    # parse filenames
    x_ref1 = parse_normal_freq_from_file(reference1)
    x_test1 = parse_normal_freq_from_file(test1)

    x_ref2, y_ref2, z_ref2 = parse_normal_displ_from_file(reference2)
    x_test2, y_test2, z_test2 = parse_normal_displ_from_file(test2)
    
    # compare test against reference
    assert np.allclose(x_ref1, x_test1, atol=1e-8)
    assert np.allclose(np.abs(x_ref2), np.abs(x_test2), atol=1e-8)
    assert np.allclose(np.abs(y_ref2), np.abs(y_test2), atol=1e-8)
    assert np.allclose(np.abs(z_ref2), np.abs(z_test2), atol=1e-8)
