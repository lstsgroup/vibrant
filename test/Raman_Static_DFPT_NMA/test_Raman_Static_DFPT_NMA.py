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

def parse_spectrum_from_file(fname:str):
    """ parse the generated spectrum """
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]
    return x, y

def test_Raman_Static_DFPT_NMA(vibrant_binary_path):
    """ run test case """
    # run vibrant calculation and test for successfull exit
    returncode, stdout = run_vibrant(vibrant_binary_path, "input.txt")
    assert returncode == 0
    
    # file names
    reference1 = "output/normal_mode_freq.txt"
    test1 = "normal_mode_freq.txt"
    
    reference2 = "output/normal_mode_displ.txt" 
    test2 = "normal_mode_displ.txt"
    
    reference3 = "output/result_static_raman.txt"
    test3 = "result_static_raman.txt"
    
    # parse filenames
    x_ref1 = parse_normal_freq_from_file(reference1)
    x_test1 = parse_normal_freq_from_file(test1)

    x_ref2, y_ref2, z_ref2 = parse_normal_displ_from_file(reference2)
    x_test2, y_test2, z_test2 = parse_normal_displ_from_file(test2)
    
    x_ref3, y_ref3 = parse_spectrum_from_file(reference3)
    x_test3, y_test3 = parse_spectrum_from_file(test3)

    # compare test against reference
    assert np.allclose(x_ref1, x_test1, atol=1e-8)
    assert np.allclose(np.abs(x_ref2), np.abs(x_test2), atol=1e-8)
    assert np.allclose(np.abs(y_ref2), np.abs(y_test2), atol=1e-8)
    assert np.allclose(np.abs(z_ref2), np.abs(z_test2), atol=1e-8)
    assert np.allclose(x_ref3, x_test3, atol=1e-8)
    assert np.allclose(y_ref3, y_test3, atol=1e-8)
