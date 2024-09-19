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

    

def parse_spectrum_from_file(fname:str):
    """ parse the generated spectrum """
    data = np.loadtxt(fname)
    x = data[:, 0]
    y = data[:, 1]
    return x, y


def test_Raman_Berry(vibrant_binary_path):
    """ run test case """
    # run vibrant calculation and test for successfull exit
    returncode, stdout = run_vibrant(vibrant_binary_path, "input.txt")
    assert returncode == 0
    
    # file names
    reference = "output/COF-1_Raman_sc_Berry.txt"
    test = "result_fft_water_lib_unpol.txt"
    
    # parse filenames
    x_ref, y_ref = parse_spectrum_from_file(reference)
    x_test, y_test = parse_spectrum_from_file(test)

    # compare test against reference
    assert np.allclose(x_ref, x_test, atol=1e-8)
    assert np.allclose(y_ref, y_test, atol=1e-8)

    