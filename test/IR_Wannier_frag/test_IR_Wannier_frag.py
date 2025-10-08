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



def parse_ir_spectrum_from_file(fname:str, hdr:int=0):
    """ parse the generated IR spectrum """
    data = np.loadtxt(fname, skiprows=hdr)
    x = data[:, 0]
    y = data[:, 1]
    return x, y


def test_IR_Wannier_whole(vibrant_binary_path):
    """ run test case IR Berry """
    # run vibrant calculation and test for successfull exit
    returncode, stdout = run_vibrant(vibrant_binary_path, "input.txt")
    assert returncode == 0
    
    # file names
    reference1 = "output/IR_spectrum_fragment_1.txt"
    test1 = "IR_spectrum_fragment_1.txt"
    
    reference2 = "output/IR_spectrum_fragment_2.txt" 
    test2 = "IR_spectrum_fragment_2.txt"

    # parse filenames
    x_ref1, y_ref1 = parse_ir_spectrum_from_file(reference1)
    x_test1, y_test1 = parse_ir_spectrum_from_file(test1, hdr=1)

    x_ref2, y_ref2 = parse_ir_spectrum_from_file(reference2)
    x_test2, y_test2 = parse_ir_spectrum_from_file(test2, hdr=1)

    # compare test against reference
    assert np.allclose(x_ref1, x_test1, atol=1e-8)
    assert np.allclose(y_ref1, y_test1, atol=1e-8)
    assert np.allclose(x_ref2, x_test2, atol=1e-8)
    assert np.allclose(y_ref2, y_test2, atol=1e-8)

