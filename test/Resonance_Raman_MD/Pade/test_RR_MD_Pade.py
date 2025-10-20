import numpy as np 
import os
import subprocess
import pytest 

@pytest.fixture
def vibrant_binary_path():
    return "../../../vibrant"

def run_vibrant(binary_path: str, input_file: str, omp_threads: int = 1):
    """Run vibrant using the specified input file as a command-line argument"""
    my_env = {**os.environ, "OMP_NUM_THREADS": str(omp_threads)}
    calc = subprocess.run([binary_path, input_file],
                          env=my_env,
                          text=True,
                          capture_output=True)
    return calc.returncode, calc.stdout


def parse_spectrum_from_file(fname:str, hdr:int=0):
    """ parse the generated spectrum """
    data = np.loadtxt(fname, skiprows=hdr)
    x = data[:, 0]
    y = data[:, 1]
    return x, y



def test_Raman_Berry(vibrant_binary_path):
    """ run test case """
    # run vibrant calculation and test for successfull exit
    returncode, stdout = run_vibrant(vibrant_binary_path, "input.txt")
    assert returncode == 0
    
    # file names
    reference = "output/resonance_raman_unpol.txt"
    test = "resonance_raman_unpol.txt"
    
    # parse filenames
    x_ref, y_ref = parse_spectrum_from_file(reference)
    x_test, y_test = parse_spectrum_from_file(test, hdr=1)

    # compare test against reference
    assert np.allclose(x_ref, x_test, atol=1e-8)
    scale = np.max(np.abs(y_ref))
    assert np.allclose(y_ref, y_test, rtol=1e-5, atol=max(1e-12, 1e-6 * scale))

    
