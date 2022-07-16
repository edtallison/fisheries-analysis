# Tests for the functions in functions.py.

from functions import solve_explicit_rk, step_ieuler, step_rk4
import numpy as np

# error tolerance
tol = 1.e-10

def dydt1(t, y):
    return t - y

def dydt2(t, y):
    return np.cos(t)

def test_step_ieuler():
    """
    Test implementation of the improved euler method for numerical solutions of ODEs.
    """

    # correct solution for single step
    correct = 3

    # call step_ieuler function
    y = step_ieuler(dydt1, 0, 1, 2)

    # compare calculated to correct
    assert abs(y - correct) < tol

def test_step_rk4():
    """
    Test implementation of the classic RK4 method for numerical solutions of ODEs.
    """
    
    # correct solution for single step
    correct = 5/3

    # call step_rk4 function
    y = step_rk4(dydt1, 0, 1, 2)

    # compare calculated to correct
    assert abs(y - correct) < tol

def test_solve_explicit_rk():
    """
    Test implementation of solving of explicit rk methods for numerical solutions of ODEs.
    This uses the mean absolute error, compared to an error tolerance.
    """

    # call the test_solve_explicit_rk function
    t, y = solve_explicit_rk(dydt2, 0, 10, 0, 1, 'rk4')

    # calculate sum of differences between function and exact for computing mean absolute error
    sum = 0
    for i in range(10):
        sum += abs(y[i] - np.sin(i))

    # mean absolute error
    mae = sum / 10
    
    # compare against error tolerance
    assert mae < 0.01
    