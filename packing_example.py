# Example/template for using packing of arguments in functions.

def solve_explicit_rk(f, t0, t1, y0, h, method='rk4', args=None):
    if args is None:
        args = []

    tk = 0.
    yk = y0
    tspan = t1 - t0
    # other commands

    # standardised ode derivative function call
    derivative = f(tk, yk, *args)


def dydt1(t, y):
    return t - y


def dydt2(t, y, a, b):
    return a * t - b * y


solve_explicit_rk(dydt1, 0., 1., 1., 0.1, 'rk4')
solve_explicit_rk(dydt2, 0., 1., 1., 0.1, 'rk4', [3., 4.])
