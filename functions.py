import numpy as np

def step_ieuler(f, tk, yk, h, args=None):
	"""	Compute solution value of initial value ODE problem after one step of the Improved Euler method.

		Parameters
		----------
		f : callable
			Derivative function.
		tk : float
			Initial value of independent variable.
		yk : float
			Initial value of solution.
		h : float
			Step size.
		args : iterable
			Optional parameters to pass into derivative function.

		Returns
		-------
		y : float
			Solution value after one step.
	"""
	
	if args is None:
		args = []
	
	# calculate the f0 and f1 derivatives using the derivative function provided
	f0 = f(tk, yk, *args)
	f1 = f(tk+h, yk+h*f0, *args)

	# calculate and return the y value after one step of the improved euler method
	y = yk + h*(f0/2 + f1/2)

	return y


def step_rk4(f, tk, yk, h, args=None):
	"""	Compute solution value of initial value ODE problem after one step of the classic RK4 method.

		Parameters
		----------
		f : callable
			Derivative function.
		tk : float
			Initial value of independent variable.
		yk : float
			Initial value of solution.
		h : float
			Step size.
		args : iterable
			Optional parameters to pass into derivative function.

		Returns
		-------
		y : float
			Solution value after one step.
	"""

	if args is None:
		args = []
	
	# calculate the f0, f1, f2, and f3 derivatives using the derivative function provided
	f0 = f(tk, yk, *args)
	f1 = f(tk+h/2, yk+(h*f0)/2, *args)
	f2 = f(tk+h/2, yk+(h*f1)/2, *args)
	f3 = f(tk+h, yk+h*f2, *args)

	# calculate and return the y value after one step of the classic RK4 method
	y = yk + h*((f0 + 2*f1 + 2*f2 + f3)/6)

	return y


def solve_explicit_rk(f, t0, t1, y0, h, method='rk4', args=None):
	"""	Compute solution of initial value ODE problem using explicit RK method.

		Parameters
		----------
		f : callable
			Derivative function.
		t0 : float
			Initial value of independent variable.
		t1 : float
			Final value of independent variable.
		y0 : float
			Initial value of solution.
		h : float
			Step size.
		method : str
			String specifying RK method, either 'rk4' or 'ieuler'. Default is 'rk4'.
		args : iterable
			Optional parameters to pass into derivative function.

		Returns
		-------
		t : array-like
			Independent variable at solution.
		y : array-like
			Solution.

		Notes
		-----
		Assumes that order of inputs to f is f(t, y, *args).
	"""
	if args is None:
		args = []

	# find number of steps to compute
	steps = (t1 - t0) / h

	# initial values
	t = [t0]
	y = [y0]

	# iterate through steps
	for i in range(int(steps)):

		# solving based on method
		if method == 'rk4':
			# perform step of classic rk4 to find next t and y values
			y.append(step_rk4(f,t[i],y[i],h,args))
			t.append(t[i]+ h)
			
		# if not rk4, we assume the method to be improved euler
		else:
			# perform step of improved euler to find next t and y values
			y[i+1] = step_ieuler(f,t[i],y[i],h,args)
			t[i+1] = t[i]+ h

	return t, y

def dndt_quota(t, n, r, k, f0):
	"""
	Derivative function for the original quota management system.
	"""

	return r*n*(1-n/k) - f0


def dndt_kaitiakitanga(t, n, r, k, f0, fr):
	"""
	Derivative function for the kaitiakitanga management system.
	"""

	return r*n*(1-n/k) - min(f0, fr*n)


def dndt_rahui(t, n, r, k, f0, x):
	"""
	Derivative function for the rahui management system.
	"""

	# assume start in period of fishing, then switch to no fishing etc
	if t % (2*x) < x:
		return r*n*(1-n/k) - f0
	else:
		return r*n*(1-n/k)
