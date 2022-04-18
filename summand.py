import numpy as np
import sympy as sp


def V_tilde(a):
    """
    V_tilde takes a numpy array 'a' and outputs
    """
    q, t = sp.symbols('q t')
    x = 0
    for i in range(len(a)):
        x += a[i] * (q ** (4 * i)) * t ** (2 * i)
    y = sp.Poly(x, t)

    return x, y, y.all_coeffs()

vec = [1, 2, 1]
print(V_tilde(vec))
