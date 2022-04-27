import numpy as np
import sympy as sp
import shift


def poly_helper(poly, q, t):
    """
    returns 3 lists; coefficients, power of q's, power of t's
    """
    qt_poly = sp.Poly(poly, q, t)
    terms = qt_poly.terms()
    q_power = []
    t_power = []
    coeffs = []
    for i in range(len(terms)):
        q_power.append(terms[i][0][0])
        t_power.append(terms[i][0][1])
        coeffs.append(terms[i][1])
    return q_power, t_power, coeffs


def V_tilde(a):
    """
    V_tilde takes a numpy array 'a' and outputs
    """
    q, t = sp.symbols('q t')
    x = 0
    for i in range(len(a)):
        x += a[i] * (q ** (4 * i)) * t ** (2 * i)
    ppp = poly_helper(x, q, t)

    return x, ppp


def VaE(a, E):
    q, t = sp.symbols('q t')
    aprime = []
    exceptional = 0
    for i in range(len(a)):
        aprime.append(a[i] - E[i])
        exceptional += E[i] * (1 + q ** 2) * (q ** (4 * i - 4)) * t ** (2 * i - 2)
    knights_move = sp.expand((1 + q ** 4 * t) * V_tilde(aprime)[0])
    output = sp.expand(knights_move * exceptional)
    ppp = poly_helper(output, q, t)
    return output, ppp


def Va(a):
    q, t = sp.symbols('q t')
    knights_move = sp.expand((1 + q ** 4 * t) * V_tilde(a)[0])
    output = sp.expand(knights_move)
    ppp = poly_helper(output, q, t)
    return output, ppp


def special_sequence(sequence_type, k):
    output = []
    if k <= 0:
        return output

    if sequence_type == "a":
        for i in range(k):
            output.append(i + 1)
            if len(output) >= k:
                return output
            output.append(i)
            if len(output) >= k:
                return output

    if sequence_type == "b":
        for i in range(k):
            output.append(i + 1)
            if len(output) >= k:
                return output
            if i <= 1:
                output.append(0)
            else:
                output.append(i - 1)
            if len(output) >= k:
                return output

    if sequence_type == "c":
        for i in range(k):
            output.append(i + 1)
            if len(output) >= k:
                return output
            output.append(i + 1)
            if len(output) >= k:
                return output

    return 0


def lower_summand(l, m, n):
    if (m != l) and (l % 2 == 1):
        sequence = special_sequence("a", l - 1)
        sequence += [int((l - 1) / 2)] * 2
        sequence += special_sequence("a", l - 1)[::-1]
        Lower_sum = Va(sequence)
        return sequence, Lower_sum

    if (m != l) and (l % 2 == 0):
        sequence = [1] + special_sequence("c", l - 4)
        sequence += [int((l - 2) / 2)] * 3 + [int(l / 2)]
        sequence += special_sequence("b", l)[::-1]
        E = [1] + [0] * (len(sequence) - 1)
        Lower_sum = VaE(sequence, E)
        return sequence, Lower_sum

    if (m == l) and (l % 2 == 1):
        sequence = special_sequence("a", l - 1) + [int((l - 1) / 2)] * 2
        sequence += special_sequence("a", l - 1)[::-1]
        sequence += [0] * (n - l) + [1]
        E = [0] * (len(sequence) - 1) + [1]
        Lower_sum = VaE(sequence, E)
        return sequence, Lower_sum

    if (m == l) and (l % 2 == 0) and (n % 2 == 1):
        sequence = [1] + special_sequence("c", l - 4)
        sequence += [int((l - 2) / 2)] * 3 + [int(l / 2)]
        sequence += special_sequence("b", l)[::-1] + [0, 1] * (int((n - l - 1) / 2))
        E = [1] + [0] * (len(sequence) - 1)
        Lower_sum = VaE(sequence, E)
        return sequence, Lower_sum

    if (m == l) and (l % 2 == 0) and (n % 2 == 0) and (n != l):
        sequence = [1] + special_sequence("c", l - 4)
        sequence += [int((l - 2) / 2)] * 3 + [int(l / 2)]
        sequence += special_sequence("b", l)[::-1] + [0, 1] * (int((n - l) / 2))
        E = [1] + [0] * (len(sequence) - 2) + [1]
        Lower_sum = VaE(sequence, E)
        return sequence, Lower_sum

    if (m == l) and (l % 2 == 0) and (n == l):
        sequence = [1] + special_sequence("c", l - 4)
        sequence += [int((l - 2) / 2)] * 3 + [int(l / 2)]
        temp = special_sequence("b", l)[0]
        temp[0] += 1
        sequence += temp[::-1]
        E = [1] + [0] * (len(sequence) - 2) + [2]
        Lower_sum = VaE(sequence, E)
        return sequence, Lower_sum


def upper_summand(l, m, n):
    g = m - l
    h = n - m

    if (l % 2 == 1) and (m % 2 == 1) and (n % 2 == 1):
        sequence = special_sequence("a", g)
        sequence += [int(g / 2)] * h
        sequence += special_sequence("a", g)[::-1]
        E = [0] * (len(sequence) - 1) + [1]
        upper_sum = VaE(sequence, E)
        return sequence, upper_sum

    if (l % 2 == 1) and (m % 2 == 1) and (n % 2 == 0):
        sequence = special_sequence("a", g)
        sequence += [int(g / 2)] * (h - 1)
        sequence += special_sequence("c", g)[::-1]
        E = [0] * (g + h - 1) + [1] + [0] * (g - 1)
        upper_sum = VaE(sequence, E)
        return sequence, upper_sum

    if (l % 2 == 1) and (m % 2 == 0) and (n % 2 == 1):
        sequence = special_sequence("a", g - 1)
        sequence += [int((g + 1) / 2), int((g - 1) / 2)] * (int((h + 1) / 2))
        sequence += special_sequence("c", g - 1)[::-1]
        E = [0] * (g - 1) + [1] + [0] * h + [0] * (g - 1)
        upper_sum = VaE(sequence, E)
        return sequence, upper_sum

    if (l % 2 == 1) and (m % 2 == 0) and (n % 2 == 0):
        sequence = special_sequence("a", g - 1)
        sequence += [int((g + 1) / 2), int((g - 1) / 2)] * (int(h / 2)) + [int((g + 1) / 2)]
        sequence += special_sequence("c", g - 1)[::-1]
        if h == 0:
            E = [0] * (g - 1) + [2] + [0] * (g - 1)
            upper_sum = VaE(sequence, E)
            return sequence, upper_sum
        else:
            E = [0] * (g - 1) + [1] + [0] * (h - 1) + [1] + [0] * (g - 2)
            upper_sum = VaE(sequence, E)
            return sequence, upper_sum

    if (l % 2 == 0) and (m % 2 == 1) and (n % 2 == 1):
        sequence = special_sequence("b", g + 1)
        sequence += [int((g + 1) / 2), int((g - 1) / 2)] * (int(h / 2))
        sequence += special_sequence("c", g - 1)[::-1]
        upper_sum = Va(sequence)
        return sequence, upper_sum

    if (l % 2 == 0) and (m % 2 == 1) and (n % 2 == 0):
        sequence = special_sequence("b", g + 1)
        sequence += [int((g + 1) / 2), int((g - 1) / 2)] * (int((h - 1) / 2))
        sequence += special_sequence("c", g)[::-1]
        E = [0] * (g + 1) + [0] * (h - 1) + [1] + [0] * (g - 1)
        upper_sum = VaE(sequence, E)
        return sequence, upper_sum

    if (l % 2 == 0) and (m % 2 == 0) and (n % 2 == 1):
        sequence = special_sequence("b", g)
        if g == 0:
            sequence += [1, 0] * (int((h + 1) / 2))
            sequence += special_sequence("c", 0)[::-1]
        else:
            sequence += [int((g + 2) / 2), int((g - 2) / 2)] * (int((h + 1) / 2))
            sequence += special_sequence("c", g - 1)[::-1]
        E = [0] * g + [1] + [0] * h + [0] * (g - 1)
        upper_sum = VaE(sequence, E)
        return sequence, upper_sum

    if (l % 2 == 0) and (m % 2 == 0) and (n % 2 == 0):
        sequence = special_sequence("b", g)
        if g == 0:
            sequence += [1, 0] * (int(h / 2)) + [1]
            sequence += special_sequence("a", 0)[::-1]
        else:
            sequence += [int((g + 2) / 2), int((g - 2) / 2)] * (int(h / 2)) + [int((g + 2) / 2)]
            sequence += special_sequence("a", g)[::-1]
        if h == 0:
            E = [0] * g + [2] + [0] * (g - 1) + [1]
            upper_sum = VaE(sequence, E)
            return sequence, upper_sum
        else:
            E = [0] * g + [1] + [0] * (h - 1) + [1] + [0] * (g - 1) + [1]
            upper_sum = VaE(sequence, E)
            return sequence, upper_sum


def Kh(l, m, n):
    if 2 <= l <= m <= n:
        q, t = sp.symbols('q t')
        s_l, t_l, s_u, t_u = shift.shift_table(l, m, n)
        part1 = (q ** int(s_l)) * (t ** int(t_l)) * lower_summand(l, m, n)[1][0]
        part2 = (q ** int(s_u)) * (t ** int(t_u)) * upper_summand(l, m, n)[1][0]
        homology_poly = sp.expand(part1 + part2)
        q_power, t_power, coeffs = poly_helper(homology_poly, q, t)

        return q_power, t_power, coeffs
    else:
        return -404, -404, -404


# vec = [1, 2, 1]
# # print(V_tilde(vec, q, )[0])
# print(VaE(vec, [0,1,0])[1])

# print(special_sequence("c", 10))
# print(lower_summand(2,3,4)[0])
print(lower_summand(2,3,4)[1])
print(Kh(2, 3, 4))
