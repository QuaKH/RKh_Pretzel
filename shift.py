def orientation(l, m, n):
    """
    orientation takes in l, m, n as input variables
    and outputs the orientation of the pretzel diagram
    """
    if l % 2 == 0:
        return "RL"
    if l % 2 == 1 and m % 2 == 1:
        if n % 2 == 1:
            return "RR"
        else:
            return "LR"
    if l % 2 == 1 and m % 2 == 0 and n % 2 == 1:
        return "LL"
    # if (l % 2 == 0 and m % 2 == 1 and n % 2 == 1):
    #     return "RL"
    # missing a case for l odd, m even, n even


def shift_table(l, m, n):
    """
    shift_table takes in l, m, n as input variables
    and outputs the corresponding polynomial shifts
    """
    s_l = 0
    t_l = 0
    s_u = 0
    t_u = 0
    orient = orientation(l, m, n)

    if orient == "RR":
        s_l = - 2 * m - 2 * n - 1
        t_l = - m - n
        s_u = 4 * l - 2 * m - 2 * n - 1
        t_u = 2 * l - m - n + 1

    elif orient == "LL":
        s_l = - 3 * l - 2 * m + n - 1
        t_l = - l - m
        s_u = l - 2 * m + n - 1
        t_u = l - m + 1

    elif orient == "LR":
        s_l = - 3 * l + m - 2 * n - 1
        t_l = - l - n
        s_u = l + m - 2 * n - 1
        t_u = l - n + 1

    else:
        s_l = m + n - 1
        t_l = 0
        s_u = 4 * l + m + n - 3
        t_u = 2 * l

    return s_l, t_l, s_u, t_u
