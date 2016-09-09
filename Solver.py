def transpose(m):
    return [list(x) for x in zip(*m)]


def rk4(func, old, initial, **kw):
    cols = len(old[0])
    rows = (len(old)+1)/2

    h = 1.0/(rows - 1)

    res = [[0 for i in xrange(cols)] for j in xrange(rows)]
    res[0] = initial

    for j in xrange(1, rows):
        K1 = func(old[2*(j-1)], res[j-1], **kw)
        K1 = [h*v for v in K1]

        K2 = func(old[2*j-1], [res[j-1][col] + K1[col]/2.0 for col in xrange(cols)], **kw)
        K2 = [h * v for v in K2]

        K3 = func(old[2*j-1], [res[j-1][col] + K2[col]/2.0 for col in xrange(cols)], **kw)
        K3 = [h * v for v in K3]

        K4 = func(old[2*j], [res[j-1][col] + K3[col] for col in xrange(cols)], **kw)
        K4 = [h * v for v in K4]

        res[j] = [res[j-1][col] + (K1[col]+ K4[col])/6.0 + (K2[col] + K3[col])/3.0 for col in xrange(cols)]

    return res


def gauss_jordan(m, eps=1.0 / (10 ** 11)):
    (h, w) = (len(m), len(m[0]))

    for y in range(0, h):
        maxRow = y
        for y2 in range(y + 1, h):
            if abs(m[y2][y]) > abs(m[maxRow][y]):
                maxRow = y2

        (m[y], m[maxRow]) = (m[maxRow], m[y])

        if abs(m[y][y]) <= eps:
            return False

        for y2 in range(y + 1, h):
            c = m[y2][y] / m[y][y] * 1.0
            for x in range(y, w):
                m[y2][x] -= m[y][x] * c

    for y in range(h - 1, 0 - 1, -1):
        c = m[y][y]
        for y2 in range(0, y):
            for x in range(w - 1, y - 1, -1):
                m[y2][x] -= m[y][x] * m[y2][y] / c*1.0

        m[y][y] /= c*1.0
        for x in range(h, w):
            m[y][x] /= c*1.0

    return True


def interpolate(old):

    oldExp = old + old
    oldExp.pop()
    jLen = len(old)
    iLen = len(old[0])

    for j in xrange(2 * jLen - 1):
        if j % 2 == 0:
            oldExp[j] = old[j / 2]
        else:
            oldExp[j] = [(old[(j - 1) / 2][i] + old[(j + 1) / 2][i]) / 2.0 for i in xrange(iLen)]

    return oldExp


def linSolve(old, initialConds, boundary, fP, fH, **kw):

    condLen = len(initialConds)
    jLen = len(old)
    iLen = len(old[0])

    z = [old for k in xrange(condLen)]
    old = interpolate(old)

    for k, cond in enumerate(initialConds):
        if k == 0:
            z[k] = rk4(fP, old, cond, **kw)
        else:
            z[k] = rk4(fH, old, cond, **kw)

    system = transpose([[zk[-1][i] for i in boundary[0]] for k, zk in enumerate(z) if k > 0] + [[boundary[1][key]-z[0][-1][i] for key, i in enumerate(boundary[0])]])
    gauss_jordan(system)
    c = [1] + [sysRow[-1] for sysRow in system]

    return [[sum([z[ck][j][i] * cv for ck, cv in enumerate(c)]) for i in xrange(iLen)] for j in xrange(jLen)]


def nonLinSolve(Ec, Sh, A, Zeta, maxIter=10, **kw):
    eps = 1.0 / (10 ** 7)

    iLen = 14
    jLen = 1001

    old = [[0.1 for i in xrange(iLen)]for j in xrange(jLen)]

    init = [[1-A, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]]

    bound = [[0, 1, 4, 6, 8, 10, 12],
             [1.0, 0, 0, 1.0/Ec, 0, 1.0/Sh, 0]]

    error = 1.0
    it = 0

    while error > eps and it < maxIter:
        new = linSolve(old, init, bound, Ec=Ec, Sh=Sh, Zeta=Zeta, **kw)

        error = sum([sum([(new[j][i] - old[j][i]) ** 2 for i in xrange(iLen)]) for j in xrange(jLen)])
        old = new
        it += 1

    if error > eps:
        raise ValueError('Convergence not achieved in allowed no. of iterations')
    else:
        return [[Zeta*old[j][1] for j in xrange(jLen)],
                [old[j][0] for j in xrange(jLen)],
                [Zeta*old[j][4] for j in xrange(jLen)],
                [Ec*(old[j][6] + (Zeta**2)*old[j][8]) for j in xrange(jLen)],
                [Sh*(old[j][10] + (Zeta**2)*old[j][12]) for j in xrange(jLen)]]
