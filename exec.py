import csv

import pygal

from Solver import nonLinSolve

from originalFunctions import originalP, originalH
from modifiedFunctions import modifiedP, modifiedH

defaults = {'Zeta': 0.22, 'Kr': 2.0, 'Gr': 4.0, 'Gm': 4.0, 'Re': 2.0, 'Du': 2.0, 'Sc': 0.22, 'Pr': 1.0, 'Psi': 0.2,
            'A': 0.2, 'R': 0.2, 'J': 0.2, 'S1': 2.0, 'S2': 10.0, 'Ha': 1.0, 'D': 2.0, 'Ec': 1.0, 'Sr': 1.5, 'Sh': 0.75,
            'Rd': 0.4, 'Alpha': 1.0}


# defaults = {'Zeta': 0.22, 'Kr': 2.0, 'Gr': 4.0, 'Gm': 4.0, 'Re': 2.0, 'Du': 0.2, 'Sc': 0.22, 'Pr': 0.71, 'Psi': 0.8,
#             'A': 0.2, 'R': 2.0, 'J': 0.2, 'S1': 2.0, 'S2': 1.0, 'Ha': 2.0, 'D': 2.0, 'Ec': 1.0, 'Sr': 0.2, 'Sh': 0.75,
#             'Rd': 0.4, 'Alpha': 1.0}

for k, v in {'Alpha': [0.2, 0.4, 0.6, 1.0]}.iteritems():
    res = []
    params = defaults.copy()
    for val in v:
        params[k] = val
        res.append(nonLinSolve(fP=modifiedP, fH=modifiedH, **params))

    # for i in xrange(len(res[0])):
    #     line = pygal.Line()
    #     for j in xrange(len(res)):
    #         line.add(k + '=' + str(v[j]), res[j][i])
    #     line.render_to_file('result_' + k + '_' + str(i) + '.svg')

    with open(k + str(v) + '.csv', 'wb') as fileDp:
        writer = csv.writer(fileDp)
        for i in xrange(len(res[0])):
            writer.writerows([res[j][i] for j in xrange(len(res))])
            writer.writerow([])
