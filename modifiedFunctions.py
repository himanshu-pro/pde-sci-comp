import math


def modifiedH(xj, s, Zeta, Kr, Gr, Gm, Re, Du, Sc, Pr, Psi, R, J, S1, S2, Ha, D, Ec, Sr, Sh, Alpha, Rd, **kw):

    return [s[1],

            s[2],

            s[3],

            (1 / (1 + R))*(-R * (S1 * (s[2] + 2 * s[4]) + J * (s[0] * xj[5] + xj[0] * s[5] - s[1] * xj[4] - s[4] * xj[1]) * math.cos(Psi)) -
            Re * (s[1] * xj[2] + s[2] * xj[1] - s[0] * xj[3] - s[3] * xj[0]) * math.cos(Psi) +
            ((1 + R) * D + (Ha ** 2)) * s[2] -
            (Ec * Gr / Zeta) * (s[7] + (Zeta ** 2) * s[9])*math.sin(Alpha) -
            (Sh * Gm / Zeta) * (s[11] + (Zeta ** 2) * s[13])*math.sin(Alpha)-
            2*Ec*Gr*s[8]*math.cos(Alpha)-
            2*Sh*Gm*s[12]*math.cos(Alpha)),

            s[5],

            S1 * (s[2] + 2 * s[4]) +
            J * (s[0] * xj[5] + s[5] * xj[0] - s[1] * xj[4] - xj[1] * s[4]) * math.cos(Psi),

            s[7],


            (-1 / (1 - Du * Sc * Sr + Rd)) * (2 * s[8] * (1 - Du * Sc * Sr) +
            Re * Pr * (8 * s[1] * xj[1] + 2 * ((1 + R) * D + (Ha ** 2)) * s[0] * xj[0] - s[0] * xj[7] - s[7] * xj[0] + (2 * S2 / Pr) * s[4] * xj[4]) * math.cos(Psi) +
            Kr * Du * s[10] +
            Du * Sc * Re * (s[0] * xj[11] + xj[0] * s[11]) * math.cos(Psi)),

            s[9],

            -(1 / (1 - Du * Sc * Sr + Rd)) * (Re * Pr * ((R / 2) * (2 * s[2] * xj[2] + 8 * s[4] * xj[4] + 4 * s[2] * xj[4] + 4 * s[4] * xj[2]) + 2 * xj[2] * s[2] + ((1 + R) * D + (Ha ** 2)) * 2 * s[1] * xj[1] + 2 * s[1] * xj[8] + 2 * xj[1] * s[8] - s[0] * xj[9] - xj[0] * s[9] + (S2 / Pr) * 2 * s[5] * xj[5]) * math.cos(Psi) +
            Kr * Du * s[12] +
            Du * Sc * Re * (s[0] * xj[13] + xj[0] * s[13] - 2 * xj[1] * s[12] - 2 * s[1] * xj[12]) * math.cos(Psi)),

            s[11],

            -2 * s[12] +
            (1 / (1 - Du * Sc * Sr + Rd)) * (Sc * Sr * Re * Pr * (8 * s[1] * xj[1] + ((1 + R) * D + (Ha ** 2)) * 2 * s[0] * xj[0] - s[0] * xj[7] - s[7] * xj[0] + (2 * S2 / Pr) * s[4] * xj[4]) * math.cos(Psi) +
            Kr * s[10] * (1 + Rd) +
            Sc * Re * (s[0] * xj[11] + xj[0] * s[11]) * (1 + Rd) * math.cos(Psi) - Sc * Sr * Rd * s[8]),

            s[13],

            (1 / (1 - Du * Sc * Sr + Rd)) * (Sc * Sr * Re * Pr * ((R / 2) * (2 * s[2] * xj[2] + 8 * s[4] * xj[4] + 4 * s[2] * xj[4] + 4 * s[4] * xj[2]) + 2 * xj[2] * s[2] + ((1 + R) * D + (Ha ** 2)) * 2 * s[1] * xj[1] + 2 * s[1] * xj[8] + 2 * xj[1] * s[8] - s[0] * xj[9] - xj[0] * s[9] + (S2 / Pr) * 2 * s[5] * xj[5]) * math.cos(Psi) +
            Kr * s[12] * (1 + Rd)+
            Sc * Re * (1 + Rd) * (s[0] * xj[13] + xj[0] * s[13] - 2 * xj[1] * s[12] - 2 * s[1] * xj[12]) * math.cos(Psi))]


def modifiedP(xj, s, Zeta, Kr, Gr, Gm, Re, Du, Sc, Pr, Psi, R, J, S1, S2, Ha, D, Ec, Sr, Sh, Alpha, Rd, **kw):
    return [s[1],

            s[2],

            s[3],

            (1 / (1 + R)) * (-R * (S1 * (s[2] + 2 * s[4]) + J * (s[0] * xj[5] + xj[0] * s[5] - s[1] * xj[4] - s[4] * xj[1] + xj[1] * xj[4] - xj[0] * xj[5]) * math.cos(Psi)) -
            Re * (s[1] * xj[2] + s[2] * xj[1] - s[0] * xj[3] - s[3] * xj[0] + xj[0] * xj[3] - xj[1] * xj[2]) * math.cos(Psi) +
            ((1 + R) * D + (Ha ** 2)) * s[2] -
            (Ec * Gr / Zeta) * (s[7] + (Zeta ** 2) * s[9])*math.sin(Alpha) -
            (Sh * Gm / Zeta) * (s[11] + (Zeta ** 2) * s[13])*math.sin(Alpha) -
            2*Ec*Gr*s[8]*math.cos(Alpha) -
            2*Sh*Gm*s[12]*math.cos(Alpha)),

            s[5],

            S1 * (s[2] + 2 * s[4]) +
            J * (s[0] * xj[5] + s[5] * xj[0] - s[1] * xj[4] - xj[1] * s[4] + xj[1] * xj[4] - xj[0] * xj[5]) * math.cos(Psi),

            s[7],

            (1 / (1 - Du * Sc * Sr + Rd)) * (-2 * s[8] * (1 - Du * Sc * Sr) - Re * Pr * (8 * s[1] * xj[1] - 4*xj[1]*xj[1] + ((1 + R) * D + (Ha ** 2)) * (2*s[0] * xj[0] - xj[0]*xj[0]) - s[0] * xj[7] - s[7] * xj[0] + xj[7]*xj[0] + (S2 / Pr) * (2*s[4] * xj[4]-xj[4]*xj[4])) * math.cos(Psi) -
            Kr * Du * s[10] -
            Du * Sc * Re * (s[0] * xj[11] + xj[0] * s[11] - xj[0]*xj[11]) * math.cos(Psi)),

            s[9],

            (-1 / (1 - Du * Sc * Sr + Rd))*(Re * Pr * ((R / 2) * (2 * s[2] * xj[2] - xj[2]*xj[2] + 8 * s[4] * xj[4] - 4*xj[4]*xj[4] + 4 * s[2] * xj[4] + 4 * s[4] * xj[2] - 4*xj[4]*xj[2]) + 2 * xj[2] * s[2] - xj[2]*xj[2] + ((1 + R) * D + (Ha ** 2)) * (2 * s[1] * xj[1] - xj[1]*xj[1]) + 2 * s[1] * xj[8] + 2 * xj[1] * s[8] - s[0] * xj[9] - xj[0] * s[9] + xj[0] * xj[9] - 2*xj[1] * xj[8] + (S2 / Pr) * (2 * s[5] * xj[5] - xj[5]*xj[5])) * math.cos(Psi) +
            Kr * Du * s[12] +
            Du * Sc * Re * (s[0] * xj[13] + xj[0] * s[13] - 2 * xj[1] * s[12] - 2 * s[1] * xj[12] + 2*xj[1] * xj[12] - xj[0] * xj[13]) * math.cos(Psi)),

            s[11],

            -2 * s[12] +
            (1 / (1 - Du * Sc * Sr + Rd)) * (Sc * Sr * Re * Pr * (8 * s[1] * xj[1] - 4*xj[1]*xj[1] + ((1 + R) * D + (Ha ** 2)) * (2*s[0] * xj[0] - xj[0]*xj[0]) - s[0] * xj[7] - s[7] * xj[0] + xj[7]*xj[0] + (S2 / Pr) * (2*s[4] * xj[4]-xj[4]*xj[4])) * math.cos(Psi) +
            Kr * s[10] * (1 + Rd)+
            Sc * Re * (1 + Rd) * (s[0] * xj[11] + xj[0] * s[11] - xj[0]*xj[11]) * math.cos(Psi) - Sc * Sr * Rd * s[8]),

            s[13],

            (1 / (1 - Du * Sc * Sr + Rd)) * (Sc * Sr * Re * Pr * ((R / 2) * (2 * s[2] * xj[2] - xj[2]*xj[2] + 8 * s[4] * xj[4] - 4*xj[4]*xj[4] + 4 * s[2] * xj[4] + 4 * s[4] * xj[2] - 4*xj[4]*xj[2]) + 2 * xj[2] * s[2]- xj[2]*xj[2] + ((1 + R) * D + (Ha ** 2)) * (2 * s[1] * xj[1] - xj[1]*xj[1]) + 2 * s[1] * xj[8] + 2 * xj[1] * s[8] - s[0] * xj[9] - xj[0] * s[9]+ xj[0] * xj[9] - 2*xj[1] * xj[8] + (S2 / Pr) * (2 * s[5] * xj[5] - xj[5]*xj[5])) * math.cos(Psi) +
            Kr * s[12] * (1 + Rd) +
            Sc * Re * (1 + Rd) * (s[0] * xj[13] + xj[0] * s[13] - 2 * xj[1] * s[12] - 2 * s[1] * xj[12] + 2*xj[1] * xj[12] - xj[0] * xj[13]) * math.cos(Psi))]