# HW6 Main File for LEO for EO
#
# Alec Lewandowski
# Last Edited 4/12/2026

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import math

def runP1():

    ae = 6378
    abar = np.array([717, 1336, 490, 782, 475, 812]) + ae
    ebar = np.array([0.0003565, 0.0007665, 0.00005855, 0.0031394, 0.0004512, 0.0205681])
    ibar = np.deg2rad(np.array([92.0194, 66.0425, 88.9691, 98.7932, 97.4358, 49.8242]))
    mu = 3.986004415e+05
    j2 = 0.001082

    ## Part A - Westward Drift of the Geographic Longitude ##
    nbar = np.sqrt(mu / abar**3)
    odotbar = (-3/2)*(nbar)*(ae/abar)**2 *j2*(1/(1-ebar**2)**2)*np.cos(ibar)
    wdotbar = (-3/4)*nbar*(ae/abar)**2 *j2*(1/(1-ebar**2)**2)*(1 - 5*np.cos(ibar)**2)
    mdotbar = nbar*(1 - (3/4)*(ae/abar)**2 *j2*(1/(1-ebar**2)**(3/2))*(1 - 3*np.cos(ibar**2)))
    udotbar = wdotbar + mdotbar
    thetadotg = 7.292115e-05
    qs = (odotbar - thetadotg) / udotbar
    dLons = np.rad2deg(2*np.pi*qs)

    print()
    print("Westward shift of the ascending crossing")
    print(f"Cryosat-2: Δλ = {dLons[0]:.4f}")
    print(f"TOPEX: Δλ = {dLons[1]:.4f}")
    print(f"GRACE-FO-1: Δλ = {dLons[2]:.4f}")
    print(f"ERS-1: Δλ = {dLons[3]:.4f}")
    print(f"Pelican-2: Δλ = {dLons[4]:.4f}")
    print(f"Starlette: Δλ = {dLons[5]:.4f}")

    ## Part B - Calculate the {m, k} repeat conditions within 45 nodal days ##
    ## repeat condition is k * Dn = m Td
    
    Dn = abs((2*np.pi)/(odotbar - thetadotg))
    Td = abs(2*np.pi / udotbar)
    pairD = []
    pairR = []
    pairSat = []
    for j in range(6):
        ratio = Td[j] / Dn[j]
        tol = 0.05
        for k in range(1, 46):
            m = k / ratio
            mR = round(m)

            if mR == 0:
                continue

            if abs(m - mR) < tol:
                if math.gcd(int(k), int(mR)) == 1:
                    days = k
                    pairD.append(days)
                    pairR.append(mR)
                    pairSat.append(j)
    
    D = np.array(pairD)
    R = np.array(pairR)
    Sat = np.array(pairSat)

    print()
    print("Primal pair results")
    for k in range(len(D)):
        if Sat[k] == 0:
            print(f"Cryosat-2 Pair: {D[k]} days and {R[k]} revs")
        if Sat[k] == 1:
            print(f"TOPEX Pair: {D[k]} days and {R[k]} revs")   
        if Sat[k] == 2:
            print(f"GRACE-FO-1 Pair: {D[k]} days and {R[k]} revs")
        if Sat[k] == 3:
            print(f"ERS-1 Pair: {D[k]} days and {R[k]} revs")
        if Sat[k] == 4:
            print(f"Pelican-2 Pair: {D[k]} days and {R[k]} revs")
        if Sat[k] == 5:
            print(f"Starlette Pair: {D[k]} days and {R[k]} revs")

    ## Part C - Defining the criterion for misclosure ##
    ## misclosure for time -> dt = k*Dn - m*Td
    print()
    print("Misclosure Results")
    for k in range(len(D)):
        if Sat[k] == 0:
            dt = Dn[0]*D[k] - R[k]*Td[0]
            dl = thetadotg * dt
            print(f"Cryosat-2 misclosure: {dt:.4f} seconds and {np.rad2deg(dl):.4f}°")
        if Sat[k] == 1:
            dt = Dn[1]*D[k] - R[k]*Td[1]
            dl = thetadotg * dt
            print(f"TOPEX misclosure: {dt:.4f} seconds and {np.rad2deg(dl):.4f}°")  
        if Sat[k] == 2:
            dt = Dn[2]*D[k] - R[k]*Td[2]
            dl = thetadotg * dt
            print(f"GRACE-FO-1 misclosure: {dt:.4f} seconds and {np.rad2deg(dl):.4f}°")
        if Sat[k] == 3:
            dt = Dn[3]*D[k] - R[k]*Td[3]
            dl = thetadotg * dt
            print(f"ERS-1 misclosure: {dt:.4f} seconds and {np.rad2deg(dl):.4f}°")
        if Sat[k] == 4:
            dt = Dn[4]*D[k] - R[k]*Td[4]
            dl = thetadotg * dt
            print(f"Pelican-2 misclosure: {dt:.4f} seconds and {np.rad2deg(dl):.4f}°")
        if Sat[k] == 5:
            dt = Dn[5]*D[k] - R[k]*Td[5]
            dl = thetadotg * dt
            print(f"Starlett misclosure: {dt:.4f} seconds and {np.rad2deg(dl):.4f}°")

    
    ## Part D - misclosure less than 0.5 degrees ##
    aFracSpace = np.linspace(0.95, 1.05, 100000)
    print()
    print("SMA Misclosure Sweep Results")
    for i in range(6):
        a = abar[i]
        inc = ibar[i]
        ecc = ebar[i]
        k = D[Sat == i][0]
        m = R[Sat == i][0]

        aSweep = a * aFracSpace
        dl = compMisclosure(aSweep, inc, ecc, k, m)

        valid = aSweep[dl < (0.5*np.pi/180)]
        if len(valid) > 0:
            print(f"Sat {i}: a in [{valid[0]:.4f}, {valid[-1]:.4f}] km  "
                f"(width = {(valid[-1]-valid[0]):.4f} km)")
        else:
            print(f"Sat {i}: no valid range found — check k,m pair or widen sweep")


def compMisclosure(a, i, e, k, m):
    ae = 6378
    mu = 3.986004415e+05
    j2 = 0.001082
    n = np.sqrt(mu / a**3)
    odotbar = (-3/2)*(n)*(ae/a)**2 *j2*(1/(1-e**2)**2)*np.cos(i)
    wdotbar = (-3/4)*n*(ae/a)**2 *j2*(1/(1-e**2)**2)*(1 - 5*np.cos(i)**2)
    mdotbar = n*(1 - (3/4)*(ae/a)**2 *j2*(1/(1-e**2)**(3/2))*(1 - 3*np.cos(i**2)))
    udotbar = wdotbar + mdotbar
    thetadotg = 7.292115e-05

    Dn = abs(2*np.pi / (odotbar - thetadotg))
    Td = abs(2*np.pi / udotbar)

    misTime = k*Dn - m*Td
    dL = abs(thetadotg*misTime)
    return dL



if __name__ == "__main__":
    runP1()