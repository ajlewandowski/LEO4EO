# HW5 Main File for LEO for EO
#
# Alec Lewandowski
# Last Edited 4/3/2026

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def runP1():
    a = np.linspace(6378, 50000, 100000)
    T = 2*np.pi*np.sqrt(a**3 /398600.4418)
    rho = np.arcsin(6378/a)
    tShadows = T * rho / np.pi
    f = rho / np.pi

    plt.plot(a, tShadows)
    plt.xlabel('Semi-Major Axis (km)')
    plt.ylabel("Time in Shadow (s)")
    plt.title('Shadow Time as a function of Semi-Major Axis')
    plt.legend()
    plt.grid(True)
    plt.show()

    plt.plot(a, f*100)
    plt.xlabel('Semi-Major Axis (km)')
    plt.ylabel("Shadow Fraction (%)")
    plt.title('Shadow Fraction as a function of Semi-Major Axis')
    plt.legend()
    plt.grid(True)
    plt.show()

def runP2():

    ae = 6378
    semis = np.array([717, 1336, 490, 782, 475, 812]) + ae
    es = np.array([0.0003565, 0.0007665, 0.00005855, 0.0031394, 0.0004512, 0.0205681])
    ibar = np.deg2rad(np.array([92.0194, 66.0425, 88.9691, 98.7932, 97.4358, 49.8242]))

    ubar = np.linspace(0, 360, 1000)

    lam1 = (12 + (ubar * np.cos(ibar[0]))/15) % 24
    lam2 = (12 + (ubar * np.cos(ibar[1]))/15) % 24
    lam3 = (12 + (ubar * np.cos(ibar[2]))/15) % 24
    lam4 = (12 + (ubar * np.cos(ibar[3]))/15) % 24
    lam5 = (12 + (ubar * np.cos(ibar[4]))/15) % 24
    lam6 = (12 + (ubar * np.cos(ibar[5]))/15) % 24

    plt.plot(ubar, lam4, label="Cryosat-2")
    plt.plot(ubar, lam2, label="TOPEX")
    plt.plot(ubar, lam5, label="GRACE-FO-1")
    plt.plot(ubar, lam1, label="ERS-1")
    plt.plot(ubar, lam3, label="Pelican")
    plt.plot(ubar, lam6, label="Starlette")
    plt.xlabel('Argument of Latitude (°)')
    plt.ylabel("Local Time at SSP (hours)")
    plt.title('SSP as a function of the orbit')
    plt.legend()
    plt.grid(True)
    plt.show()

def runP3():


    abar = 6378 + 1336
    ebar = 0.0007665
    ibar = np.deg2rad(66.0425)

    epochUTC = "1992-12-21T14:47:00"
    RAAN = 198.651
    betaP = 39.262

    dayFrac = 133 + 14/24 + 38/3600 + 53/86400
    RAANapprox = (116.56 - 2.08*(dayFrac)) % 360

    print()
    print(f"Approximated RAAN value: {RAANapprox:.3f}°")
    print(f"RAAN value from hmwk: {RAAN:.3f}°")
    print(f"Error in the approximated: {abs(RAAN - RAANapprox):.3f}°")
    print()

    RAAN = np.deg2rad(RAANapprox)
    betaP = np.deg2rad(betaP)

    eh = np.array([np.sin(ibar)*np.sin(RAAN), -np.sin(ibar)*np.cos(RAAN), np.cos(ibar)])
    n = np.array([np.cos(RAAN), np.sin(RAAN), 0.0])
    esun = np.cos(betaP)*n + np.sin(betaP)*eh
    A = 2.0 * 1361.0

    uVals = np.linspace(0, 2*np.pi, 10000)
    sunSat = np.zeros(len(uVals))
    pR = np.zeros(len(uVals))
    pT = np.zeros(len(uVals))
    pH = np.zeros(len(uVals))
    for k in range(len(uVals)):
        er = np.array([np.cos(RAAN)*np.cos(uVals[k]) - np.sin(RAAN)*np.sin(uVals[k])*np.cos(ibar),
                       np.sin(RAAN)*np.cos(uVals[k]) + np.cos(RAAN)*np.sin(uVals[k])*np.cos(ibar),
                       np.sin(ibar)*np.sin(uVals[k])])
        temp = np.clip(esun @ er, -1.0, 1.0)

        et = np.cross(eh, er)
        et /= np.linalg.norm(et)
        ssT = np.arccos(np.clip(np.dot(esun, et), -1.0, 1.0))
        sunSat[k] = np.arccos(temp)
        
        pR[k] = A * np.cos(sunSat[k])
        pT[k] = A * np.cos(ssT)
        pH[k] = A * np.cos(np.deg2rad(50.74))


    shadowSSA = []
    shadowU = []
    for j in range(len(uVals)):
        if 90.0 < np.rad2deg(sunSat[j]) <= 180:
            r = abar * (1 - ebar**2) / (1 + ebar*np.cos(uVals[j]))
            if r*np.sin(np.pi - sunSat[j]) < abar:
                shadowSSA.append(sunSat[j])
                shadowU.append(uVals[j])
                pR[j] = 0
                pT[j] = 0
                pH[j] = 0
    sSSA = np.array(shadowSSA)
    sU = np.array(shadowU)

    for k in range(len(uVals)):
        if pR[k] < 0:
            pR[k] = 0.0
        if pT[k] < 0:
            pT[k] = 0.0
        if pH[k] < 0:
            pH[k] = 0.0

    print(pH)

    plt.plot(np.rad2deg(uVals), np.rad2deg(sunSat), label="TOPEX")
    plt.plot(np.rad2deg(sU), np.rad2deg(sSSA), label="Shadow")
    plt.xlabel('Argument of Latitude (°)')
    plt.ylabel("Sun-Satellite Angle (°)")
    plt.title('SSA as a function of the orbit')
    plt.legend()
    plt.grid(True)
    plt.show()

    
    plt.plot(np.rad2deg(uVals), pR, label="Radial")
    plt.plot(np.rad2deg(uVals), pT, label="Transverse")
    plt.plot(np.rad2deg(uVals), pH, label="Normal")
    plt.xlabel('Argument of Latitude (°)')
    plt.ylabel("Available Power (W)")
    plt.title('Power as a function of the orbit')
    plt.legend()
    plt.grid(True)
    plt.show()

def runP4():
    alt = 350
    j2 = 0.001082
    inc = np.deg2rad(35.0)
    # epoch 1999-01-21T20:43:47
    # RO14 cross at 22:20
    geoLon = -5.157
    ebar = 0.0

    ae = 6378
    a = ae + alt
    n = np.sqrt(398600.4418 / a**3)

    thetadotg = 7.292115e-5
    odotbar = (-3/2)*(n)*(ae/a)**2 *j2*(1/(1-ebar**2)**2)*np.cos(inc)
    wdotbar = (-3/4)*n*(ae/a)**2 *j2*(1/(1-ebar**2)**2)*(1 - 5*np.cos(inc)**2)
    mdotbar = n*(1 - (3/4)*(ae/a)**2 *j2*(1/(1-ebar**2)**(3/2))*(1 - 3*np.cos(inc)**2))
    udotbar = wdotbar + mdotbar
    thetadotg = 7.292115e-05
    Td = (2*np.pi) / udotbar

    dLonN = np.rad2deg((odotbar - thetadotg) * Td)

    LMT_dot = (np.degrees(-1.5 * n * j2 * (ae/a)**2 * np.cos(inc)) * 86400 - 360/365.25) / 15
    epoch, LMT0 = datetime(1999,1,21,20,43,47), 20+43/60+47/3600 + (-5.157/15)
    target = 22 + 20/60
    
    dates = [epoch + timedelta(days=(target - LMT0 + 24*k) / LMT_dot)
            for k in range(-10, 50)
            if datetime(1999,1,1) <= epoch + timedelta(days=(target - LMT0 + 24*k) / LMT_dot) <= datetime(1999,12,31)]

    print(*dates, sep='\n')

        
if __name__ == "__main__":
    # runP1()
    # runP2()
    runP3()
    # runP4()