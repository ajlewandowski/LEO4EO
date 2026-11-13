# HW3 Main File for LEO for EO
#
# Alec Lewandowski
# Last Edited 2/21/2026

import numpy as np
import matplotlib.pyplot as plt
from twoBodyClass import twoBodyClass

def runP1():
    ## HW 3 Problem 1
    print()
    print()
    print("="*70)
    print("Homework 3 Problem 1")
    print("="*70)
    epochGiven = "2023-07-31T23:59:42"
    hrTemp = 42/3600 + 59/60 + 23
    JD0 = 2460156.500000
    JD = JD0 + hrTemp/24
    print("GMST Calculations for Part (a)")
    print(f"    JD0: {JD0:.8f} and JD: {JD:.8f}")
    Dtt = JD - 2451545.0
    Dut = JD0 - 2451545.0
    print(f"    Dtt: {Dtt:.8f} and Dut: {Dut:.8f}")
    T = Dtt / 36525
    print(f"    T: {T:.8f}")
    GMSThrs = (6.697375 + 0.065707485828*Dut + 1.0027379*hrTemp + 0.0854103*T + 0.0000258*T**2) % 24
    print(f"    GMST in hrs: {GMSThrs:.8f}")
    GMSTdeg = GMSThrs * 360 / 24
    print(f"    GMST in deg from USNO site: {GMSTdeg:.8f}")

    print("GAST Calculation for Part (b)")
    Omega = 125.04 - 0.052954*Dtt
    L = 280.47 + 0.98565*Dtt
    eps = 23.4393 - 0.0000004*Dtt
    dPsi = -0.000319*np.sin(Omega) - 0.000024*np.sin(2*L)
    GASThrs = GMSThrs + dPsi*np.cos(eps)
    print(f"    GAST in hours: {GASThrs:.8f}")
    GASTdeg = GASThrs * 360 / 24
    print(f"    GAST in degrees: {GASTdeg:.8f}")

    print("GMST and GAST Calculation for Part (c)")
    # Input time is going to be a JD in UTC
    inputTime = 2461089.163889 # Feb 17th, 2026 at 9:56 AM Central Time
    GMST, GAST = getGMSTandGAST(inputTime)

    print(f"    Results for Feb. 17th, 2026 at 9:56:00AM CT (JD = {inputTime:.8f})")
    print(f"        GMST in hours: {GMST[0]:.8f}")
    print(f"        GMST in degrees: {GMST[1]:.8f}")
    print(f"        GAST in hours: {GAST[0]:.8f}")
    print(f"        GAST in degrees: {GAST[1]:.8f}")


def runP2():
    ## HW 3 Problem 2
    print()
    print()
    print("="*70)
    print("Homework 3 Problem 2")
    print("="*70)
    JDepoch = 0.246015749979166666666666666667E+07
    stateECI = np.array([6093861.926471250131726, -1840393.570987288141623, 2590090.032104292418808,-2709.870119294875167, 975.605814756081486, 7047.604436212530345])
    stateECEF = np.array([5259791.807165756821632, 3575614.384020817000419, 2603908.595081204082817, -2211.985280491541744, -1889.418848876843867, 7041.442349260451920])
    r = np.array([stateECEF[0], stateECEF[1], stateECEF[2]])
    v = np.array([stateECEF[3], stateECEF[4], stateECEF[5]])

    rECI = np.array([stateECI[0], stateECI[1], stateECI[2]])
    vECI = np.array([stateECI[3], stateECI[4], stateECI[5]])

    GMST, _ = getGMSTandGAST(JDepoch)

    rECEF, vECEF = eci2ecef(GMST[1], stateECI)

    latApprox, lonApprox = getLatLon(rECEF)
    latECEF, lonECEF = getLatLon(r)

    errLat = np.rad2deg(abs(latApprox - latECEF))
    errLon = np.rad2deg(abs(lonApprox - lonECEF))

    print("Latitude and Longitude Calculations for Part (a)")
    print(f"    Latitude approximation: {np.rad2deg(latApprox):.8f}°")
    print(f"    Longitude approximation: {np.rad2deg(lonApprox):.8f}°")

    print("Latitude and Longitude Error for Part (b)")
    print(f"    Latitude conversion from given ECEF epoch: {np.rad2deg(latECEF):.8f}°")
    print(f"    Longitude conversion from given ECEF epoch: {np.rad2deg(lonECEF):.8f}°")
    print(f"    Latitude Error: {errLat:.8f}°")
    print(f"    Longitude Error: {errLon:.8f}°")

    latTriApprox, lonTriApprox = getLatLonViaTriangle(rECI, vECI, GMST[1])
    latTriECEF, lonTriECEF = getLatLonViaTriangle(rECI, vECI, GMST[1])

    errTriLat = np.rad2deg(abs(latTriApprox - latTriECEF))
    errTriLon = np.rad2deg(abs(lonTriApprox - lonTriECEF))

    print("Latitude and Longitude Error for Part (c)")
    print(f"    Latitude conversion from given ECEF epoch: {np.rad2deg(latTriECEF):.8f}°")
    print(f"    Longitude conversion from given ECEF epoch: {np.rad2deg(lonTriECEF):.8f}°")
    print(f"    Latitude conversion approximation from (a): {np.rad2deg(latTriApprox):.8f}°")
    print(f"    Longitude conversion approximation from (a): {np.rad2deg(lonTriApprox):.8f}°")
    print(f"    Latitude Error: {errTriLat:.8f}°")
    print(f"    Longitude Error: {errTriLon:.8f}°")

def getGMSTandGAST(inputTime):
    # Using the modular to back out how many hours are in this given JD
    dayFracRem = (inputTime % 1.0) + 0.5
    hrsRem = dayFracRem * 24

    JD = inputTime
    JD0 = inputTime - dayFracRem

    Dtt = JD - 2451545.0
    Dut = JD0 - 2451545.0

    T = Dtt / 36525

    GMSThrs = (6.697375 + 0.065707485828*Dut + 1.0027379*hrsRem + 0.0854103*T + 0.0000258*T**2) % 24
    GMSTdeg = GMSThrs * 360 / 24

    Omega = 125.04 - 0.052954*Dtt
    L = 280.47 + 0.98565*Dtt
    eps = 23.4393 - 0.0000004*Dtt
    dPsi = -0.000319*np.sin(Omega) - 0.000024*np.sin(2*L)

    GASThrs = GMSThrs + dPsi*np.cos(eps)
    GASTdeg = GASThrs * 360 / 24

    GMST = np.array([GMSThrs, GMSTdeg])
    GAST = np.array([GASThrs, GASTdeg])

    return GMST, GAST

def eci2ecef(GMSTdeg, stateECI):
    x, y, z, xdot, ydot, zdot = stateECI
    GMST = np.deg2rad(GMSTdeg)
    
    R3 = np.array([[np.cos(-GMST), -np.sin(-GMST), 0], [np.sin(-GMST), np.cos(-GMST), 0], [0, 0, 1]])

    rECI = np.array([x, y, z])
    vECI = np.array([xdot, ydot, zdot])
    
    rECEF = R3 @ rECI
    vECEF = R3 @ vECI

    return rECEF, vECEF

def getLatLon(r):

    lon = np.atan2(r[1], r[0])

    p = np.sqrt(r[0]**2 + r[1]**2)
    lat = np.atan2(r[2], p)

    return lat, lon

def getLatLonViaTriangle(r, v, GMSTdeg):
    GMST = np.deg2rad(GMSTdeg)

    JDepoch = 0
    eMu = 3.986004415E+14
    eRad = 6378136.3
    eAngVel = 7.292115E-05
    j2 = 0.001082
    orb = twoBodyClass(eMu, eRad, eAngVel, j2, JDepoch)

    a, e, i, aop, RAAN, M, u = orb.rv2oe(r, v) # in radians

    phi = np.asin(np.sin(i) * np.sin(u))
    lam = np.arctan2(np.cos(i) * np.sin(u), np.cos(phi))

    lon = RAAN + lam - GMST
    lat = phi

    return lat, lon


if __name__ == "__main__":
    runP1()
    runP2()