# HW4 Main File for LEO for EO
#
# Alec Lewandowski
# Last Edited 3/14/2026

import numpy as np
import matplotlib.pyplot as plt

def runP1():
    epoch = "2024-04-08T18:37:00"

    ## Part (a)
    hrTemp = 18 + 37/60
    JD = 2460408.5 + hrTemp/24
    print(JD)

    n = JD - 2451545.0
    L = 280.460 + 0.9856474*n # Mean longitude of sun
    g = 357.528 + 0.9856003*n # Mean anomaly

    GMSTdeg = (280.46061837 + 360.98564736629*n) % 360

    L = L % 360
    g = g % 360

    eclLon = L + 1.915*np.sin(np.deg2rad(g)) + 0.020*np.sin(np.deg2rad(2*g))
    eclLat = 0
    oblEcl = 23.439 - 0.0000004*n

    rAscen = np.rad2deg(np.arctan(np.cos(np.deg2rad(oblEcl)) * np.tan(np.deg2rad(eclLon))))
    delt = np.rad2deg(np.arcsin(np.sin(np.deg2rad(oblEcl)) * np.sin(np.deg2rad(eclLon))))

    geoLat = delt
    geoLon = rAscen - GMSTdeg

    print()
    print("Problem 1.a.i")
    print(f"   Ecliptic Latitude: {eclLat:.4f}°")
    print(f"   Ecliptic Longitude: {eclLon:.4f}°")
    print("Problem 1.a.ii")
    print(f"   Ecliptic Right Ascension: {rAscen:.4f}°")
    print(f"   Solar Declination: {delt:.4f}°")
    print("Problem 1.a.iii")
    print(f"   Geographical Latitude: {geoLat:.4f}°")
    print(f"   Geographical Longitude: {geoLon:.4f}°")

    ## Part (c)
    rAscenMoon = (15 + 11.36/4)
    decMoon = 7.37
    obsLon = -97.74
    obsLat = 30.27

    moonLon = rAscenMoon - GMSTdeg
    moonLat = decMoon
    print()
    print("Problem 2.c")
    print(f"   Geographical Latitude: {moonLat:.4f}°")
    print(f"   Geographical Longitude: {moonLon:.4f}°")


def runP2():
    ae = 6378.1363
    mu = 3.986004415e+05
    cryo = np.array([ae + 717, 0.0003565, 92.0194])
    topex = np.array([ae + 1336, 0.0007665, 66.0425])
    grace = np.array([ae + 490, 0.0005855, 88.9691])
    ERS = np.array([ae + 782, 0.0031394, 98.7932])
    pel = np.array([ae + 475, 0.0004512, 97.4358])
    star = np.array([ae + 812, 0.0205681, 49.8242])
    j2 = 0.001082

    abar = np.array([cryo[0], topex[0], grace[0], ERS[0], pel[0], star[0]])
    Ts = (2*np.pi)*(np.sqrt(abar**3 /mu))
    nbar = np.sqrt(mu / abar**3)
    ebar = np.array([cryo[1], topex[1], grace[1], ERS[1], pel[1], star[1]])
    ibar = np.deg2rad(np.array([cryo[2], topex[2], grace[2], ERS[2], pel[2], star[2]]))
    odotbar = (-3/2)*(nbar)*(ae/abar)**2 *j2*(1/(1-ebar**2)**2)*np.cos(ibar)
    wdotbar = (-3/4)*nbar*(ae/abar)**2 *j2*(1/(1-ebar**2)**2)*(1 - 5*np.cos(ibar)**2)
    mdotbar = nbar*(1 - (3/4)*(ae/abar)**2 *j2*(1/(1-ebar**2)**(3/2))*(1 - 3*np.cos(ibar**2)))
    udotbar = wdotbar + mdotbar
    thetadotg = 7.292115e-05
    omegaSun = 0.9856476*(np.pi/180)/86400
    Ta = (2*np.pi) / mdotbar
    Td = (2*np.pi) / udotbar
    Dn = (2*np.pi) / abs(odotbar - thetadotg)
    Cs = (2*np.pi) / abs(odotbar - (omegaSun))

    print()
    print("           Cryosat-2, TOPEX, GRACE, ERS, Pelican, Starlette")
    print(f"Values in order for Keplerian period: {Ts[0]:.4f}, {Ts[1]:.4f}, {Ts[2]:.4f}, {Ts[3]:.4f}, {Ts[4]:.4f}, and {Ts[5]:.4f}")
    print(f"Values in order for Anomalistic period: {Ta[0]:.4f}, {Ta[1]:.4f}, {Ta[2]:.4f}, {Ta[3]:.4f}, {Ta[4]:.4f}, and {Ta[5]:.4f}")
    print(f"Values in order for Draconitic period: {Td[0]:.4f}, {Td[1]:.4f}, {Td[2]:.4f}, {Td[3]:.4f}, {Td[4]:.4f}, and {Td[5]:.4f}")
    print(f"Values in order for nodal day: {Dn[0]:.4f}, {Dn[1]:.4f}, {Dn[2]:.4f}, {Dn[3]:.4f}, {Dn[4]:.4f}, and {Dn[5]:.4f}")
    print(f"Values in order for sun cycle: {Cs[0]/86400:.4f}, {Cs[1]/86400:.4f}, {Cs[2]/86400:.4f}, {Cs[3]/86400:.4f}, {Cs[4]/86400:.4f}, and {Cs[5]/86400:.4f}")
    print()

    oblSun = np.radians(23.44)

    tC = np.linspace(0, 4*Cs[0], 500)
    tT = np.linspace(0, 4*Cs[1], 500)
    tG = np.linspace(0, 4*Cs[2], 500)
    tS = np.linspace(0, 4*Cs[5], 500)

    omegaC = odotbar[0] * tC
    omegaT = odotbar[1] * tT
    omegaG = odotbar[2] * tG
    omegaS = odotbar[5] * tS

    LC = np.deg2rad(280.460 + 0.9856476*tC/86400)
    LT = np.deg2rad(280.460 + 0.9856476*tT/86400)
    LG = np.deg2rad(280.460 + 0.9856476*tG/86400)
    LS = np.deg2rad(280.460 + 0.9856476*tS/86400)

    dotProdC = 0.5*np.sin(ibar[0])*(1 - np.cos(oblSun))*np.sin(omegaC + LC) + 0.5*np.sin(ibar[0])*(1 + np.cos(oblSun))*np.sin(omegaC - LC) + np.cos(ibar[0])*np.sin(oblSun)*np.sin(LC)
    dotProdT = 0.5*np.sin(ibar[1])*(1 - np.cos(oblSun))*np.sin(omegaT + LT) + 0.5*np.sin(ibar[1])*(1 + np.cos(oblSun))*np.sin(omegaT - LT) + np.cos(ibar[1])*np.sin(oblSun)*np.sin(LT)
    dotProdG = 0.5*np.sin(ibar[2])*(1 - np.cos(oblSun))*np.sin(omegaG + LG) + 0.5*np.sin(ibar[2])*(1 + np.cos(oblSun))*np.sin(omegaG - LG) + np.cos(ibar[2])*np.sin(oblSun)*np.sin(LG)
    dotProdS = 0.5*np.sin(ibar[5])*(1 - np.cos(oblSun))*np.sin(omegaS + LS) + 0.5*np.sin(ibar[5])*(1 + np.cos(oblSun))*np.sin(omegaS - LS) + np.cos(ibar[5])*np.sin(oblSun)*np.sin(LS)

    betaC = 90 - np.degrees(np.arccos(dotProdC))
    betaT = 90 - np.degrees(np.arccos(dotProdT))
    betaG = 90 - np.degrees(np.arccos(dotProdG))
    betaS = 90 - np.degrees(np.arccos(dotProdS))

    plt.plot(tC/86400, betaC, label="Cryo-Sat2")
    plt.xlabel('Days')
    plt.ylabel("β' (degrees)")
    plt.title('Beta Prime')
    plt.legend()
    plt.grid(True)
    # plt.show()

    plt.plot(tT/86400, betaT, label="TOPEX")
    plt.xlabel('Days')
    plt.ylabel("β' (degrees)")
    plt.title('Beta Prime')
    plt.legend()
    plt.grid(True)
    # plt.show()

    plt.plot(tG/86400, betaG, label="GRACE")
    plt.xlabel('Days')
    plt.ylabel("β' (degrees)")
    plt.title('Beta Prime')
    plt.legend()
    plt.grid(True)
    # plt.show()

    plt.plot(tS/86400, betaS, label="Starlett")
    plt.xlabel('Days')
    plt.ylabel("β' (degrees)")
    plt.title('Beta Prime')
    plt.legend()
    plt.grid(True)
    # plt.show()

    crossings = []
    for i in range(len(betaC)-1):
        if betaC[i] * betaC[i+1] < 0:  # sign change
            # linear interpolation for more accurate crossing time
            t_cross = tC[i] - betaC[i] * (tC[i+1]-tC[i])/(betaC[i+1]-betaC[i])
            crossings.append(t_cross)
    timesC = np.array(t_cross)

    crossings = []
    for i in range(len(betaT)-1):
        if betaT[i] * betaT[i+1] < 0:  # sign change
            # linear interpolation for more accurate crossing time
            t_cross = tT[i] - betaT[i] * (tT[i+1]-tT[i])/(betaT[i+1]-betaT[i])
            crossings.append(t_cross)
    timesT = np.array(t_cross)

    crossings = []
    for i in range(len(betaG)-1):
        if betaG[i] * betaG[i+1] < 0:  # sign change
            # linear interpolation for more accurate crossing time
            t_cross = tG[i] - betaG[i] * (tG[i+1]-tG[i])/(betaG[i+1]-betaG[i])
            crossings.append(t_cross)
    timesG = np.array(t_cross)

    crossings = []
    for i in range(len(betaS)-1):
        if betaS[i] * betaS[i+1] < 0:  # sign change
            # linear interpolation for more accurate crossing time
            t_cross = tS[i] - betaS[i] * (tS[i+1]-tS[i])/(betaS[i+1]-betaS[i])
            crossings.append(t_cross)
    timesS = np.array(t_cross)

    print(f"Cryosat consecutive crossings: {timesC/(4*86400):.4f}")
    print(f"TOPEX consecutive crossings: {timesT/(4*86400):.4f}")
    print(f"GRACE consecutive crossings: {timesG/(4*86400):.4f}")
    print(f"Starlette consecutive crossings: {timesS/(4*86400):.4f}")

def runP3():
    ae = 3396.2
    mu = 4.2828e+04
    j2 = 0.00195545

    rps = np.array([100, 200, 400, 600, 1000, 2000, 2000]) + ae
    ras = np.array([100, 100, 100, 100, 100, 100, 2000]) + ae
    es = np.array([0.0, 0.008, 0.023, 0.037, 0.065, 0.128, 0.0])
    a = (rps + ras) / 2
    p = a * (1 - es**2)
    n = np.sqrt(mu / a**3)

    i = np.linspace(0, np.pi, 10000)

    omegadot1 = (-3/2)*((n[0]*j2*ae**2)/(p[0]**2))*np.cos(i)*(180/np.pi)*86400
    omegadot2 = (-3/2)*((n[1]*j2*ae**2)/(p[1]**2))*np.cos(i)*(180/np.pi)*86400
    omegadot3 = (-3/2)*((n[2]*j2*ae**2)/(p[2]**2))*np.cos(i)*(180/np.pi)*86400
    omegadot4 = (-3/2)*((n[3]*j2*ae**2)/(p[3]**2))*np.cos(i)*(180/np.pi)*86400
    omegadot5 = (-3/2)*((n[4]*j2*ae**2)/(p[4]**2))*np.cos(i)*(180/np.pi)*86400
    omegadot6 = (-3/2)*((n[5]*j2*ae**2)/(p[5]**2))*np.cos(i)*(180/np.pi)*86400
    omegadot7 = (-3/2)*((n[6]*j2*ae**2)/(p[6]**2))*np.cos(i)*(180/np.pi)*86400


    rps = np.array([100, 500, 1000, 2000, 3000, 4000, 4000]) + ae
    ras = np.array([100, 100, 100, 100, 100, 100, 2000]) + ae
    es = np.array([0.0, 0.030, 0.065, 0.128, 0.183, 0.231, 0.0])
    a = (rps + ras) / 2
    p = a * (1 - es**2)
    n = np.sqrt(mu / a**3)

    aopdot1 = (-3/2)*((n[0]*j2*ae**2)/(p[0]**2))*((5/2)*(np.sin(i))**2 - 2)*(180/np.pi)*86400
    aopdot2 = (-3/2)*((n[1]*j2*ae**2)/(p[1]**2))*((5/2)*(np.sin(i))**2 - 2)*(180/np.pi)*86400
    aopdot3 = (-3/2)*((n[2]*j2*ae**2)/(p[2]**2))*((5/2)*(np.sin(i))**2 - 2)*(180/np.pi)*86400
    aopdot4 = (-3/2)*((n[3]*j2*ae**2)/(p[3]**2))*((5/2)*(np.sin(i))**2 - 2)*(180/np.pi)*86400
    aopdot5 = (-3/2)*((n[4]*j2*ae**2)/(p[4]**2))*((5/2)*(np.sin(i))**2 - 2)*(180/np.pi)*86400
    aopdot6 = (-3/2)*((n[5]*j2*ae**2)/(p[5]**2))*((5/2)*(np.sin(i))**2 - 2)*(180/np.pi)*86400
    aopdot7 = (-3/2)*((n[6]*j2*ae**2)/(p[6]**2))*((5/2)*(np.sin(i))**2 - 2)*(180/np.pi)*86400

    plt.plot(np.degrees(i), omegadot1, label="100x100")
    plt.plot(np.degrees(i), omegadot2, label="100x200")
    plt.plot(np.degrees(i), omegadot3, label="100x400")
    plt.plot(np.degrees(i), omegadot4, label="100x600")
    plt.plot(np.degrees(i), omegadot5, label="100x1000")
    plt.plot(np.degrees(i), omegadot6, label="100x2000")
    plt.plot(np.degrees(i), omegadot7, label="100x2000")
    plt.xlabel('Inclination (°)')
    plt.ylabel("RAAN Dot (°/day)")
    plt.title('Figure 9.4 Recreation for Mars')
    plt.legend()
    plt.grid(True)
    # plt.show()


    plt.plot(np.degrees(i), aopdot1, label="100x100")
    plt.plot(np.degrees(i), aopdot2, label="100x500")
    plt.plot(np.degrees(i), aopdot3, label="100x1000")
    plt.plot(np.degrees(i), aopdot4, label="100x2000")
    plt.plot(np.degrees(i), aopdot5, label="100x3000")
    plt.plot(np.degrees(i), aopdot6, label="100x4000")
    plt.plot(np.degrees(i), aopdot7, label="4000x4000")
    plt.xlabel('Inclination (°)')
    plt.ylabel("AOP Dot (°/day)")
    plt.title('Figure 9.6 Recreation for Mars')
    plt.legend()
    plt.grid(True)
    # plt.show()
    
    

if __name__ == "__main__":
    runP1()
    runP2()
    runP3()