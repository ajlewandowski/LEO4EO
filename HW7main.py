# HW7 Main File for LEO for EO
#
# Alec Lewandowski
# Last Edited 4/24/2026
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def getTombstone(alt, e, i, offset):
    ae = 6378
    abar = alt + ae
    ebar = e
    ibar = np.deg2rad(i)
    mu = 3.986004415e+05
    j2 = 0.001082
    thetadotg = 7.292115e-05

    GMST0 = (18.697375 / 24) * 2*np.pi

    nbar = np.sqrt(mu / abar**3)
    odg = (-3/2)*(nbar)*(ae/abar)**2 *j2*(1/(1-ebar**2)**2)*np.cos(ibar)
    wdg = (-3/4)*nbar*(ae/abar)**2 *j2*(1/(1-ebar**2)**2)*(1 - 5*np.cos(ibar)**2)
    mdg = nbar*(1 - (3/4)*(ae/abar)**2 *j2*(1/(1-ebar**2)**(3/2))*(1 - 3*np.cos(ibar)**2))
    udg = mdg + wdg

    towerLat = np.deg2rad(30 + 17/60 + 10.27/3600)
    towerLon = np.deg2rad(-(97 + 44/60 + 21.79/3600))
    print("########## LEO4EO HW7 P1 ##########")
    print(f"UT Tower Coordinates ({np.rad2deg(towerLat):.6f}, {np.rad2deg(towerLon):.6f})")
    print()
    print("Epoch: January 1st, 2000 12:00 UTC")
    print("Chosen Satellite: GRACE-FO-1")

    tVec = np.linspace(0,86400, 86401)

    x = ae * np.cos(towerLat) * np.cos(towerLon)
    y = ae * np.cos(towerLat) * np.sin(towerLon)
    z = ae * np.sin(towerLat)
    rTowerECEF = np.array([x, y, z])

    sma = np.zeros(len(tVec))
    ecc = np.zeros(len(tVec))
    inc = np.zeros(len(tVec))
    raan = np.zeros(len(tVec))
    aop = np.zeros(len(tVec))
    anom = np.zeros(len(tVec))
    arglat = np.zeros(len(tVec))
    satLats = np.zeros(len(tVec))
    satLons = np.zeros(len(tVec))
    vis = np.zeros(len(tVec))

    aop0 = np.deg2rad(0.0)
    raan0 = np.deg2rad(0.0)
    u0 = np.deg2rad(0.0 + offset)

    for k in range(len(tVec)):
        sma[k] = abar
        ecc[k] = ebar
        inc[k] = ibar
        aop[k] = aop0 + wdg*tVec[k]
        raan[k] = raan0 + odg*tVec[k]
        arglat[k] = u0 + udg*tVec[k]
        anom[k] = arglat[k] - aop[k]
        GMST = GMST0 + thetadotg*tVec[k]

        oesNow = np.array([abar, ebar, ibar, raan[k], aop[k], anom[k]])

        lat = np.arcsin(np.sin(ibar) * np.sin(arglat[k]))
        Lam = np.arcsin((np.sin(arglat[k])/np.cos(lat))*np.cos(ibar))
        lon = raan[k] + Lam - GMST

        satLats[k] = np.rad2deg(lat)
        satLons[k] = np.rad2deg(lon)

        r, _ = oe2rv(oesNow)

        eci2ecef = np.array([[np.cos(-GMST), np.sin(-GMST), 0], [-np.sin(-GMST), np.cos(-GMST), 0], [0, 0, 1]]).T

        rECEF = eci2ecef @ r

        rangeVec = rECEF - rTowerECEF
        lonNow = towerLon
        latNow = np.pi/2 - towerLat
        r3lon = np.array([[np.cos(lonNow), np.sin(lonNow), 0], [-np.sin(lonNow), np.cos(lonNow), 0], [0, 0, 1]])
        r2lat = np.array([[np.cos(latNow), 0, -np.sin(latNow)],[0, 1, 0],[np.sin(latNow), 0, np.cos(latNow)]])
        rangeSEZ = r2lat @ r3lon @ rangeVec

        if rangeSEZ[2] >= 0.0:
             vis[k] = 1.0
        else:
             vis[k] = 0.0

    return vis, satLats, satLons

def plotGroundTrack(lats1, lons1, lats2, lons2, lats3, lons3):
        fig = plt.figure(figsize=(14, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        # Add map features
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
        
        # Plot ground track
        ax.plot(lons1, lats1, 'r-', linewidth=2, transform=ccrs.Geodetic(), label='GRACE')
        ax.plot(lons1[0], lats1[0], 'ro', markersize=15, transform=ccrs.Geodetic(), label='Start')

        ax.plot(lons2, lats2, 'g-', linewidth=2, transform=ccrs.Geodetic(), label='GRACE 120')
        ax.plot(lons2[0], lats2[0], 'go', markersize=15, transform=ccrs.Geodetic(), label='Start')

        ax.plot(lons3, lats3, 'b-', linewidth=2, transform=ccrs.Geodetic(), label='GRACE 240')
        ax.plot(lons3[0], lats3[0], 'bo', markersize=15, transform=ccrs.Geodetic(), label='Start')
        
        # Set extent to show full track
        margin = 10
        ax.set_global()
        
        plt.legend(loc='upper right')
        plt.tight_layout()
        
        return fig

def oe2rv(oe):
    a, e, i, RAAN, aop, nu = oe

    p = a * (1 - e**2)
    rPQW = (p / (1 + e * np.cos(nu))) * np.array([np.cos(nu), np.sin(nu), 0])
    vPQW = np.sqrt(3.986004415e+05 / p) * np.array([-np.sin(nu), e + np.cos(nu), 0])

    # Rotation to ECI
    R = pqw2eci(RAAN, aop, i)
    r = R @ rPQW
    v = R @ vPQW

    return r, v

def pqw2eci(RAAN, aop, i):
    r3RAAN = np.array([[np.cos(RAAN), -np.sin(RAAN), 0], [np.sin(RAAN), np.cos(RAAN), 0], [0, 0, 1]])
    r1I = np.array([[1, 0, 0], [0, np.cos(i), -np.sin(i)], [0, np.sin(i), np.cos(i)]])
    r3aop = np.array([[np.cos(aop), -np.sin(aop), 0], [np.sin(aop), np.cos(aop), 0], [0, 0, 1]])

    return r3RAAN @ r1I @ r3aop

def getElongation(FOV, abar):
     a = abar
     ae = 6378
     alt = a - ae

     maxEta = np.arcsin(ae / a)
     eta = np.linspace(-maxEta, maxEta, 1000)

     rho = a*np.cos(eta) - np.sqrt(ae**2 - ((np.sin(eta))**2 * a**2))
     rhoToe = a*np.cos(eta + FOV/2) - np.sqrt(ae**2 - ((np.sin(eta + FOV/2))**2 * a**2))
     rhoHeel = a*np.cos(eta - FOV/2) - np.sqrt(ae**2 - ((np.sin(eta - FOV/2))**2 * a**2))

     lamHeel = np.arccos(1 + (alt - rhoHeel*np.cos(eta - FOV/2))/(ae))
     lamToe = np.arccos(1 + (alt - rhoToe*np.cos(eta + FOV/2))/(ae))

     LF = abs(ae*(lamToe - lamHeel))
     WF = abs(ae*np.arcsin(rho*np.sin(FOV)/ ae))
     A = (np.pi)/(4) * LF * WF

     return LF, WF, eta, A

if __name__ == "__main__":
    visGRACE, sLat1, sLon1 = getTombstone(490, 0.00005855, 88.9691, 0.0)
    print()
    print()
    visGRACEoff120, sLat2, sLon2 = getTombstone(490, 0.00005855, 88.9691, 120.0)
    print()
    print()
    visGRACEoff240, sLat3, sLon3 = getTombstone(490, 0.00005855, 88.9691, 240.0)

    tVec = np.linspace(0,86400, 86401)/60.0

    plt.plot(tVec, visGRACE)
    plt.xlabel('Time past 2000-01-01T12:00:00 UTC (min)')
    plt.ylabel("Visibility (1 visibile, 0 not)")
    plt.title('GRACE-FO-1 Tombstone Plot for Jan 1st, 2000 12hr UTC')
    plt.grid(True)
    plt.show()
    print()
    print()
    plt.plot(tVec, visGRACE, label="0.0° offset")
    plt.plot(tVec, visGRACEoff120,label="120.0° offset")
    plt.plot(tVec, visGRACEoff240, label="240.0° offset")
    plt.xlabel('Time past 2000-01-01T12:00:00 UTC (min)')
    plt.ylabel("Visibility (1 visibile, 0 not)")
    plt.title('GRACE-FO-1 Tombstone Plot for Jan 1st, 2000 12hr UTC')
    plt.legend()
    plt.grid(True)
    plt.show()
    print()
    print()
    plt.plot(tVec, visGRACE + visGRACEoff120 + visGRACEoff240)
    plt.xlabel('Time past 2000-01-01T12:00:00 UTC (min)')
    plt.ylabel("# of Sats Visible")
    plt.title('GRACE-FO-1 and Offsets Tombstone Plot for Jan 1st, 2000 12hr UTC')
    plt.grid(True)
    plt.show()
    print()
    print()

    fig2 = plotGroundTrack(sLat1, sLon1, sLat2, sLon2, sLat3, sLon3)
    plt.show()

    LF, WF, eta, A = getElongation(np.deg2rad(4), 490+6378)
    plt.plot(np.rad2deg(eta), LF)
    plt.xlabel('Off-Nadir Angle (deg)')
    plt.ylabel("Elongation (km)")
    plt.title(f'Elongation of a Conical FOV - Maximum Eta {np.rad2deg(eta[-1]):.3f}')
    plt.grid(True)
    plt.show()
    print()
    print()
    plt.plot(np.rad2deg(eta), WF)
    plt.xlabel('Off-Nadir Angle (deg)')
    plt.ylabel("Width (km)")
    plt.title(f'Width of a Conical FOV - Maximum Eta {np.rad2deg(eta[-1]):.3f}')
    plt.grid(True)
    plt.show()
    print()
    print()
    plt.plot(np.rad2deg(eta), A)
    plt.xlabel('Off-Nadir Angle (deg)')
    plt.ylabel("Area (km^2)")
    plt.title(f'Width of a Conical FOV - Maximum Eta {np.rad2deg(eta[-1]):.3f}')
    plt.grid(True)
    plt.show()
    print()
    print()