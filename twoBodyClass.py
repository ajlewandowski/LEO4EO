# This file is the Keplerian orbit propagator class for the satellite/earth problem

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

class twoBodyClass:

    def __init__(self, mu, eRad, eAngVel, j2, epochStart):
        self.mu = mu
        self.eRad = eRad
        self.eAngVel = eAngVel
        self.j2 = j2
        
        self.epochStart = epochStart

    def oe2rv(self, oe):
        """ Converting classical orbital elements to position and velocity """
        
        a, e, i, aop, RAAN, M, u = oe

        # True anomaly
        nu = self.mean2nu(M, e)

        # Perifocal frame
        p = a * (1 - e**2)
        rPQW = (p / (1 + e * np.cos(nu))) * np.array([np.cos(nu), np.sin(nu), 0])
        vPQW = np.sqrt(self.mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0])

        # Rotation to ECI
        R = self.pqw2eci(RAAN, aop, i)
        r = R @ rPQW
        v = R @ vPQW

        return r, v
    
    def rv2oe(self, r, v):
        r_norm = np.linalg.norm(r)
        v_norm = np.linalg.norm(v)
        h = np.cross(r, v)
        h_norm = np.linalg.norm(h)
        mu = self.mu
        
        # Semi-major axis and eccentricity
        energy = v_norm**2 / 2 - mu / r_norm
        a = -mu / (2*energy)
        e_vec = (1/mu)*np.cross(v,h) - r/r_norm
        e = np.linalg.norm(e_vec)
        
        # Inclination
        i = np.arccos(h[2]/h_norm)
        
        # RAAN
        n = np.cross([0,0,1], h)
        n_norm = np.linalg.norm(n)
        RAAN = np.arccos(n[0]/n_norm)
        if n[1] < 0: RAAN = 2*np.pi - RAAN
        
        # Argument of periapsis
        aop = np.arccos(np.dot(n,e_vec)/(n_norm*e))
        if e_vec[2] < 0: aop = 2*np.pi - aop
        
        # True anomaly
        nu = np.arccos(np.dot(e_vec,r)/(e*r_norm))
        if np.dot(r,v) < 0: nu = 2*np.pi - nu
        
        # Mean anomaly
        E = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(nu/2))
        M = E - e*np.sin(E)
        M = M % (2*np.pi)

        u = aop + nu
        u = u % (2*np.pi)
        
        return [a, e, i, aop, RAAN, M, u]

    def pqw2eci(self, RAAN, aop, i):
        """ Rotation matrix from perifocal frame to ECI frame """

        r3RAAN = np.array([[np.cos(RAAN), -np.sin(RAAN), 0],[np.sin(RAAN), np.cos(RAAN), 0],[0, 0, 1]])
        r1I = np.array([[1, 0, 0],[0, np.cos(i), -np.sin(i)],[0, np.sin(i), np.cos(i)]])
        r3aop = np.array([[np.cos(aop), -np.sin(aop), 0],[np.sin(aop), np.cos(aop), 0],[0, 0, 1]])

        return r3RAAN @ r1I @ r3aop

    def mean2nu(self, M, e, tol=1e-10):
        """ Converting mean anomaly to true anomaly """

        E = M
        for k in range(50):
            f = E - e * np.sin(E) - M
            fp = 1 - e * np.cos(E)
            Enew = E - f / fp
            if abs(Enew - E) < tol:
                E = Enew
                break
            E = Enew

        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

        return nu

    def orbitEOMs(self, t, y):
        """ This function produces the numerical keplerian equations of motions """  
        # unpack state
        rx, ry, rz, vx, vy, vz = y
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])

        # norm of radius vector
        normR = np.linalg.norm(r)

        # two body acceleration
        ax, ay, az = -r * self.mu / normR**3     
        
        return np.array([vx, vy, vz, ax, ay, az])
    
    def propagate(self, y0, tF, dtSec):

        print(f"Propagating orbit")
        print(f"    Duration: {tF/60} minutes")
        print(f"    Time Step: {dtSec} seconds")
        nP = int(tF / dtSec) + 1
        tSpan = np.linspace(0, tF, nP)

        # Integrate
        sol = solve_ivp(self.orbitEOMs, [0, tF], y0, method='RK45', t_eval=tSpan, rtol=1e-12, atol=1e-12)

        return sol.t, sol.y
    
    def keplerProp(self, epochR, epochV, dt):
        print(f"Solving Kepler's problem to propagate orbit")
        print(f"    Duration: {dt} seconds")

        normR = np.linalg.norm(epochR)
        normV = np.linalg.norm(epochV)

        zeta = normV**2 / 2 - (self.mu / normR)
        alpha = -(normV**2 / self.mu) + (2 / normR)
        a = 1 / alpha
        chi = np.sqrt(self.mu) * dt * alpha

        diffTol = 1.0

        for _ in range (250):
            psi = alpha * (chi**2)
            c2, c3 = self.findc2c3(psi)
            r = (c2 * chi**2) + (np.dot(epochR, epochV))/np.sqrt(self.mu) *chi*(1 - psi*c3) + normR * (1 - psi*c2)
            chiN = chi + ((np.sqrt(self.mu)*dt - c3*chi**3) - ((np.dot(epochR, epochV)/np.sqrt(self.mu))*c2*chi**2) - (normR*chi*(1 - psi*c3))) / r
            diffTol = abs(chi - chiN)
            chi = chiN
            if diffTol < 1E-06:
                break

        f = 1 - (c2*((chi**2) / normR))
        g = dt - (c3*((chi**3)/np.sqrt(self.mu)))

        fdot = np.sqrt(self.mu)/(r*normR)*chi*(psi*c3 - 1)
        gdot = 1 - (chi**2 / r * c2)

        rF = f*epochR + g*epochV
        vF = fdot*epochR + gdot*epochV

        checkVal = f*gdot - fdot*g
        return rF, vF, checkVal
    
    def findc2c3(self, psi):
        
        c2 = (1 - np.cos(np.sqrt(psi))) / psi
        c3 = (np.sqrt(psi) - np.sin(np.sqrt(psi))) / (np.sqrt(psi**3))

        return c2, c3
    
    def rotMatP5(self, ang3, ang2):

        r3 = np.array([[np.cos(ang3), -np.sin(ang3), 0],[np.sin(ang3), np.cos(ang3), 0],[0, 0, 1]])
        r2 = np.array([[np.cos(ang2), 0, np.sin(ang2)],[0, 1, 0],[-np.sin(ang2), 0, np.cos(ang2)]])

        return r3, r2
    
    def plotOEs(self, t, a, aAvg, eltitle,glabel,ylabel):

        fig = plt.figure(figsize=(10, 6))
        plt.plot(t, a, 'b-', linewidth=2, label=glabel)
        plt.axhline(y=aAvg,color='r',linestyle='--',linewidth=2, label=f'Mean = {aAvg:.2f}')
        plt.title(eltitle)
        plt.xlabel('Time since Epoch (s)')
        plt.ylabel(ylabel)
        plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=2)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
    
        return fig
    
    def getFFT(self, et, a, eltitle):
        nP = len(et)
        dt = et[1] - et[0]
        
        # data_detrended = a - np.mean(a)
        
        yf = fft(a)
        xf = fftfreq(nP, dt)
        
        positive_freq_idx = xf > 0
        xf_pos = xf[positive_freq_idx]
        magnitude = np.abs(yf[positive_freq_idx])
        magmax = max(magnitude)
        magmin = min(magnitude)
        normFac = (magmax - magmin) / 2
        amps = magnitude / normFac

        idgraph = np.where(xf_pos > 0.001)[0]
        xftest = xf_pos[0:idgraph[0]]
        magnew = amps[0:idgraph[0]]

        # Find top N dominant frequencies
        top_5_idx = np.argsort(magnew)[::-1][:5]
        
        print(f"\n{'='*70}")
        print(f"Top {5} Dominant Frequencies for {eltitle}:")
        print(f"{'='*70}")
        print(f"{'Rank':<6} {'Frequency (Hz)':<20} {'Period (s)':<18} {'Amplitude':<15}")
        print(f"{'-'*70}")
        
        for k in range(5):
            idx = top_5_idx[k]
            freq = xftest[idx]
            period = 1/freq
            amp = magnew[idx]
            print(f"{k+1:<6} {freq:<20.10f} {period:<18.4f} {amp:<15.6e}")
        
        print()

        fig = plt.figure(figsize=(10, 6))
        plt.plot(xftest, magnew, 'r-', linewidth=2, label='Signal')
        plt.title(eltitle)
        plt.xlabel('Frequencies (Hz)')
        plt.ylabel(eltitle)
        plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=1)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        return fig
    
    def plotLongOEs(self, t, a, eltitle, glabel, ylabel):

        fig = plt.figure(figsize=(10, 6))
        plt.plot(t, a, 'b-', linewidth=2, label=glabel)
        plt.title(eltitle)
        plt.xlabel('Orbital Revolutions')
        plt.ylabel(ylabel)
        plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=2)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
    
        return fig
    
    def eci2rtn(self, r, RAAN, i, u):

        r3RAAN = np.array([[np.cos(RAAN), -np.sin(RAAN), 0],[np.sin(RAAN), np.cos(RAAN), 0],[0, 0, 1]])
        r1I = np.array([[1, 0, 0],[0, np.cos(i), -np.sin(i)],[0, np.sin(i), np.cos(i)]])
        r3u = np.array([[np.cos(u), -np.sin(u), 0],[np.sin(u), np.cos(u), 0],[0, 0, 1]])

        rNew = r3u @ r1I @ r3RAAN @ r
        
        return rNew
    
    def plotRTNDiff(self, t, rDiff, eltitle, ylabel):

        fig = plt.figure(figsize=(10, 6))
        plt.plot(t, rDiff[:,0], 'r-', linewidth=2, label='Radial')
        plt.plot(t, rDiff[:,1], 'b-', linewidth=2, label='Transverse')
        plt.plot(t, rDiff[:,2], 'g-', linewidth=2, label='Normal')
        plt.title(eltitle)
        plt.xlabel('Time since Epoch (s)')
        plt.ylabel(ylabel)
        plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=3)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
    
        return fig
    
    def plotaBar(self, rev, aData, aBar):

        fig = plt.figure(figsize=(10, 6))
        plt.plot(rev, aData, 'r-', linewidth=1, label='SMA')
        plt.axhline(y=aBar,color='b',linestyle='--',linewidth=2, label=f'Mean = {aBar:.2f}')
        plt.title('Oscullating Twice-per-Rev Amplitude')
        plt.xlabel('Orbit #')
        plt.ylabel('SMA Estimate (m)')
        plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=2)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
    
        return fig