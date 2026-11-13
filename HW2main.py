# HW2 Main File for LEO for EO
#
# Alec Lewandowski
# Last Edited 1/27/2026

import numpy as np
import matplotlib.pyplot as plt
from twoBodyClass import twoBodyClass

# # Loading in the data for the GRACE Follow-On Mission for the first time (stored as array file later for quicker access)
# file_path = r'C:\Users\alecl\OneDrive\Desktop\LEO for EO\HW2\hw2.data\1yr.005s.ECI.txt'
# dataMat = np.loadtxt(file_path, unpack=True)
# np.save('epochDataHW2P1.npy', dataMat)
# file_path = r'C:\Users\alecl\OneDrive\Desktop\LEO for EO\HW2\hw2.data\pred.1d.005s.ECI.txt'
# dataMat = np.loadtxt(file_path, unpack=True)
# np.save('predDataHW2P2.npy', dataMat)
# file_path = r'C:\Users\alecl\OneDrive\Desktop\LEO for EO\HW2\hw2.data\osc.1d.005s.ECI.txt'
# dataMat = np.loadtxt(file_path, unpack=True)
# np.save('oscDataHW2P2.npy', dataMat)

epochMat = np.load('epochDataHW2P1.npy')
et, eMJD, ex, ey, ez, exdot, eydot, ezdot = epochMat
oeMat = np.load('oeDataHW2P1.npy')
a, e, i, aop, RAAN, M, u = oeMat.T
predMat = np.load('predDataHW2P2.npy')
pt, pMJD, px, py, pz, pxdot, pydot, pzdot = predMat
oscMat = np.load('oscDataHW2P2.npy')
ot, oMJD, ox, oy, oz, oxdot, oydot, ozdot = oscMat

# Starting epoch is the GRACE Follow-on mission
epochStart = "2023-08-01T00:00:00" # in UTC time
eMu = 3.986004415E+14
eRad = 6378136.3
eAngVel = 7.292115E-05
j2 = 0.001082

# Initializing propagator class
prop = twoBodyClass(eMu, eRad, eAngVel, j2, epochStart)

# # One time run to store oe data, makes file run much faster calculating OEs once
# n = len(ex)
# oeData = np.zeros((n, 7))
# for k in range(n):
#     r = np.array([ex[k], ey[k], ez[k]])
#     v = np.array([exdot[k], eydot[k], ezdot[k]])
#     oes = prop.rv2oe(r, v)
#     oeData[k,:] = oes
# np.save('oeDataHW2P1.npy', oeData)

## HW 2 Problem 1
print()
print("="*70)
print("Homework 2 Problem 1")
print("="*70)
print("Initial Orbital Elements:")
print(f"  a = {a[0]:.7f} m")
print(f"  e = {e[0]:.8f}")
print(f"  i = {np.degrees(i[0]):.7f}°")
print(f"  ω = {np.degrees(aop[0]):.7f}°")
print(f"  Ω = {np.degrees(RAAN[0]):.7f}°")
print(f"  M = {np.degrees(M[0]):.7f}°")
print(f"  u = {np.degrees(u[0]):.7f}°")
print()

# First finding the indices of the time-series where the equator crossing happens
eqZDiffs = np.diff(np.sign(ez))
eqCrossInd = np.where(eqZDiffs > 0)[0] + 1
eqCrossInd = np.insert(eqCrossInd, 0, 0)
eqCrossInd = np.append(eqCrossInd, len(et))

# Now calculating the mean orbital elements between each crossing (per revolution)
aAvg = []
eAvg = []
iAvg = []
aopAvg = []
RAANAvg = []
MAvg = []
uAvg = []
revNum = []
tAvg = []
for k in range(len(eqCrossInd)-1):
    aAvg.append(np.mean(a[eqCrossInd[k]:(eqCrossInd[k+1])]))
    eAvg.append(np.mean(e[eqCrossInd[k]:(eqCrossInd[k+1])]))
    iAvg.append(np.mean(i[eqCrossInd[k]:(eqCrossInd[k+1])]))
    aopAvg.append(np.mean(aop[eqCrossInd[k]:(eqCrossInd[k+1])]))
    RAANAvg.append(np.mean(RAAN[eqCrossInd[k]:(eqCrossInd[k+1])]))
    MAvg.append(np.mean(M[eqCrossInd[k]:(eqCrossInd[k+1])]))
    uAvg.append(np.mean(u[eqCrossInd[k]:(eqCrossInd[k+1])]))
    revNum.append(k+1)
    tAvg.append(et[eqCrossInd[k]])

# Printing TLE data to check with posted average values
print()
print('Arithmetic average of orbital elements for first 5 orbits:')
print(f"{eqCrossInd[1]-eqCrossInd[0]:.0f}   {aAvg[0]:.4f}   {eAvg[0]:.8f}    {np.degrees(iAvg[0]):.7f}     {np.degrees(aopAvg[0]):.7f}     {np.degrees(RAANAvg[0]):.7f}     {np.degrees(MAvg[0]):.7f}    {np.degrees(uAvg[0]):.7f}")
print(f"{eqCrossInd[2]-eqCrossInd[1]:.0f}   {aAvg[1]:.4f}   {eAvg[1]:.8f}    {np.degrees(iAvg[1]):.7f}     {np.degrees(aopAvg[1]):.7f}     {np.degrees(RAANAvg[1]):.7f}     {np.degrees(MAvg[1]):.7f}    {np.degrees(uAvg[1]):.7f}")
print(f"{eqCrossInd[3]-eqCrossInd[2]:.0f}   {aAvg[2]:.4f}   {eAvg[2]:.8f}    {np.degrees(iAvg[2]):.7f}     {np.degrees(aopAvg[2]):.7f}     {np.degrees(RAANAvg[2]):.7f}     {np.degrees(MAvg[2]):.7f}    {np.degrees(uAvg[2]):.7f}")
print(f"{eqCrossInd[4]-eqCrossInd[3]:.0f}   {aAvg[3]:.4f}   {eAvg[3]:.8f}    {np.degrees(iAvg[3]):.7f}     {np.degrees(aopAvg[3]):.7f}     {np.degrees(RAANAvg[3]):.7f}     {np.degrees(MAvg[3]):.7f}    {np.degrees(uAvg[3]):.7f}")
print(f"{eqCrossInd[5]-eqCrossInd[4]:.0f}   {aAvg[4]:.4f}   {eAvg[4]:.8f}    {np.degrees(iAvg[4]):.7f}     {np.degrees(aopAvg[4]):.7f}     {np.degrees(RAANAvg[4]):.7f}     {np.degrees(MAvg[4]):.7f}    {np.degrees(uAvg[4]):.7f}")


# Plotting each of the orbital elements 
dayInd = np.where(et == 86400)[0]
fig1 = prop.plotOEs(et[0:dayInd[0]], a[0:dayInd[0]], np.mean(a[0:dayInd[0]]), eltitle='Semi-Major Axis (SMA) Evolution', glabel='SMA', ylabel='SMA (m)')
fig2 = prop.plotOEs(et[0:dayInd[0]], e[0:dayInd[0]], np.mean(e[0:dayInd[0]]), eltitle='Eccentricity (Ecc) Evolution', glabel='Ecc', ylabel='Eccentricity')
fig3 = prop.plotOEs(et[0:dayInd[0]], np.degrees(i[0:dayInd[0]]), np.degrees(np.mean(i[0:dayInd[0]])), eltitle='Inclination (Inc) Evolution', glabel='Inc', ylabel='Inc (°)')
fig4 = prop.plotOEs(et[0:dayInd[0]], np.degrees(aop[0:dayInd[0]]), np.degrees(np.mean(aop[0:dayInd[0]])), eltitle='Argument of Perigee (AOP) Evolution', glabel='AOP', ylabel='AOP (°)')
fig5 = prop.plotOEs(et[0:dayInd[0]], np.degrees(RAAN[0:dayInd[0]]), np.degrees(np.mean(RAAN[0:dayInd[0]])), eltitle='Right Acenscion of the Ascending Node (RAAN) Evolution', glabel='RAAN', ylabel='RAAN (°)')
fig6 = prop.plotOEs(et[0:dayInd[0]], np.degrees(M[0:dayInd[0]]), np.degrees(np.mean(M[0:dayInd[0]])), eltitle='Mean Anomaly (M) Evolution', glabel='M', ylabel='M (°)')
fig7 = prop.plotOEs(et[0:dayInd[0]], np.degrees(u[0:dayInd[0]]), np.degrees(np.mean(u[0:dayInd[0]])), eltitle='Argument of Latitude (u) Evolution', glabel='u', ylabel='u (°)')

# Using a fourier transform to get frequencies of orbital elements
fig24 = prop.getFFT(et[0:dayInd[0]], a[0:dayInd[0]], 'FFT of Semi-Major Axis')
fig25 = prop.getFFT(et[0:dayInd[0]], e[0:dayInd[0]], 'FFT of Eccentricity')
fig26 = prop.getFFT(et[0:dayInd[0]], np.degrees(i[0:dayInd[0]]), 'FFT of Inclination')
fig27 = prop.getFFT(et[0:dayInd[0]], np.degrees(aop[0:dayInd[0]]), 'FFT of Argument of Perigee')
fig28 = prop.getFFT(et[0:dayInd[0]], np.degrees(RAAN[0:dayInd[0]]), 'FFT of Right Ascension of the Ascending Node')
fig29 = prop.getFFT(et[0:dayInd[0]], np.degrees(M[0:dayInd[0]]), 'FFT of Mean Anomaly')
fig30 = prop.getFFT(et[0:dayInd[0]], np.degrees(u[0:dayInd[0]]), 'FFT of Argument of Latitude')

# Plotting over the entire year
fig8 = prop.plotLongOEs(revNum[1:-2], aAvg[1:-2], eltitle='Semi-Major Axis (SMA) Year Evolution', glabel='SMA', ylabel='SMA (m)')
fig9 = prop.plotLongOEs(revNum[1:-2], eAvg[1:-2], eltitle='Eccentricity (Ecc) Year Evolution', glabel='Ecc', ylabel='Eccentricity')
fig10 = prop.plotLongOEs(revNum[1:-2], np.degrees(iAvg[1:-2]), eltitle='Inclination (Inc) Year Evolution', glabel='Inc', ylabel='Inc (°)')
fig11 = prop.plotLongOEs(revNum[1:-2], np.degrees(aopAvg[1:-2]), eltitle='Argument of Perigee (AOP) Year Evolution', glabel='AOP', ylabel='AOP (°)')
fig12 = prop.plotLongOEs(revNum[1:-2], np.degrees(RAANAvg[1:-2]), eltitle='Right Acenscion of the Ascending Node (RAAN) Year Evolution', glabel='RAAN', ylabel='RAAN (°)')
fig13 = prop.plotLongOEs(revNum[1:-2], np.degrees(MAvg[1:-2]), eltitle='Mean Anomaly (M) Year Evolution', glabel='M', ylabel='M (°)')
fig14 = prop.plotLongOEs(revNum[1:-2], np.degrees(uAvg[1:-2]), eltitle='Argument of Latitude (u) Year Evolution', glabel='u', ylabel='u (°)')

## HW 2 Problem 2
print()
print("="*70)
print("Homework 2 Problem 2")
print("="*70)
diffMat = predMat - oscMat
dt, dMJD, dx, dy, dz, dxdot, dydot, dzdot = diffMat

# Calculating OEs from the position and velocity vectors (prediction and oscullating)
n = len(dt)
oeDiff = np.zeros((n, 7))
oePred = np.zeros((n,7))
oeOsc = np.zeros((n, 7))
for k in range(n):
    rD = np.array([dx[k], dy[k], dz[k]])
    vD = np.array([dxdot[k], dydot[k], dzdot[k]])
    oesD = prop.rv2oe(rD, vD)
    oeDiff[k,:] = oesD

    rO = np.array([ox[k], oy[k], oz[k]])
    vO = np.array([oxdot[k], oydot[k], ozdot[k]])
    oesO = prop.rv2oe(rO, vO)
    oeOsc[k,:] = oesO

    rP = np.array([px[k], py[k], pz[k]])
    vP = np.array([pxdot[k], pydot[k], pzdot[k]])
    oesP = prop.rv2oe(rP, vP)
    oePred[k,:] = oesP

aDel = oeOsc[:,0] - oePred[:,0]
eDel = oeOsc[:,1] - oePred[:,1]
iDel = oeOsc[:,2] - oePred[:,2]
aopDel = oeOsc[:,3] - oePred[:,3]
RAANDel = oeOsc[:,4] - oePred[:,4]
MDel = oeOsc[:,5] - oePred[:,5]
uDel = oeOsc[:,6] - oePred[:,6]

# Plotting difference - oscullating
fig15 = prop.plotLongOEs(pt, aDel, eltitle='Difference of Oscullating/Predicted Semi-Major Axis (SMA)', glabel='SMA', ylabel='del(SMA) (m)')
fig16 = prop.plotLongOEs(pt, eDel, eltitle='Difference of Oscullating/Predicted Eccentricity (Ecc)', glabel='Ecc', ylabel='del(Ecc) Ecc')
fig17 = prop.plotLongOEs(pt, np.degrees(iDel), eltitle='Difference of Oscullating/Predicted Inclination (Inc)', glabel='Inc', ylabel='del(Inc) (°)')
fig18 = prop.plotLongOEs(pt, np.degrees(aopDel), eltitle='Difference of Oscullating/Predicted Argument of Perigee (AOP)', glabel='AOP', ylabel='del(AOP) (°)')
fig19 = prop.plotLongOEs(pt, np.degrees(RAANDel), eltitle='Difference of Oscullating/Predicted Right Acenscion of the Ascending Node (RAAN)', glabel='RAAN', ylabel='del(RAAN) (°)')
fig20 = prop.plotLongOEs(pt, np.degrees(MDel), eltitle='Difference of Oscullating/Predicted Mean Anomaly (M)', glabel='M', ylabel='del(M) (°)')
fig21 = prop.plotLongOEs(pt, np.degrees(uDel), eltitle='Difference of Oscullating/Predicted Argument of Latitude (u)', glabel='u', ylabel='del(u) (°)')

# Metric Used for determing prediction quality - arithmetic average deviation divided by arithmetic average over timespan *100
aDelAvg = 100 * np.mean(aDel) / np.mean(oeOsc[:,0])
eDelAvg = 100 * np.mean(eDel) / np.mean(oeOsc[:,1])
iDelAvg = 100 * np.mean(iDel) / np.mean(oeOsc[:,2])
aopDelAvg = 100 * np.mean(aopDel) / np.mean(oeOsc[:,3])
RAANDelAvg = 100 * np.mean(RAANDel) / np.mean(oeOsc[:,4])
MDelAvg = 100 * np.mean(MDel) / np.mean(oeOsc[:,5])
uDelAvg = 100 * np.mean(uDel) / np.mean(oeOsc[:,6])
print()
print('Arithmetic average of Oscullating - Predicted:')
print(f'SMA: {aDelAvg:.8f}')
print(f'Ecc: {eDelAvg:.8f}')
print(f'Inc: {iDelAvg:.8f}')
print(f'AOP: {aopDelAvg:.8f}')
print(f'RAAN: {RAANDelAvg:.8f}')
print(f'M: {MDelAvg:.8f}')
print(f'u: {uDelAvg:.8f}')

# Converting to RTN
stateOscRTN = np.zeros((len(dt), 3))
statePredRTN = np.zeros((len(dt), 3))
normDiffVal = np.zeros((len(dt), 1))
dataRMS = np.zeros((len(dt), 1))

for k in range(len(dt)):
    rO = np.array([ox[k], oy[k], oz[k]])
    rON = prop.eci2rtn(rO, oeOsc[k,4], oeOsc[k,2], oeOsc[k,6])
    stateOscRTN[k,:] = rON

    rP = np.array([px[k], py[k], pz[k]])
    rPN = prop.eci2rtn(rP, oePred[k,4], oePred[k,2], oePred[k,6])
    statePredRTN[k,:] = rPN

    normDiffVal[k] = np.linalg.norm(rON - rPN)
    dataRMS[k] = np.linalg.norm(rON)

rDiff = stateOscRTN - statePredRTN
fig22 = prop.plotRTNDiff(pt, rDiff, eltitle='Position Differences of Oscullating/Predicted in RTN', ylabel='Position (m)')

rDataOsc = np.array([np.mean(statePredRTN[0,:]), np.mean(statePredRTN[1,:]), np.mean(statePredRTN[2,:])])
rMeanDiff = np.array([np.mean(rDiff[0,:]), np.mean(rDiff[1,:]), np.mean(rDiff[2,:])])
rMeanDiffTotal = np.mean(normDiffVal)
rDevDiff = np.array([np.std(rDiff[0,:]), np.std(rDiff[1,:]), np.std(rDiff[2,:])])
rDevDiffTotal = np.std(normDiffVal)

print("RTN Osculating mean from data: ")
print(rDataOsc)
print(f"Mean RMS: {np.mean(dataRMS):.4f}")
print("Mean vector of the difference: ")
print(rMeanDiff)
print(f"Mean RMSE: {rMeanDiffTotal:.4f}")
print("Standard deviation vector of the difference: ")
print(rDevDiff)
print(f"Standard deviation of the RMSE: {rDevDiffTotal:.4f}")


## HW 2 Problem 3
print()
print("="*70)
print("Homework 2 Problem 3")
print("="*70)

# To estimate abar2, going to take the amplitude from peak-to-peak variation with one orbit's data
abar2guess = np.zeros(len(eqCrossInd)-1)
orbs = np.zeros(len(eqCrossInd)-1)
for k in range(len(eqCrossInd)-1):
    a1orb = a[eqCrossInd[k]:eqCrossInd[k+1]]
    a1orbmax = max(a1orb)
    a1orbmin = min(a1orb)
    abar2guess[k] = (a1orbmax - a1orbmin) / 2
    orbs[k] = k + 1
abar2P1 = np.mean(abar2guess)
a1bar2P1 = np.mean(abar2guess[0:15])
a2bar2P1 = np.mean(abar2guess[len(abar2guess)-15:len(abar2guess)])
fig23 = prop.plotaBar(orbs, abar2guess, abar2P1)
print(f'The mean estimate for abar2 over entire year (the twice-per-rev osc. element): {abar2P1:.8f} m')

# Calculating abar2 for first and last day given analytic equations
a1avg = np.mean(aAvg[0:15])
e1avg = np.mean(eAvg[0:15])
i1avg = np.mean(iAvg[0:15])
a2avg = np.mean(aAvg[len(aAvg)-15:len(aAvg)])
e2avg = np.mean(eAvg[len(aAvg)-15:len(aAvg)])
i2avg = np.mean(iAvg[len(aAvg)-15:len(aAvg)])

n1bar = np.sqrt(eMu / (a1avg**3))
n2bar = np.sqrt(eMu / (a2avg**3))

omegabardot1 = -0.75 * n1bar * (eRad / a1avg)**2 * j2 * (1/(1 - e1avg**2)**2) * (1 - 5*(np.cos(i1avg))**2)
omegabardot2 = -0.75 * n2bar * (eRad / a2avg)**2 * j2 * (1/(1 - e2avg**2)**2) * (1 - 5*(np.cos(i2avg))**2)

Mbardot1 = n1bar * (1 - 0.75*(eRad/a1avg)**2 * j2 * (1/(1 - e1avg**2)**1.5) * (1 - 3*(np.cos(i1avg))**2))
Mbardot2 = n2bar * (1 - 0.75*(eRad/a2avg)**2 * j2 * (1/(1 - e2avg**2)**1.5) * (1 - 3*(np.cos(i2avg))**2))

ubardot1 = omegabardot1 + Mbardot1
ubardot2 = omegabardot2 + Mbardot2

a1bar2 = 0.75 * np.sqrt(eMu) * j2 * (eRad**2) * (((a1avg**(-5/2))*(1 - np.cos(2*i1avg)))/(ubardot1))
a2bar2 = 0.75 * np.sqrt(eMu) * j2 * (eRad**2) * (((a2avg**(-5/2))*(1 - np.cos(2*i2avg)))/(ubardot2))

print(f"First day abar, ebar, ibar: {a1avg:.4f}, {e1avg:.4f}, {np.rad2deg(i1avg):.4f}")
print(f"First day abar2 guess: {a1bar2P1:.4f}")
print(f"First day abar2 error: {(a1bar2P1 - a1bar2):.4f}")
print(f"Last day abar, ebar, ibar: {a2avg:.4f}, {e2avg:.4f}, {np.rad2deg(i2avg):.4f}")
print(f"Last day abar2 guess: {a2bar2P1:.4f}")
print(f"Last day abar2 error: {(a2bar2P1 - a2bar2):.4f}")


##### Comment out if you want to SHOW the plot #####
## Problem 1 Plots
plt.close(fig1)
plt.close(fig2)
plt.close(fig3)
plt.close(fig4)
plt.close(fig5)
plt.close(fig6)
plt.close(fig7)
plt.close(fig8)
plt.close(fig9)
plt.close(fig10)
plt.close(fig11)
plt.close(fig12)
plt.close(fig13)
plt.close(fig14)
plt.close(fig24)
plt.close(fig25)
plt.close(fig26)
plt.close(fig27)
plt.close(fig28)
plt.close(fig29)
plt.close(fig30)
## Problem 2 Plots
plt.close(fig15)
plt.close(fig16)
plt.close(fig17)
plt.close(fig18)
plt.close(fig19)
plt.close(fig20)
plt.close(fig21)
plt.close(fig22)
# Problem 3 Plots
plt.close(fig23)
plt.show()

