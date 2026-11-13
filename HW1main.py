# HW1 Main File for LEO for EO
#
# Alec Lewandowski
# Last Edited 1/21/2026

import numpy as np
import matplotlib.pyplot as plt
from twoBodyClass import twoBodyClass

# Starting epoch is the GRACE Follow-on mission
epochStart = "2023-07-31T23:59:42" # in UTC time
epochR = np.array([0.60938619264713E+07, -0.18403935709873E+07,  0.25900900321043E+07]) # in meters in ECI
epochV = np.array([-0.27098701192949E+04,  0.97560581475608E+03,  0.70476044362125E+04]) # in meters per second in ECI

# Defining the initial Earth constants as described in problem statement
eMu = 3.986004415E+14
eRad = 6378136.3
eAngVel = 7.292115E-05

# Initializing propagator class
prop = twoBodyClass(eMu, eRad, eAngVel, epochStart)

## HW 1 Problem 1
print()
print("="*70)
print("Homework 1 Problem 1")
print("="*70)
oes0 = prop.rv2oe(epochR, epochV)
print("Initial Orbital Elements:")
print(f"  a = {oes0[0]:.6f} m")
print(f"  e = {oes0[1]:.6f}")
print(f"  i = {np.degrees(oes0[2]):.6f}°")
print(f"  Ω = {np.degrees(oes0[3]):.6f}°")
print(f"  ω = {np.degrees(oes0[4]):.6f}°")
print(f"  M = {np.degrees(oes0[5]):.6f}°")
print()

## HW 1 Problem 2
print()
print("="*70)
print("Homework 1 Problem 2")
print("="*70)
r0, v0 = prop.oe2rv(oes0)
print("Initial State Vector in ECI Frame converted back to R, V:")
print(f"  rx = {r0[0]:.14e} m")
print(f"  ry = {r0[1]:.14e} m")
print(f"  rz = {r0[2]:.14e} m")
print(f"  vx = {v0[0]:.14e} m/s")
print(f"  vy = {v0[1]:.14e} m/s")
print(f"  vz = {v0[2]:.14e} m/s")
numDiffs = np.array([r0[0] - epochR[0], r0[1] - epochR[1], r0[2] - epochR[2], v0[0] - epochV[0], v0[1] - epochV[1], v0[2] - epochV[2]])
print("Part (a): Initial State Vector Numerical Differences in ECI Frame:")
print(f"  rx = {numDiffs[0]:.14e} m")
print(f"  ry = {numDiffs[1]:.14e} m")
print(f"  rz = {numDiffs[2]:.14e} m")
print(f"  vx = {numDiffs[3]:.14e} m/s")
print(f"  vy = {numDiffs[4]:.14e} m/s")
print(f"  vz = {numDiffs[5]:.14e} m/s")
rNormDiff = np.linalg.norm(r0 - epochR)
vNormDiff = np.linalg.norm(v0 - epochV)
print("Part (b): Initial State Vector Numerical Difference Norms in ECI Frame:")
print(f"  rNorm = {rNormDiff:.14e} m")
print(f"  vNorm = {vNormDiff:.14e} m/s")
print()

## HW 1 Problem 3
print()
print("="*70)
print("Homework 1 Problem 3")
print("="*70)
dtSec = 5.0
tF = 60 * 60 # sec
y0 = np.array([epochR[0], epochR[1], epochR[2], epochV[0], epochV[1], epochV[2]])
# Method 1 - Anomaly solving for propagating orbit
oes = prop.rv2oe(epochR, epochV)
n = np.sqrt(prop.mu / oes[0]**3)
Mnew = oes[5] + n*tF
oes[5] = Mnew
rnew, vnew = prop.oe2rv(oes)
print("Method 1:")
print("Final State Vector in ECI Frame after 60 minutes via Solving Kepler:")
print(f"  rx = {rnew[0]:.14e} m")
print(f"  ry = {rnew[1]:.14e} m")
print(f"  rz = {rnew[2]:.14e} m")
print(f"  vx = {vnew[0]:.14e} m/s")
print(f"  vy = {vnew[1]:.14e} m/s")
print(f"  vz = {vnew[2]:.14e} m/s")

# Method 1 - f and g functions for propagating orbit
rF, vF, checkVal = prop.keplerProp(epochR, epochV, tF)
print("Method 3:")
print("If this check value is equal to 1, then Kepler Prop was a success:")
print(f"  f * gdot - g * fdot = {checkVal:.14f}")
print("Final State Vector in ECI Frame after 60 minutes via Solving Kepler:")
print(f"  rx = {rF[0]:.14e} m")
print(f"  ry = {rF[1]:.14e} m")
print(f"  rz = {rF[2]:.14e} m")
print(f"  vx = {vF[0]:.14e} m/s")
print(f"  vy = {vF[1]:.14e} m/s")
print(f"  vz = {vF[2]:.14e} m/s")

# Method 2 for propagating orbit
t, y = prop.propagate(y0, tF, dtSec)
print("Part (b):")
print("Final State Vector in ECI Frame after 60 minutes via Propagating EoMs:")
print(f"  rx = {y[0,-1]:.14e} m")
print(f"  ry = {y[1,-1]:.14e} m")
print(f"  rz = {y[2,-1]:.14e} m")
print(f"  vx = {y[3,-1]:.14e} m/s")
print(f"  vy = {y[4,-1]:.14e} m/s")
print(f"  vz = {y[5,-1]:.14e} m/s")
numDiffs12 = np.array([y[0,-1] - rF[0], y[1,-1] - rF[1], y[2,-1] - rF[2], y[3,-1] - vF[0], y[4,-1] - vF[1], y[5,-1] - vF[2]])
numDiffs23 = np.array([y[0,-1] - rnew[0], y[1,-1] - rnew[1], y[2,-1] - rnew[2], y[3,-1] - vnew[0], y[4,-1] - vnew[1], y[5,-1] - vnew[2]])
print("Part (c):")
print("Initial State Vector Numerical Differences in ECI Frame:")
print(f"  rx = {numDiffs[0]:.14e} m")
print(f"  ry = {numDiffs[1]:.14e} m")
print(f"  rz = {numDiffs[2]:.14e} m")
print(f"  vx = {numDiffs[3]:.14e} m/s")
print(f"  vy = {numDiffs[4]:.14e} m/s")
print(f"  vz = {numDiffs[5]:.14e} m/s")
rNormDiff12 = np.linalg.norm(numDiffs12[:2])
vNormDiff12 = np.linalg.norm(numDiffs12[3:])
rNormDiff23 = np.linalg.norm(numDiffs23[:2])
vNormDiff23 = np.linalg.norm(numDiffs23[3:])
print("Initial State Vector Numerical Difference Norms in ECI Frame between Methods 2 and 3:")
print(f"  rNorm = {rNormDiff23:.14e} m")
print(f"  vNorm = {vNormDiff23:.14e} m/s")
print("Initial State Vector Numerical Difference Norms in ECI Frame between Methods 1 and 2:")
print(f"  rNorm = {rNormDiff12:.14e} m")
print(f"  vNorm = {vNormDiff12:.14e} m/s")


## HW 1 Problem 5
print()
print("="*70)
print("Homework 1 Problem 5")
print("="*70)
zAng = np.deg2rad(144.0)
yAng = np.deg2rad(23.5)
print("Part (a): GRACE Follow-On Epoch position and velocity in X1Y1Z1 Frame:")
r3, r2 = prop.rotMatP5(zAng, yAng)
r0new = r2 @ (r3 @ epochR)
v0new = r2 @ (r3 @ epochV)
oes0 = prop.rv2oe(r0new, v0new)
print("Initial Orbital Elements:")
print(f"  a = {oes0[0]:.6f} m")
print(f"  e = {oes0[1]:.6f}")
print(f"  i = {np.degrees(oes0[2]):.6f}°")
print(f"  Ω = {np.degrees(oes0[3]):.6f}°")
print(f"  ω = {np.degrees(oes0[4]):.6f}°")
print(f"  M = {np.degrees(oes0[5]):.6f}°")
print("Part (c): Only a Z-axis rotation:")
r0new = (r3 @ epochR)
v0new = (r3 @ epochV)
oes0 = prop.rv2oe(r0new, v0new)
print("Initial Orbital Elements:")
print(f"  a = {oes0[0]:.6f} m")
print(f"  e = {oes0[1]:.6f}")
print(f"  i = {np.degrees(oes0[2]):.6f}°")
print(f"  Ω = {np.degrees(oes0[3]):.6f}°")
print(f"  ω = {np.degrees(oes0[4]):.6f}°")
print(f"  M = {np.degrees(oes0[5]):.6f}°")
print("Part (d): Only a Y-axis rotation:")
r0new = r2 @ epochR
v0new = r2 @ epochV
oes0 = prop.rv2oe(r0new, v0new)
print("Initial Orbital Elements:")
print(f"  a = {oes0[0]:.6f} m")
print(f"  e = {oes0[1]:.6f}")
print(f"  i = {np.degrees(oes0[2]):.6f}°")
print(f"  Ω = {np.degrees(oes0[3]):.6f}°")
print(f"  ω = {np.degrees(oes0[4]):.6f}°")
print(f"  M = {np.degrees(oes0[5]):.6f}°")
print()

## HW 1 Problem 6
print()
print("="*70)
print("Homework 1 Problem 6")
print("="*70)
normER = np.linalg.norm(epochR)
normEV = np.linalg.norm(epochV)
instAngVel = np.cross(epochV, epochR) / normER**2
print("Part (a): Instantaneous Angular Velocity in ECI Frame:")
print(f"  wx = {instAngVel[0]:.14e} rad/s")
print(f"  wy = {instAngVel[1]:.14e} rad/s")
print(f"  wz = {instAngVel[2]:.14e} rad/s")
print(" Instantaneous Angular Velocity in Orbital Plane Frame:")
print(f"  w = {np.linalg.norm(instAngVel):.14e} rad/s")
print("Part (b): Mean Angular Velocity over One Revolution (mean motion IS mean ang velo):")
print(f"  wmean = n = {n:.14e} rad/s")
print("The difference between the instantaneous and mean:")
print(f"  diff = {abs(n - np.linalg.norm(instAngVel)):.14e} rad/s")
print()

## HW 1 Problem 7
print()
print("="*70)
print("Homework 1 Problem 7")
print("="*70)
specE = normEV**2 / 2 - prop.mu / normER
hVec = np.cross(epochR, epochV)
specH = np.linalg.norm(hVec)
eVec = (1/prop.mu)*np.cross(epochV,hVec) - epochR / normER
print("Specific Energy of System:")
print(f" \N{GREEK SMALL LETTER EPSILON} = {specE:.14e} m²/s²")
print("Specific Angular Momentum of System:")
print(f" h = {specH:.14e} m²/s²")
print("Laplace (Eccentricity) Vector of System:")
print(f" e = ({eVec[0]:.14e}, {eVec[1]:.14e}, {eVec[2]:.14e}) ")
normVF = np.linalg.norm(vnew)
normRF = np.linalg.norm(rnew)
specEnew = normVF**2 / 2 - prop.mu / normRF
hVecnew = np.cross(rnew, vnew)
specHnew = np.linalg.norm(hVecnew)
eVecnew = (1/prop.mu)*np.cross(vnew,hVecnew) - rF / normRF
print("The values below are calculated with the position and velocity vectors after 60 minutes:")
print("Specific Energy of System:")
print(f" \N{GREEK SMALL LETTER EPSILON} = {specEnew:.14e} m²/s²")
print("Specific Angular Momentum of System:")
print(f" h = {specHnew:.14e} m²/s²")
print("Laplace (Eccentricity) Vector of System:")
print(f" e = ({eVecnew[0]:.14e}, {eVecnew[1]:.14e}, {eVecnew[2]:.14e}) ")
print("The differences in values for reference:")
print(f" \N{GREEK SMALL LETTER EPSILON} = {abs(specEnew - specE):.14e} m²/s²")
print(f" h = {abs(specHnew-specH):.14e} m²/s²")
print(f" e norm = {np.linalg.norm(eVec - eVecnew):.14e} ")
print("Part (c): To check 5 of the 7 IoMs independence, check linear independence at epoch and 60 minutes later")
print("   The 7 integrals are:")
print("       1-3: Angular momentum components (hx, hy, hz)")
print("       4-6: Laplace vector components (ex, ey, ez)")
print("       7: Energy (ε)")
print("   But only 5 are independent because:")
print("       h dotted with e is the zero vector and magnitude of e can be derived from h and energy")
print("                       h·e = 0")
print("                       |e|² = 1 + 2ε|h|²/μ²")
print("   If 2 algebraic constraints are proven, then therefore 7 IoMs - 2 Constraints = 5 independent IoMs")
print("Checking the two constraints at epoch and epoch + 60 min")
print(f"at epoch: h·e = {np.dot(hVec, eVec):.14e}")
print(f"          |e|² = {(np.linalg.norm(eVec)**2):.14e} = 1 + 2ε|h|²/μ² = {(1 + (2*specE*specH**2 / prop.mu**2)):.14e}")
print(f"at epoch: h·e = {np.dot(hVecnew, eVecnew):.14e}")
print(f"          |e|² = {(np.linalg.norm(eVecnew)**2):.14e} = 1 + 2ε|h|²/μ² = {(1 + (2*specEnew*specHnew**2 / prop.mu**2)):.14e}")
print('Because both dot products are about zero, and the values of the seconds equation are all equal, 2 constraints means 5 independent IoMs!')
print()


