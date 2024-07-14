import math
import numpy as np
import matplotlib.pyplot as plt
import celestial_constants as cc
from celestial_constants import timestep

def Times(steps: int) -> np.ndarray:
    times = np.zeros((steps+1,1))
    for i in range(steps+1):
        times[i] = i * timestep
    return times
def Force(m1: float, m2: float, rj: np.ndarray, ri: np.ndarray) -> np.ndarray:
    assert len(rj) == len(ri) == 2
    return -cc.G * m1 * m2 * (rj-ri) / (((rj[0]-ri[0])**2 + (rj[1]-ri[1])**2)**(3/2))

def Positions(Perihelion: float, PlanetMass: float, SunMass: float, Velocity: float, steps: int) -> np.array:
    positions = np.zeros((steps+1, 2))
    positions[0] = np.array([Perihelion, 0])
    positions[1] = np.array([Perihelion - ((1/2) * ((cc.G * SunMass) / (Perihelion**2)) * (timestep**2)), Velocity * timestep])
    
    # Second order difference equation
    for l in range(1, steps):
        positions[l+1,:] = (2 * positions[l,:]) - (positions[l-1,:]) + (Force(PlanetMass, SunMass, positions[l,:], [0,0]) * timestep**2 / PlanetMass)
    
    return positions

def Velocity_0(ra: float,rp: float, period: float) -> float:
    return (math.pi * (ra+rp) * math.sqrt(ra*rp))/(rp*period)

def Accelerations(Mj: float, Mi: float, positions_1: np.ndarray, positions_2: np.ndarray = None) -> np.array:
    if positions_2 is None:
        positions_2 = np.zeros((len(positions_1), 2))
        
    assert positions_1.shape == positions_2.shape
    
    accelerations = np.zeros((len(positions_1), 2))
    
    for i in range(len(positions_1)):
        accelerations[i] = Force(Mj, Mi, positions_1[i], positions_2[i]) / Mj
        
    return accelerations

def Differences(positions_1: np.ndarray, positions_2: np.ndarray) -> np.ndarray:
    assert positions_1.shape == positions_2.shape
    return positions_1 - positions_2

def simulate_orbits(positionsE: np.ndarray, positionsM: np.ndarray, positionsJ: np.ndarray, simulate = True) -> None:
    plt.figure(figsize=(6.5, 6.5))
    plt.plot(positionsE[:, 0], positionsE[:, 1], 'b-', label='Earth Orbit')
    plt.plot(positionsM[:, 0], positionsM[:, 1], 'r-', label='Mars Orbit')
    plt.plot(positionsJ[:, 0], positionsJ[:, 1], 'g-', label='Jupiter Orbit')
    plt.plot(0, 0, 'ro', markersize=18, markerfacecolor='yellow', markeredgecolor='black', label='Sun')  # Marker for the Sun
    
    h1, = plt.plot([], [], 'bo', markersize=8, markerfacecolor='blue', markeredgecolor='black', label='Earth') # this is for
    h2, = plt.plot([], [], 'ro', markersize=8, markerfacecolor='red', markeredgecolor='black', label='Mars')
    h3, = plt.plot([], [], 'go', markersize=12, markerfacecolor='green', markeredgecolor='black', label='Jupiter')
    
    plt.axis((-8.5e11, 8.5e11, -8.5e11, 8.5e11))
    plt.grid(True)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Orbits of Earth, Mars, and Jupiter around the Sun')
    plt.legend()
        
    h1.set_data([positionsE[0, 0]], [positionsE[0, 1]])
    h2.set_data([positionsM[0, 0]], [positionsM[0, 1]])
    h3.set_data([positionsJ[0, 0]], [positionsJ[0, 1]])  
        
    # Simulate the motion by updating the marker's position
    if simulate:
        plt.pause(2)
        for k in range(1, len(positionsE)):
            h1.set_data([positionsE[k, 0]], [positionsE[k, 1]])
            h2.set_data([positionsM[k, 0]], [positionsM[k, 1]])
            h3.set_data([positionsJ[k, 0]], [positionsJ[k, 1]])
            
            plt.pause(0.01)  # Pause for a short duration to simulate motion
        
    plt.show()
    
def main() -> None:
    ...
    Esteps = math.ceil(cc.EarthO/timestep)
    Msteps = math.ceil(cc.MarsO/timestep)
    Jsteps = math.ceil(cc.JupiterO/timestep)
    
    VEarth = Velocity_0(cc.EarthA, cc.EarthP, cc.EarthO)
    VMars = Velocity_0(cc.MarsA, cc.MarsP, cc.MarsO)
    VJupiter = Velocity_0(cc.JupiterA, cc.JupiterP, cc.JupiterO)
    
    Etimes = Times(Esteps) 
    Mtimes = Times(Msteps)
    Jtimes = Times(Jsteps)
    
    positionsE = Positions(cc.EarthP, cc.EarthMass, cc.SunMass, VEarth, Esteps)
    positionsM = Positions(cc.MarsP, cc.MarsMass, cc.SunMass, VMars, Msteps)
    positionsJ = Positions(cc.JupiterP, cc.JupiterMass, cc.SunMass, VJupiter, Jsteps)
    
    positions3yearsE = Positions(cc.EarthP, cc.EarthMass, cc.SunMass, VEarth, Esteps*3)
    positions3yearsM = Positions(cc.MarsP, cc.MarsMass, cc.SunMass, VMars, Esteps*3)
    positionns3yearsJ = Positions(cc.JupiterP, cc.JupiterMass, cc.SunMass, VJupiter, Esteps*3)
    
    positions36yearsE = Positions(cc.EarthP, cc.EarthMass, cc.SunMass, VEarth, Esteps*36)
    positions36yearsM = Positions(cc.MarsP, cc.MarsMass, cc.SunMass, VMars, Esteps*36)
    positions36yearsJ = Positions(cc.JupiterP, cc.JupiterMass, cc.SunMass, VJupiter, Esteps*36)
    
    accelerationsE = Accelerations(cc.EarthMass, cc.SunMass, positionsE)
    accelerationsM = Accelerations(cc.MarsMass, cc.SunMass, positionsM)
    accelerationsJ = Accelerations(cc.JupiterMass, cc.SunMass, positionsJ)
    accelerationsEJ = accelerationsE + Accelerations(cc.EarthMass, cc.JupiterMass, positionsE, positionsJ[:Esteps+1])
    accelerationsMJ = accelerationsM + Accelerations(cc.MarsMass, cc.JupiterMass, positionsM, positionsJ[:Msteps+1])
    
    differencesEJ = Differences(positionsE, positionsJ[:Esteps+1])
    differencesMJ = Differences(positionsM, positionsJ[:Msteps+1])
    
    # simulate_orbits(positions3yearsE, positions3yearsM, positions3yearsJ, simulate = False)
    # simulate_orbits(positions36yearsE, positions36yearsM, positions36yearsJ, simulate = False)    
    
    # Simulation of the orbits
    # simulate_orbits(positions3yearsE, positions3yearsM, positions3yearsJ)
    simulate_orbits(positions36yearsE, positions36yearsM, positions36yearsJ)
    

if __name__ == '__main__':
    main()
    