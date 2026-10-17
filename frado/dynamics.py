import numpy as np

class Dynamics:
    """State order: x,vx,y,vy,z,vz. Distances in rg; time in seconds."""
    def __init__(self, scales, radiation):
        self.scales, self.radiation = scales, radiation

    def acceleration(self, state, dusty):
        x,_,y,_,z,_ = state
        r = np.sqrt(x*x+y*y+z*z)
        if r <= 0: raise ValueError('Cloud reached the singular origin')
        grav = -self.scales.gravitational_coefficient * np.array([x,y,z])/r**3
        if dusty:
            rad, _ = self.radiation.evaluate(state)
            return grav+rad
        return grav

    def derivative(self, t, state, dusty):
        ax,ay,az = self.acceleration(state,dusty)
        return np.array([state[1],ax,state[3],ay,state[5],az])
