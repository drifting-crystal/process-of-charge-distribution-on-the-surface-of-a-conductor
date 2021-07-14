# Finally put the segmentation function in the C directory
%matplotlib notebook
import math
import random as rand
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# This the initial condition
fig = plt.figure()

C = {
    'mu' : 1e-6, # kq1q2/m1
    'dT': 1,
    'eT': 0,
    'tLim': 300,
    'err' : 1e-4,
    'xbound': 2,
    'e-num': 2, # electrons/particles number
    'vx_max': 1e-4,
    'vy_max': 1e-4
}



aList = []  # active
dList = []


def cullList():
    temp = []
    for a in aList:
        if a.active == True:
            temp.append(a)
        else:
            dList.append(a)
    aList.clear()
    for a in temp:
        aList.append(a)

        
        
def derivative(ffunc,a,method='central',h=1e-4):
    if method == 'central':
        return (ffunc(a + h) - ffunc(a - h))/(2*h)
    elif method == 'forward':
        return (ffunc(a + h) - ffunc(a))/h
    elif method == 'backward':
        return (ffunc(a) - ffunc(a - h))/h
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")
        
        
def anyline(x, afunc):
    axes.plot(x, afunc, 'k-')
    
# Electronic particles
# Then we begin to set the Electrostatic Force


class body:
    def __init__(self, xy, v):
        self.xy = np.array(xy, dtype=float)
        self.v = np.array(v, dtype=float)
        self.rad = 1e-15 * C['xbound']
        self.active = True
        self.xLog = []
        self.yLog = []

    def fE(self):
        if self.active == True:
            for b in aList:
                if b.active == True and b != self:
                    dist = b.xy - self.xy
                    r = np.hypot(dist[0], dist[1])
                    if r > (self.rad + b.rad):
                        aE = -(C['mu'] / r**2) * (dist / r)
                        self.v += aE * C['dT']
                    elif r <= abs(self.rad + b.rad):
                        b.active = False
                        self.active = False
                        print("Collision")

    def translate(self):
        if self.active == True:
            self.xy += self.v * C['dT']
            self.xLog.append(self.xy[0])
            self.yLog.append(self.xy[1])  
                      
    def boundcheck(self, ffunc):
        pdist = abs(self.xy[1] - ffunc(self.xy[0]))
        if pdist <= C['err']:
            self.xy[1] = ffunc(self.xy[0])
            slope = derivative(ffunc, self.xy[0])
            theta1 = math.atan(slope)
            theta3 = math.atan(abs(self.v[1] / self.v[0]))
            theta2 = (theta3 - theta1)
            if self.v[1] == 0:
                self.v[0] = (self.v[0]) * math.cos(theta2) * math.cos(theta1) 
                self.v[1] = (self.v[0]) * math.cos(theta2) * math.sin(theta1) 
            elif self.v[1] != 0:
                self.v[0] = (self.v[1] / (math.sin(theta3))) * math.cos(theta2) * math.cos(theta1) 
                self.v[1] = (self.v[1] / (math.sin(theta3))) * math.cos(theta2) * math.sin(theta1)

            
fig,axes = plt.subplots(1)
line, = axes.plot([], [])

# Line No.1
x1 = np.linspace(-C['xbound'], C['xbound'], 100)
afunc1 = (1 - 0.25 *(x1)**2)**(1/2)
anyline(x1, afunc1)
ffunc1 = lambda x: (1 - 0.25*(x)**2)**(1/2)

# Line No.2 
x2 = np.linspace(-C['xbound'], C['xbound'], 100)
afunc2 = - (1 - 0.25*(x2)**2)**(1/2)
anyline(x2, afunc2)
ffunc2 = lambda x: - (1-0.25*(x)**2)**(1/2)
axes.set_aspect(1)

# set the electrons
for i in range(0, C['e-num']):    
    xOrigin = rand.uniform(-C['xbound'], C['xbound'])
    yOrigin = rand.uniform(ffunc2(xOrigin), ffunc1(xOrigin))
    vx = rand.uniform(-C['vx_max'], C['vx_max'])
    vy = rand.uniform(-C['vy_max'], C['vy_max'])
    aList.append(body([xOrigin, yOrigin], [vx, vy]))
    plt.plot(xOrigin, yOrigin, "b.")
    
# # Test Code
# # No.1
# aList.append(body([0.5, 0.5], [0,0.3]))
# plt.plot(0.5, 0.5, "b.")
# # No.2
# aList.append(body([-0.5, 0.5], [0,0.3]))
# plt.plot(-0.5, 0.5, "b.")





def init():
    """initialize animation"""
    line.set_data([], [])
    return line

def animate(frame):
    i = 0
    while i <= 1/C['dT'] and len(aList) > 1: # fps
        i += 1
        """perform animation step"""
        cullList()
        for a in aList:
            a.fE()
        for a in aList:
            a.translate()
            a.boundcheck(ffunc1)
            a.boundcheck(ffunc2)
            
    for a in aList:
        x = a.xLog
        y = a.yLog
        plt.plot(x, y)
        line.set_data(x[:frame], y[:frame])
    for a in dList:
        x = a.xLog
        y = a.yLog
        plt.plot(x, y, ls='--')   

ani = FuncAnimation(fig, animate, frames=300,
                              interval=20, blit=True, init_func=init)        
        
from IPython.display import HTML
HTML(ani.to_jshtml()) # javascript-HTML Embed animation        
