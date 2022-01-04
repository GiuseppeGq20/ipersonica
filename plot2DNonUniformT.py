import matplotlib.pyplot as plt
import numpy as np

x=np.linspace(-5,5,200)
y=np.linspace(-5,5,200)
X,Y=np.meshgrid(x,y)

def f(t):
    return np.exp(-(X-Y*t-1)**2 -(Y+X*t-1)**2)

#initial condition
plt.contourf(X,Y,f(0),cmap='jet',levels=30)
plt.show()

times=np.linspace(0,5,200)
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
for time in times:
    plt.cla()
    #ax.plot_surface(X,Y,f(time))
    plt.contourf(X,Y,f(time), cmap='jet',levels=30)
    plt.draw()
    plt.pause(0.001)
