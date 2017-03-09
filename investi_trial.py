##Script for Assignment 2##
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

def dx_dt(x,t,omega=1.788):
    """function which returns a list of derivatives, dx/dt and dy/dt, for y=dx/dt for the general Van der Pol oscillator."""
    return[x[1], a*(1-x[0]**2)*x[1]-x[0]+b*np.cos(omega*t)]

def lorenz(vec,t,sigma=10.0, rho=28.0,beta=8.0/3):
    """function which returns a list of derivatives, dx/dt, dy/dt and dz/dt for the Lorenz attractor."""
    x=vec[0];y=vec[1];z=vec[2]
    return[sigma*(y-x), (x*(rho-z))-y, (x*y)-(beta*z)]


ts=np.linspace(0.0,40.0,1000)

##Part 2a##
# Graphs showing the evolution of displacement of the Van der Pol oscillator 
# over time t with initial values of a=0.5 and b=0, varying ICs

fig1=plt.figure()
xics=np.array([0.1,1.0,2.0,3.0]); yic=0.0;
a=0.5; b=0.0;
for i in xics:
    xs=odeint(dx_dt,(i,yic),ts)
    plt.plot(ts,xs[:,0],'-') 
plt.xlabel("$t$",fontsize=16); plt.ylabel("$x(t)$",fontsize=16);
fig1.savefig('part_2a.jpg')



##Part 2b)
# Graph showing the phase portrait of the Van der Pol oscillator 
# with initial conditions values of a and b as above, varying ICs

fig2=plt.figure()
for i in xics:
    xs=odeint(dx_dt,(i,yic),ts)
    plt.plot(xs[:,0],xs[:,1])
    plt.scatter(xs[0,0],xs[0,1])
plt.xlabel("$x(t)$",fontsize=16); plt.ylabel("$y(t)$",fontsize=16);

fig2.savefig('part_2b.jpg')



##Part 2c)
# Phase portrait for fixed ICs, differing values of a

fig3=plt.figure()
a_array=np.array([0.1,1.0,2.0,3.0])
for j in a_array:
    a=j
    xs=odeint(dx_dt,(3.0,yic),ts)
    plt.plot(xs[int(len(ts)/2):,0],xs[int(len(ts)/2):,1],label="$a=$"+str(j))
plt.xlabel("$x(t)$",fontsize=16); plt.ylabel("$y(t)$",fontsize=16);
plt.legend()
fig3.savefig('part_2c.jpg')


##Part 3a)
# Time-displacement evolution graph of the Van der Pol oscillator, 
# now with a forcing term in the differential equation
# fixed ICs, a and b values

fig4=plt.figure()
x0=0.0 ; y0=0.0   #initial conditions 
ts=np.linspace(0.0,40.0,10000)
a=3; b=5; omega=1.788
xs=odeint(dx_dt,(x0,y0),ts)
plt.plot(ts,xs[:,0])
plt.xlabel("$t$",fontsize=20); plt.ylabel("$x(t)$",fontsize=20);
fig4.savefig('part_3a.jpg')

# phase portrait for same ODE, ICs and a and b values

fig5=plt.figure()
plt.plot(xs[:,0],xs[:,1])
plt.scatter(xs[0,0],xs[0,1])
plt.xlabel("$x(t)$",fontsize=20); plt.ylabel("$y(t)$",fontsize=20);
fig5.savefig('part_3a_2.jpg')



##Part 3b)
# Poincare section for the the above forced VdP oscillator, same ICs, a and b values

fig6=plt.figure()
k_ary=(2*np.pi*np.arange(500.0,5001.0))/omega ##gives the required time values we need to plot for (xk,yk)
x0=0.0; y0=0.0;
xs=odeint(dx_dt,(x0,y0),k_ary)
plt.plot(xs[:,0],xs[:,1],'o',markersize=2)
plt.xlabel("$x(t)$",fontsize=16); plt.ylabel("$y(t)$",fontsize=16);
fig6.savefig('part_3b.jpg')



# Graph of the evolution of the Lorenz attractor over time with given ICs x0, y0, z0

fig8=plt.figure()
ts=np.linspace(0.0,40.0,100000)
ax=Axes3D(fig8)
x0=0.1; y0=0.0; z0=0.0;
xs=odeint(lorenz,(x0,y0,z0),ts)
ax.plot(xs[:,0],xs[:,1],xs[:,2])
ax.scatter(xs[0,0],xs[0,1],xs[0,2])
ax.set_xlabel("$x(t)$",fontsize=20); ax.set_ylabel("$y(t)$",fontsize=20); ax.set_zlabel("$z(t)$",fontsize=20);
fig8.savefig('lorenz_attractor.jpg')



# Partial time-displacement evolution curves for x and y values separately
 
fig9=plt.figure()
plt.plot(ts,xs[:,0], label="$x(t)$"); plt.plot(ts,xs[:,1],label="$y(t)$"); plt.plot(ts,xs[:,2],label="$z(t)$");
plt.xlabel("$t$",fontsize=16)
plt.legend()
fig9.savefig('lor_time_displacements.jpg')

plt.show()