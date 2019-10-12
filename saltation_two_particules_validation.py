# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 20:38:21 2019

@author: Karin
"""

from matplotlib.pylab import *


# Unidades base SI (m, kg, s)
_m = 1.
_kg = 1.
_s = 1.
_mm = 1e-3*_m
_cm = 1e-2*_m
_gr = 1e-3*_kg 

# velocidades iniciales Particual 1
vfx = 10.0*_m/_s
vfy = 0.1*_m/_s

# velocidades iniciales Particual 2
vfxdos = 5.0*_m/_s
vfydos = 0.1*_m/_s

# Posicion y velocidad actual de la particula 1
x0 = array([0., 1.*_mm], dtype=double)
v0 = array([1., 1.], dtype=double) 

# Posicion y velocidad actual de la particula 2
x0dos = array([0., 20.*_mm], dtype=double) # Parte un instante despues
v0dos = array([1., 2.], dtype=double)

# Particula 1
xi = x0      # posicion actual
vi = v0      # velocidad actual
xim1 = zeros(2, dtype=double)   # posicion siguiente
vim1 = zeros(2, dtype=double)   # velocidad siguiente

# Particula 2
xidos = x0dos      # posicion actual
vidos = v0dos     # velocidad actual
xim1dos = zeros(2, dtype=double)   # posicion siguiente
vim1dos = zeros(2, dtype=double)   # velocidad siguiente


g = 9.81*_m/(_s**2)       # gravedad
d = 1*_mm			# Diametro de la particula
rho_agua = 1000.*_kg/(_m**3)    # densidad del agua
rho_particula = 2650*_kg/(_m**3)      # densidad de la particula, considerando que sea arena 
Cd = 0.47                   # coeficiente de Drag para particula esferica

A = pi*(d/2)**2			# Area de la particula
V = (4./3.)*pi*(d/2)**3	# Volumen de la particula
m = rho_particula*V  # masa de la particula

dt = 0.001*_s    # paso de tiempo 
tmax = 2.*_s    # tiempo maximo de simulacion
ti = 0.*_s     # tiempo actual

W = array([0., -m*g])  #Fuerza de peso (junto con gravedad)
fB = array([0,rho_agua*V*g]) # Fuerza boyante

t = arange(0,tmax,dt)
Nt = len(t)

norm = lambda v: sqrt(dot(v,v))


k_penal = 1000.*0.5*Cd*rho_agua*A*norm(v0)/(1*_mm)


lista_particula1=[]



# Funcion para la particula 1
def particula(z,t):
	xi = z[:2]
	vi = z[2:]
	vf = array([vfx, vfy])
	vrel = vf - vi
	fD = (0.5*Cd*rho_agua*norm(vrel)*A)*vrel 
	# fL = 3.0/4.0*alpha*Cd*()
	Fi = W + fD + fB

	if xi[1] < 0:
		Fi[1] += -k_penal*xi[1] 



	zp = zeros(4)
	zp[:2] = vi


	zp[2:] = Fi/m
	lista_particula1.append(xi)
	return zp

# Funcion para la particula 2
def particula2(z,t):
	xidos = z[:2]
	vidos = z[2:]
	vfdos = array([vfxdos, vfydos])
	vreldos = vfdos - vidos
	fD = (0.5*Cd*rho_agua*norm(vreldos)*A)*vreldos 
	# fL = 3.0/4.0*alpha*Cd*()
	Fi = W + fD + fB

	if xidos[1] < 0:
		Fi[1] += -k_penal*xidos[1] 

	zpdos = zeros(4)
	zpdos[:2] = vidos


	zpdos[2:] = Fi/m
	return zpdos	

# el vector z tiene la posicion y la velocidad de 
# la particula

from scipy.integrate import odeint

# Variables para particula 2
z0dos = zeros(4)
z0dos[:2] = x0dos
z0dos[2:] = v0dos
zdos = odeint(particula2, z0dos, t)
xdos = zdos[:,:2] # definir x igual a z
vdos = zdos[:,2:]

# Particula 1
z0 = zeros(4)
z0[:2] = x0
z0[2:] = v0
z = odeint(particula, z0, t)
x = z[:,:2] # definir x igual a z
v = z[:,2:]

# Figura de movimiento de ambas particulas superpuestas, moviendose tanto en el eje X como en el eje Y
figure()
plot(x[:,0],x[:,1], label="Particula1")
ylim([0,10*_mm])
plot(xdos[:,0],xdos[:,1], label="Particula2")
ylim([0,10*_mm])
plt.xlabel("Avance Particula (mm)")
plt.ylabel("Altura Particula (mm)")
plt.title("Movimiento de 2 particulas")
plt.legend()

# Movimiento de la particula en el tiempo
figure()
subplot(2,1,1)
plot(t,x[:,0], label="x1")
plot(t,x[:,1], label="y1")
plot(t,xdos[:,0], label="x2")
plot(t,xdos[:,1], label="y2")
plt.title("Movimiento y Velocidad de 2 Particulas en funcion del tiempo")
plt.xlabel("Tiempo (s)")
plt.ylabel("Distancia (mm)")
plt.legend()

# Velocidad de la particula en el tiempo, Aceleracion
subplot(2,1,2)
plt.plot(t,v[:,0], label="vx1")
plt.plot(t,v[:,1], label="vy1")
plot(t,vdos[:,0], label="vx2")
plot(t,vdos[:,1], label="vy2")
plt.xlabel("Tiempo (s)")
plt.ylabel("Velocidad (m/s)")
legend()
show()


© 2019 GitHub, Inc.