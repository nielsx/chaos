# Solving a system of 4 ODE using 4th-Order Runge-Kutta
# (Cf. Chapra, "Numerical Methods..." 6th. Ed. p.733,738)
# Application to Kepler Problem Perturbed
# (Cf. Goldstein, "Classical Mechanics" 3rd. Ed. p.495)
# Version: 2019-02-14 20:34
# Running in Python 3
# ns137.wordpress.com

import numpy as np  # cmd.exe as admin: pip install numpy
import matplotlib.pyplot as plt  # cmd.exe as admin: pip install numpy
from matplotlib.ticker import MultipleLocator

###################### FUNCTIONS

# calculates y1, y2, y3, y4 using Runge-Kutta
# e.g. in kepler problem with perturbation, using Hamiltonian formalism
# positions (x,y) and momenta (px,py) forms the set (x, y, px, py),
# and in Python we are going to use: (y1, y2, y3, y4) respectively.
def solve_ode(formula1, formula2, formula3, formula4, x_range,\
			  y10, y20, y30, y40):
	h = x_range[1]-x_range[0]
	y1_range = [None] * (len(x_range))
	y2_range = [None] * (len(x_range))
	y3_range = [None] * (len(x_range))
	y4_range = [None] * (len(x_range))
	y1_range[0] = y10
	y2_range[0] = y20
	y3_range[0] = y30
	y4_range[0] = y40
	for i in range(1,len(x_range)):
		xi = x_range[i-1]
		y1i = y1_range[i-1]
		y2i = y2_range[i-1]
		y3i = y3_range[i-1]
		y4i = y4_range[i-1]
		x = xi
		y1 = y1i
		y2 = y2i
		y3 = y3i
		y4 = y4i
		ka1 = eval(formula1)
		ka2 = eval(formula2)
		ka3 = eval(formula3)
		ka4 = eval(formula4)
		x = xi + h/2
		y1 = y1i + ka1*h/2
		y2 = y2i + ka2*h/2
		y3 = y3i + ka3*h/2
		y4 = y4i + ka4*h/2
		kb1 = eval(formula1)
		kb2 = eval(formula2)
		kb3 = eval(formula3)
		kb4 = eval(formula4)
		x = xi + h/2
		y1 = y1i + kb1*h/2
		y2 = y2i + kb2*h/2
		y3 = y3i + kb3*h/2
		y4 = y4i + kb4*h/2
		kc1 = eval(formula1)
		kc2 = eval(formula2)
		kc3 = eval(formula3)
		kc4 = eval(formula4)
		x = xi + h
		y1 = y1i + kc1*h
		y2 = y2i + kc2*h
		y3 = y3i + kc3*h
		y4 = y4i + kc4*h
		kd1 = eval(formula1)
		kd2 = eval(formula2)
		kd3 = eval(formula3)
		kd4 = eval(formula4)
		y1_range[i]= y1_range[i-1] + (h/6)*(ka1+2*kb1+2*kc1+kd1)
		y2_range[i]= y2_range[i-1] + (h/6)*(ka2+2*kb2+2*kc2+kd2)
		y3_range[i]= y3_range[i-1] + (h/6)*(ka3+2*kb3+2*kc3+kd3)
		y4_range[i]= y4_range[i-1] + (h/6)*(ka4+2*kb4+2*kc4+kd4)
	return y1_range, y2_range, y3_range, y4_range

def add_plot(x_range, y_range):
    plt.plot(np.array(x_range), np.array(y_range))
    
def add_plotf(x_range, formula):  
    x = np.array(x_range)  
    y = eval(formula)
    plt.plot(x, y, linestyle='--', color='red')

def config_plot(title,label_x,label_y,lim_x,lim_y,labels):
    if(len(labels)>0): plt.legend(labels, loc='lower left')
    fig.suptitle(title, fontsize=10)
    plt.xlabel(label_x, fontsize=10)
    plt.ylabel(label_y, fontsize=10)
    #fig.savefig('plot-function.png')
    axes = plt.gca()
    axes.set_xlim(lim_x)
    axes.set_ylim(lim_y)
    plt.grid()
#    axes.xaxis.set_minor_locator(MultipleLocator(10))
#    axes.yaxis.set_minor_locator(MultipleLocator(10))
#    axes.grid(which = 'minor', color='gray', linestyle=':')
    plt.gcf().subplots_adjust(bottom=0.13, left=0.17)
    plt.show()
    
def fix(numero, numero_decimales = 0): # redondea un numero
	return ("%3."+str(numero_decimales)+"f")%numero

def print_table(a,b,c): # forma una tabla en Latex
	for i in range(0,len(a)):
		print(fix(a[i],4),"&",fix(b[i],4),"&",fix(c[i],4),"\\\\")

###################### MAIN ROUTINE

fig = plt.figure(figsize=(5, 4))

m = 1e2 # earth mass [earth mass]
k = 1e-2 # ?
r1 = 0.1 # perihelio [AU]
E = 1e-2 # energy

X = np.arange(0,50,0.1) # (a,b,c) array [a,b[ step c
Y10 = r1 # initial position x
Y20 = 0 # initial position y
Y30 = 0 # initial momentum px
Y40 = np.sqrt(2*m*E+2*m*k*(+r1)) # initial momentum py
eq1 = 'y3/m'
eq2 = 'y4/m'
eq3 = '-k*y1*(y1**2+y2**2)**(-3/2.0)'
eq4 = '-k*y2*(y1**2+y2**2)**(-3/2.0)'
Y1,Y2,Y3,Y4 = solve_ode(eq1, eq2, eq3, eq4, X, Y10, Y20, Y30, Y40)

print_table(Y2,Y1,Y3)

add_plot(Y1,Y2)
# add_plotf(X,'23.0*(x/19.3)-g*(x/19.3)**2/2')

#legends=[r'$y=y(c,m,D,x)$',
#	r'$y=v_{0y}\left( \dfrac{x}{v_{0x}} \right)-\frac{1}{2}\,
#		g\left( \dfrac{x}{v_{0x}} \right)^2$']

#legends=[r'$y=y(x)$']

config_plot('Perturbed Kepler Problem '+\
	'(m='+str(m)+', k='+str(k)+', r1='+str(r1)+', E='+str(E)+')',\
            r'x | time = np.arange(0,50,0.1)',r'$y$',\
#           [min(Y1),max(Y1)],[min(Y2),max(Y2)],legends)
           [min(Y1),max(Y1)],[min(Y2),max(Y2)],'')
