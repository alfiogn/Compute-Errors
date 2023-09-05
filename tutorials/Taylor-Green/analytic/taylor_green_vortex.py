from sympy import *
import re
import subprocess

x,y,z = symbols('x,y,z', real=True)

variables = (x, y, z)

def gradient(f):
    try:
        if f.is_scalar:
            return Matrix([f.diff(v) for v in variables])
    except:
        if f.is_Matrix:
            m = f.shape[0]
            return Matrix([[f[i].diff(v) for v in variables] for i in range(m)])

def divergence(f):
    try:
        if f.is_scalar:
            return sum(Matrix([f.diff(v) for v in variables]))
    except:
        if f.is_Matrix:
            m, n = f.shape
            if n == 1:
                return sum(Matrix([f[i].diff(variables[i]) for i in range(m)]))
            else:
                return Matrix([sum([f[i, j].diff(variables[j]) for j in range(n)]) for i in range(m)])

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

subs = {'cos':'Foam::cos',
        'sin':'Foam::sin',
        'exp':'Foam::exp',
        'log':'Foam::log',
        'pow':'Foam::pow'}

def stringFunc(fun):
    s = ccode(fun, standard='c99')
    return replace_all(s, subs)

x0, x1 = -0.25, 0.25
y0, y1 = -0.25, 0.25

ux = -cos(2*pi*x)*sin(2*pi*y)
uy = sin(2*pi*x)*cos(2*pi*y)
uz = 0
U = Matrix([ux, uy, uz])

# # integrate pressure gradient for F=0
# Gp = -simplify(-divergence(gradient(U) + gradient(U).T))
# a, b = Gp[0], Gp[1]
# A = integrate(a, x)  # int_x p_x  =>  p=k3 + By
# Byy = b - diff(A, y)  # By'=p_y - k3_y=
# By = integrate(Byy, y)  # By=int_y By'
# p = A + By
p = -1/4*(cos(4*pi*x) + cos(4*pi*y))


F = simplify(-divergence((gradient(U) + gradient(U).T)) + (U.T*gradient(U)).T + gradient(p))

gradU = gradient(U)
gradp = gradient(p)

fld = 'include/'
subprocess.call('mkdir -p ' + fld, shell=True)

for i,vi in enumerate(variables):
    open(fld + 'f' + str(vi) + '.H', 'w').write('scalar f' + str(vi) + ' = ' + stringFunc(F[i]) + ';')
    open(fld + 'u' + str(vi) + '.H', 'w').write('scalar u' + str(vi) + ' = ' + stringFunc(U[i]) + ';')
    open(fld + 'gradp' + str(vi) + '.H', 'w').write('scalar gradp' + str(vi) + ' = ' + stringFunc(gradp[i]) + ';')

for i,vi in enumerate(variables):
    for j,vj in enumerate(variables):
        open(fld + 'gradu' + str(vi) + str(vj) + '.H', 'w').write('scalar gradu' + str(vi) + str(vj) + ' = ' + stringFunc(gradU[i, j]) + ';')


open(fld + 'p.H' , 'w').write('scalar p = ' + stringFunc(p) + ';')



# Plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sympy.utilities.lambdify import lambdastr

fig = plt.figure(figsize=(15, 6))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132, projection='3d')
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax1.set_title('Velocity')
ax3.set_aspect('equal')
ax3.set_title('Source')

xx = np.linspace(x0, x1, 21)
yy = np.linspace(y0, y1, 21)
Xm, Ym = np.meshgrid(xx, yy)
X, Y = Xm.ravel(), Ym.ravel()

exec("uxl = " + lambdastr((x, y), ux).replace('math', 'np'))
exec("uyl = " + lambdastr((x, y), uy).replace('math', 'np'))
ax1.quiver(X, Y, uxl(X, Y), uyl(X, Y))
exec("pl = " + lambdastr((x, y), p).replace('math', 'np'))
ax2.plot_surface(Xm, Ym, pl(Xm, Ym), cmap=cm.coolwarm)
exec("fxl = " + lambdastr((x, y), F[0]).replace('math', 'np'))
exec("fyl = " + lambdastr((x, y), F[1]).replace('math', 'np'))
ax3.quiver(X, Y, fxl(X, Y), fyl(X, Y))


plt.show()


