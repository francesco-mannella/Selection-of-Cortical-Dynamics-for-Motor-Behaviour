from pylab import *
    
X  = array([])
Y  = array([])

b=1.0
c=1.0
r=0.2
T=48

XX0 = b + r*cos(linspace(-pi*2/3.,pi*2/3.,T/2))
YY0 = c + r*sin(linspace(-pi*2/3.,pi*2/3.,T/2))

r1 = r*sin(pi*2/3.)
b1 = b+r*cos(pi*2/3.)
XX1 = b1 + r1*cos(linspace(pi/2.,-pi/2.,T/2))
YY1 = c + r1*sin(linspace(pi/2.,-pi/2.,T/2))

X = hstack([XX0,XX1])
Y = hstack([YY0,YY1])

plot(X,Y)

W = eye(2)*0.5
X,Y = dot(W,vstack([X,Y]))

plot(X,Y)


