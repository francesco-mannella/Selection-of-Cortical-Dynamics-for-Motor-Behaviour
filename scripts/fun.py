from pylab import *

ion()

N=480
K= 30
scale = [2.0,3.0,5.0,4.0]

x = linspace(0,2*pi,480)

def fun(x,A=K/3, amp = .5, center= 0.5) :
    f = sin(A*x)
    # f+= cos(4*A*x) 
    # f/= 2.0
    f*= amp/2.0
    f+= center
    return f


def plot_fun(y_1,y_2) :
    f = figure()
    f.clear()

    plot(y_1)
    plot(y_2)

    f.canvas.draw()

q = 0.1

NN =  1
NN1 = 1
TT = 298

y1 = (
        fun(x, amp= 0.8, center= 0.5) +  
        fun(x+pi/4., A=20, amp= 0.8, center= 0.5) +
        fun(x+pi/6., A=30, amp= 0.6, center= 0.6)         
        )/3.
y1 = fun(x, amp= 0.8, center= 0.5) 
y2 = fun(x, A=K, amp= 0.5, center= 0.3 + q*scale[0]) 

savetxt("FUN_1",vstack( [ 
    hstack([zeros([TT,N]),zeros([TT,N]),zeros([TT,N])]),
    vstack([hstack([y1,y1,y1]) for t in xrange(NN) ]),
    vstack([hstack([y2,y2,y2]) for t in xrange(NN1) ]) ]))

plot_fun(y1,y2)

y2 = fun(x, A=K,  amp= 0.5, center= 0.3 + q*scale[1]) 

savetxt("FUN_2",vstack( [ 
    hstack([zeros([TT,N]),zeros([TT,N]),zeros([TT,N])]),
    vstack([hstack([y1,y1,y1]) for t in xrange(NN) ]),
    vstack([hstack([y2,y2,y2]) for t in xrange(NN1) ]) ]))

plot_fun(y1,y2)

y2 = fun(x, A=K,  amp= 0.5, center= 0.3 + q*scale[2]) 

savetxt("FUN_3",vstack( [ 
    hstack([zeros([TT,N]),zeros([TT,N]),zeros([TT,N])]),
    vstack([hstack([y1,y1,y1]) for t in xrange(NN) ]),
    vstack([hstack([y2,y2,y2]) for t in xrange(NN1) ]) ]))

plot_fun(y1,y2)

y2 = fun(x, A=K,  amp= 0.5, center= 0.3 + q*scale[3]) 

savetxt("FUN_4",vstack( [ 
    hstack([zeros([TT,N]),zeros([TT,N]),zeros([TT,N])]),
    vstack([hstack([y1,y1,y1]) for t in xrange(NN) ]),
    vstack([hstack([y2,y2,y2]) for t in xrange(NN1) ]) ]))

plot_fun(y1,y2)

raw_input()

