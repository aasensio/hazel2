import numpy as np

def hsra_continuum(x):
    
    c1 = (-4.906054765549e13,1.684734544039e11,1.507254517567e7,-7561.242976546,0.,0.,0.,0.,0.,0.,0.)
    c2 = (-4.4650822755e14,6.1319780351059e11,-9.350928003805e7,0.,0.,0.,0.,0.,0.,0.,0.)
    c3 = (-1.025961e15,1.3172859e12,-3.873465e8,46486.541,-2.049,0.,0.,0.,0.,0.,0.)
    c4 = (4.861821e15,-2.2589885e12,4.3764376e8,-39279.61444,1.34388,0.,0.,0.,0.,0.,0.)
    c5 = (1.758394e15,-3.293986e11,1.6782617e7,0.,0.,0.,0.,0.,0.,0.,0.)
    c6 = (1.61455557e16,-6.544209e12,1.0159316e9,-70695.58136,1.852022,0.,0.,0.,0.,0.,0.)
    c7 = (7.97805136e14,-1.16906597e11,5.315222e6,-4.57327954,-3.473452e-3,0.,0.,0.,0.,0.,0.)

    if (x < 3644.15):        
        conhsra=c1[0]+x*(c1[1]+x*(c1[2]+x*(c1[3]+x*c1[4])))
    elif (x < 3750.):
        conhsra=c2[0]+x*(c2[1]+x*(c2[2]+x*(c2[3]+x*c2[4])))
    elif (x < 6250.):
        conhsra=c3[0]+x*(c3[1]+x*(c3[2]+x*(c3[3]+x*c3[4])))
    elif (x < 8300.):
        conhsra=c4[0]+x*(c4[1]+x*(c4[2]+x*(c4[3]+x*c4[4])))
    elif (x < 8850.):
        conhsra=c5[0]+x*(c5[1]+x*(c5[2]+x*(c5[3]+x*c5[4])))
    elif (x < 10000.):
        conhsra=c6[0]+x*(c6[1]+x*(c6[2]+x*(c6[3]+x*c6[4])))
    else:
        conhsra=c7[0]+x*(c7[1]+x*(c7[2]+x*(c7[3]+x*c7[4])))

    factor = 29979245800.0 / (x*1e-8)**2
    
    return conhsra / factor
