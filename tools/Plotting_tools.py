import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import fsolve

#############################################
#Functions for calculatoins
#############################################

#Define some masses that I might want.
me = 0.000511;
mmu = 0.10566;
mN = 0.93892;
mU = 0.9314941024;
m_1H = 1.007825031898 * mU - me;
m_2H = 2.01410178 * mU - me;
m_3H = 3.01604928199 * mU - me;
m_3He = 3.0160293 * mU - 2*me;

#Short hand for 1/sqrt(3).
i3= 1 / np.sqrt(3)

#The calculation is made in terms of
#q and omega. We need these in terms
#of Q2 and xB. This is given here.
def omega(Q2,xB):
    return Q2/(2*mN*xB)

def q(Q2,xB):
    return np.sqrt(np.power(omega(Q2,xB),2) + Q2)


#This function returns an array of
#xB values for an array of momentum
#values. We also have to set Q2
#and the configuration with f1, f2,
#and f3. We are specifically
#calculating the upper boundary of
#xB.
def get_Upper_xB_Array(Q2,ps,if1,if2,if3):
    xBs=[]
    for p in ps:
        def func(xB):
            pTot = p/if1
            return -m_3He -omega(Q2,xB) +np.sqrt(mN*mN + np.power(if2*pTot,2)) +np.sqrt(mN*mN + np.power(if3*pTot,2)) + np.sqrt(mN*mN + np.power(p-q(Q2,xB),2))
        x_guess= 1 if len(xBs)==0 else xBs[-1]
        xm = fsolve(func,x_guess)
        xBs.append(xm)
    return xBs

#Same as previous function. But
#we return the lower limit of xB.
def get_Lower_xB_Array(Q2,ps,if1,if2,if3):
    xBs=[]
    for p in ps:
        def func(xB):
            pTot = p/if1
            return -m_3He -omega(Q2,xB) +np.sqrt(mN*mN + np.power(if2*pTot,2)) +np.sqrt(mN*mN + np.power(if3*pTot,2)) + np.sqrt(mN*mN + np.power(p+q(Q2,xB),2))
        x_guess= 1 if len(xBs)==0 else xBs[-1]
        xm = fsolve(func,x_guess)
        xBs.append(xm)
    return xBs

#Retrun a color corresponding to the
#max value of the upper curve
def get_xB_Color(Q2,ps,if1,if2,if3):
    maximum = max(get_Upper_xB_Array(Q2,ps,if1,if2,if3))[0]
    c = (maximum-1.0)/2.0
    return plt.cm.jet(c)
#############################################
#Functions for ternary plots
#############################################

#Define functions to change
#coordinate systems.
def lin_to_tern(X,Y):
    f1 = (1 - Y) / 2
    f2 = (np.sqrt(3)*X + Y + 1)/4
    return [f1,f2]

def tern_to_lin(f1,f2):
    Y = 1 - 2*f1
    X = (4*f2 - Y - 1)/np.sqrt(3)
    return [X,Y]

#Adding the configurations to
#the plot.
#Adding the configurations to
#the plot.
def plot_conf(ax,X,Y,f1,f2,mycolor,ct,t1,t2,t3):
    #ct=0.01
    lt = 0.1

    #Drawing the circles and fill
    circle1 = plt.Circle((X+ct, Y), ct, color=mycolor,fill=True)
    circle2 = plt.Circle((X-ct/2, Y+ct*np.sqrt(3)/2), ct, color=mycolor,fill=True)
    circle3 = plt.Circle((X-ct/2, Y-ct*np.sqrt(3)/2), ct, color=mycolor,fill=True)

    circle1F = plt.Circle((X+ct, Y), ct*0.95, color="white",fill=True)
    circle2F = plt.Circle((X-ct/2, Y+ct*np.sqrt(3)/2), ct*0.95, color="white",fill=True)
    circle3F = plt.Circle((X-ct/2, Y-ct*np.sqrt(3)/2), ct*0.95, color="white",fill=True)

    ax.add_patch(circle2)
    ax.add_patch(circle2F)
    ax.add_patch(circle3)
    ax.add_patch(circle3F)
    ax.add_patch(circle1)
    ax.add_patch(circle1F)

    ax.text(X+ct,Y,t1,ha='center',va='center')
    ax.text(X-ct/2,Y+ct*np.sqrt(3)/2,t2,ha='center',va='center')
    ax.text(X-ct/2,Y-ct*np.sqrt(3)/2,t3,ha='center',va='center')
    
    #Drawing the arrows
    f3 = 1 - f1 - f2
    theta12 = math.pi - math.acos((f1*f1 + f2*f2 - f3*f3)/(2*f1*f2))
    #Distance from the center of the 3N shape
    ed=0.5
    f1x=f1
    f1y=0
    f2x=f2*math.cos(theta12)
    f2y=f2*math.sin(theta12)
    f3x=-(f1+f2*math.cos(theta12))
    f3y=-f2*math.sin(theta12)
    ax.arrow(X+((ct)*1)+ct*ed*(f1x/f1),Y                                   ,lt*f1x,lt*f1y,color=mycolor,linewidth=2)
    ax.arrow(X-((ct)/2)+ct*ed*(f2x/f2),Y+(ct)*(np.sqrt(3)/2)+ct*ed*(f2y/f2),lt*f2x,lt*f2y,color=mycolor,linewidth=2)
    ax.arrow(X-((ct)/2)+ct*ed*(f3x/f3),Y-(ct)*(np.sqrt(3)/2)+ct*ed*(f3y/f3),lt*f3x,lt*f3y,color=mycolor,linewidth=2)
    
def plot_conf_2N(ax,X,Y,mycolor,ct,t1,t2):
    lt = 0.1

    #Drawing the circles
    circle1 = plt.Circle((X+(ct*3/4), Y), ct, color=mycolor,fill=True)
    circle2 = plt.Circle((X-(ct*3/4), Y), ct, color=mycolor,fill=True)
    
    circle1F = plt.Circle((X+(ct*3/4), Y), ct*0.95, color="white",fill=True)
    circle2F = plt.Circle((X-(ct*3/4), Y), ct*0.95, color="white",fill=True)
    
    ax.add_patch(circle2)
    ax.add_patch(circle2F)
    ax.add_patch(circle1)
    ax.add_patch(circle1F)

    ax.text(X+(ct*3/4),Y,t1,ha='center',va='center')
    ax.text(X-(ct*3/4),Y,t2,ha='center',va='center')


    ed=0.5
    
    ax.arrow(X+(ct)+ct*ed*0.5,Y, lt*0.5,0,color=mycolor,linewidth=2)
    ax.arrow(X-(ct)-ct*ed*0.5,Y,-lt*0.5,0,color=mycolor,linewidth=2)


def plot_conf_on_tern(ax,f1,f2):
    [X,Y] = tern_to_lin(f1,f2)
    ct = 0.01
    lt = 8 * ct

    #Drawing the circles
    circle1 = plt.Circle((X+ct, Y), ct, color='r',fill=False)
    circle2 = plt.Circle((X-ct/2, Y+ct*np.sqrt(3)/2), ct, color='g',fill=False)
    circle3 = plt.Circle((X-ct/2, Y-ct*np.sqrt(3)/2), ct, color='k',fill=False)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)

    #Drawing the arrows
    f3 = 1 - f1 - f2
    theta12 = math.pi - math.acos((f1*f1 + f2*f2 - f3*f3)/(2*f1*f2))
    ax.arrow(X+(ct)*1,Y+(ct)*0           , lt*f1                       , 0                      ,color="red")
    ax.arrow(X-(ct)/2,Y+(ct)*np.sqrt(3)/2, lt*(f2*math.cos(theta12))   , lt*f2*math.sin(theta12),color='green')
    ax.arrow(X-(ct)/2,Y-(ct)*np.sqrt(3)/2,-lt*(f1+f2*math.cos(theta12)),-lt*f2*math.sin(theta12))
    
#Construct the axes of the plot
def add_axes(ax):

    ax.set_xlim(-i3-0.05,i3+0.05)
    ax.set_ylim(-0.05,1.05)

    ax.axis('off')
    ax.plot([-i3,0,i3,-i3],[0,1,0,0],'k',linewidth=2)
    ds=0.05

    ##Plotting Labels
    ax.text(0.00,-3*ds,'Proton Momentum Fraction',horizontalalignment='center',verticalalignment='center')
    ax.text(0.75*i3,0.75,'Neutron \n Momentum \n Fraction',horizontalalignment='center',verticalalignment='center')
    ax.text(-1*i3,0.75,'Proton \n Momentum \n Fraction \n $p_{1}/(p_{1}+p_{2}+p_{3})$',horizontalalignment='center',verticalalignment='center')
    
    ##Plotting Ticks
    ax.text(-1.5*ds-0.00*i3,1.00,'0%',horizontalalignment='right',verticalalignment='center')
    ax.text(-1.5*ds-0.25*i3,0.75,'12.5%',horizontalalignment='right',verticalalignment='center')
    ax.text(-1.5*ds-0.50*i3,0.50,'25%',horizontalalignment='right',verticalalignment='center')
    ax.text(-1.5*ds-0.75*i3,0.25,'37.5%',horizontalalignment='right',verticalalignment='center')
    ax.text(-1.5*ds-1.00*i3,0.00,'50%',horizontalalignment='right',verticalalignment='center')
    
    ax.text(0.00*i3+ds*0.5*1.5,1.00+ds*0.5*np.sqrt(3)*1.5,'50%',horizontalalignment='center',verticalalignment='bottom')
    ax.text(0.25*i3+ds*0.5*1.5,0.75+ds*0.5*np.sqrt(3)*1.5,'37.5%',horizontalalignment='center',verticalalignment='bottom')
    ax.text(0.50*i3+ds*0.5*1.5,0.50+ds*0.5*np.sqrt(3)*1.5,'25%',horizontalalignment='center',verticalalignment='bottom')
    ax.text(0.75*i3+ds*0.5*1.5,0.25+ds*0.5*np.sqrt(3)*1.5,'12.5%',horizontalalignment='center',verticalalignment='bottom')
    ax.text(1.00*i3+ds*0.5*1.5,0.00+ds*0.5*np.sqrt(3)*1.5,'0%',horizontalalignment='center',verticalalignment='bottom')
    
    ax.text(-1.00*i3+ds*0.5*1.5,-ds*0.5*np.sqrt(3)*1.5,'0%',horizontalalignment='center',verticalalignment='top')
    ax.text(-0.50*i3+ds*0.5*1.5,-ds*0.5*np.sqrt(3)*1.5,'12.5%',horizontalalignment='center',verticalalignment='top')
    ax.text( 0.00*i3+ds*0.5*1.5,-ds*0.5*np.sqrt(3)*1.5,'25%',horizontalalignment='center',verticalalignment='top')
    ax.text( 0.50*i3+ds*0.5*1.5,-ds*0.5*np.sqrt(3)*1.5,'37.5%',horizontalalignment='center',verticalalignment='top')
    ax.text( 1.00*i3+ds*0.5*1.5,-ds*0.5*np.sqrt(3)*1.5,'50%',horizontalalignment='center',verticalalignment='top')
    
    ##Plotting Grid Lines
    ax.plot([-ds-0.00*i3,0.00*i3],[1.00,1.00],'m:',linewidth=0.75)
    ax.plot([-ds-0.25*i3,0.25*i3],[0.75,0.75],'m:',linewidth=0.75)
    ax.plot([-ds-0.50*i3,0.50*i3],[0.50,0.50],'m--',linewidth=0.75)
    ax.plot([-ds-0.75*i3,0.75*i3],[0.25,0.25],'m:',linewidth=0.75)
    ax.plot([-ds-1.00*i3,1.00*i3],[0.00,0.00],'m:',linewidth=0.75)
    
    ax.plot([-1.00*i3,0.00*i3+ds*0.5],[0.00,1.00+ds*0.5*np.sqrt(3)],'m:',linewidth=0.75)
    ax.plot([-0.50*i3,0.25*i3+ds*0.5],[0.00,0.75+ds*0.5*np.sqrt(3)],'m:',linewidth=0.75)
    ax.plot([ 0.00*i3,0.50*i3+ds*0.5],[0.00,0.50+ds*0.5*np.sqrt(3)],'m--',linewidth=0.75)
    ax.plot([ 0.50*i3,0.75*i3+ds*0.5],[0.00,0.25+ds*0.5*np.sqrt(3)],'m:',linewidth=0.75)
    ax.plot([ 1.00*i3,1.00*i3+ds*0.5],[0.00,0.00+ds*0.5*np.sqrt(3)],'m:',linewidth=0.75)
    
    ax.plot([-1.00*i3+ds*0.5,-1.00*i3],[-ds*0.5*np.sqrt(3),0.00],'m:',linewidth=0.75)
    ax.plot([-0.50*i3+ds*0.5,-0.75*i3],[-ds*0.5*np.sqrt(3),0.25],'m:',linewidth=0.75)
    ax.plot([ 0.00*i3+ds*0.5,-0.50*i3],[-ds*0.5*np.sqrt(3),0.50],'m--',linewidth=0.75)
    ax.plot([ 0.50*i3+ds*0.5,-0.25*i3],[-ds*0.5*np.sqrt(3),0.75],'m:',linewidth=0.75)
    ax.plot([ 1.00*i3+ds*0.5,-0.00*i3],[-ds*0.5*np.sqrt(3),1.00],'m:',linewidth=0.75)
