"""
given : sigmaX, sigmaY, TauXY, rotation angle = theta

sigmaX' = (sigmaX + sigmaY)/2 + (sigmaX - sigmaY)/2*cos(2theta) + TauXY*sin(2theta)

sigmaY' = (sigmaX + sigmaY)/2 - (sigmaX - sigmaY)/2*cos(2*theta) - TauXY*sin(2*theta)

TauX'Y'  = - (sigmaX - sigmaY)/2*sin(2theta) + TauXY*cos(2theta)


MOHR'S CIRCLE

A(sigmaX,TauXY)
B(sigmaY,-TauXY)

center of circle = C = (sigmaX + sigmaY) / 2
radius of circle = Tau(max)

PRINCIPAL STRESSES

sigma(1) = (sigmaX + sigmaY)/2 + (((sigmaX - sigmaY)/2)**2 + (TauXY)**2)**0.5

sigma(2) = (sigmaX + sigmaY)/2 - (((sigmaX - sigmaY)/2)**2 + (TauXY)**2)**0.5

tan(2theta-p-) = 2 * TauXY/(sigmaX - sigmaY)


MAXIMUM SHEAR STRESS
Tau(max,min) = +- (((sigmaX - sigmaY)/2)**2 + (TauXY)**2)**0.5

tan(2theta-s-) = -(sigmaX - sigmaY) / (2 * TauXY)

"""


from math import sin, pi, cos, atan
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import matplotlib as mpl


sigmaX = float(input("SigmaX : "))
sigmaY = float(input("SigmaY : "))
TauXY = float(input("TauXY : "))
theta = float(input("Rotation Angle : "))


def r2d(theta):
    theta = theta * pi /180
    return theta

def d2r(theta):
    theta = theta * 180 / pi
    return theta


def sigmaXprime(sigmaX, sigmaY, theta, TauXY):
    new_sigmaX = (sigmaX + sigmaY)/2 + (sigmaX - sigmaY)/2*cos(r2d(2*theta)) + TauXY*sin(r2d(2*theta)) # convert to degree
    return new_sigmaX


def sigmaYprime(sigmaX, sigmaY, theta, TauXY):
    new_sigmaY=(sigmaX + sigmaY)/2 - (sigmaX - sigmaY)/2*cos(r2d(2*theta)) - TauXY*sin(r2d(2*theta))
    return new_sigmaY


def TauXprimeYprime(sigmaX, sigmaY, theta, TauXY):
    new_TauXY = - (sigmaX - sigmaY)/2*sin(r2d(2*theta)) + TauXY*cos(r2d(2*theta))
    return new_TauXY


def sigmaAVG(sigmaX,sigmaY):
    sigmaAVG = (sigmaX + sigmaY) / 2
    return sigmaAVG


def sigmaPrincipal(sigmaX, sigmaY, TauXY):
    sigma1 = (sigmaX + sigmaY)/2 + (((sigmaX - sigmaY)/2)**2 + (TauXY)**2)**0.5
    sigma2 = (sigmaX + sigmaY)/2 - (((sigmaX - sigmaY)/2)**2 + (TauXY)**2)**0.5
    return sigma1, sigma2

def sigmaPrincipalOrientation(sigmaX, sigmaY, TauXY):
    try:
        theta  = d2r((0.5)*(atan((2 * TauXY/(sigmaX - sigmaY)))))
    except ZeroDivisionError:
        theta = 45.0

    while theta < 0:
        theta += 90
    angles = {"theta1":theta, "theta2":theta + 90}
    return angles

def MaxMinShear(sigmaX, sigmaY, TauXY):
    MaxTauXY = (((sigmaX - sigmaY)/2)**2 + (TauXY)**2)**0.5
    a = {"Max":MaxTauXY, "Min":-1 * MaxTauXY}
    return a

def ShearOrientation(sigmaX,sigmaY,TauXY):
    try:
        theta  = d2r((0.5)*(atan(-(sigmaX - sigmaY)/(2 * TauXY))))
    except ZeroDivisionError:
        theta = 45.0
    while theta < 0:
        theta += 90
    angles = {"theta1":theta, "theta2":theta + 90}
    return angles



sigma_x_prime = sigmaXprime(sigmaX, sigmaY, theta, TauXY)

sigma_y_prime = sigmaYprime(sigmaX, sigmaY, theta, TauXY)

tau_x_prime_y_prime= TauXprimeYprime(sigmaX, sigmaY, theta, TauXY)

sigma_avg= sigmaAVG(sigmaX,sigmaY)

sigma_1,sigma_2 = sigmaPrincipal(sigmaX, sigmaY, TauXY)

sigma_angles = sigmaPrincipalOrientation(sigmaX, sigmaY, TauXY)

shear = MaxMinShear(sigmaX,sigmaY,TauXY)


shear_angles = ShearOrientation(sigmaX, sigmaY, TauXY)


print("Sigma X' : ",round(sigma_x_prime, 2),"MPa\n")
print("Sigma Y' : ",round(sigma_y_prime, 2),"MPa\n")
print("TauX'Y' : ",round(tau_x_prime_y_prime, 2),"MPa\n")
print("Average Sigma : ",round(sigma_avg, 2), "MPa\n")
print("Principal Sigma 1 : ",round(sigma_1, 2),"MPa\n")
print("Principal Sigma 2 : ",round(sigma_2, 2),"MPa\n")
print("Principal Normal Theta 1 : ",round(sigma_angles["theta1"], 2)," degree\n")
print("Principal Normal Theta 2 :",round(sigma_angles["theta2"], 2)," degree\n")
print("Maximum Tau : ",round(shear["Max"], 2),"MPa\n")
print("Minimum Tau : ",round(shear["Min"], 2),"MPa\n")
print("Principal Shear Theta 1 : ",round(shear_angles["theta1"], 2)," degree\n")
print("Principal Shear Theta 2 : ",round(shear_angles["theta2"], 2)," degree\n")

def MohrsCirclePlot(sigmaX, sigmaY, theta, TauXY):
    radians = np.linspace(0,360,361) * (2 * pi /360)
    sigma_points = sigma_avg + shear["Max"] * np.cos(radians)
    tau_points = shear["Max"] * np.sin(radians)
    plt.figure(figsize = [5,5])
    plt.plot(sigma_points,tau_points,)
    plt.plot([sigmaX,sigmaY],[TauXY,-TauXY],color = "r",label = "Before Rotating")
    plt.plot([sigma_avg],[0],color= "black",marker = "o")
    plt.fill_between(sigma_points, tau_points, alpha=0.2)
    plt.title("Mohr's Circle")
    plt.xlabel(r'$\sigma$')
    plt.ylabel(r'$\tau$')
    #plt.axhline(color = 'k')
    #plt.axvline(color = 'k')
    
    plt.plot([sigma_x_prime, sigma_y_prime],[tau_x_prime_y_prime, -tau_x_prime_y_prime],color = "green",label = "After Rotating")
    
    plt.text(sigma_1, 0, round(sigma_1,2), va = 'bottom', ha = 'right', fontsize = 12)
    plt.text(sigma_2, 0, round(sigma_2,2), va = 'bottom', ha = 'left', fontsize = 12)
    plt.text(sigma_avg, shear["Max"], round(shear["Max"],2), va = 'top', ha = 'center', fontsize = 12)
    plt.text(sigma_avg, shear["Min"], round(shear["Min"],2), va = 'bottom', ha = 'center', fontsize = 12)
    plt.legend(loc="lower right")
    plt.grid()
    plt.show()
    return 


MohrsCirclePlot(sigmaX, sigmaY, theta, TauXY)

"""
Plane Transformation
"""

def plane_transformation_no_rotating():

    fig = plt.figure()
    ax = fig.add_subplot(111)


    rectangle = patches.Rectangle((0,0), 5, 5, color="black",  alpha=1)

    ax.add_patch(rectangle)

    plt.xlim(-5, 10)
    plt.ylim(-5, 10)

    plt.grid(True)


    plt.title("Before Rotating")
    if sigmaX > 0:

        ax.annotate("", xytext=(5, 2.5), xy=(8, 2.5),
                    arrowprops=dict(arrowstyle="->"))  ##sigma_x
        ax.text(8,2.5,str(sigmaX)+" Mpa")


        ax.annotate("", xytext=(0, 2.5), xy=(-3, 2.5),
            arrowprops=dict(arrowstyle="->")) ##sigma_negative_x
        ax.text(-5,2.5,str(sigmaX)+" MPa")


    else:
        ax.annotate("", xy=(5, 2.5), xytext=(8, 2.5),
                    arrowprops=dict(arrowstyle="->"))  ##sigma_x
        ax.text(8,2.5,str(sigmaX)+" Mpa")


        ax.annotate("", xy=(0, 2.5), xytext=(-3, 2.5),
            arrowprops=dict(arrowstyle="->")) ##sigma_negative_x
        ax.text(-5,2.5,str(sigmaX)+" MPa")

    
    if sigmaY > 0:
        ax.annotate("", xytext=(2.5, 5), xy=(2.5, 8),
                arrowprops=dict(arrowstyle="->")) ##sigma_y
        ax.text(2,8,str(sigmaY)+" MPa")


        ax.annotate("", xytext=(2.5, 0), xy=(2.5, -3),
                arrowprops=dict(arrowstyle="->")) ##sigma_negative_y 
        ax.text(2,-3.5,str(sigmaY)+" MPa")

    else:
        ax.annotate("", xytext=(2.5, 5), xy=(2.5, 8),
                arrowprops=dict(arrowstyle="->")) ##sigma_y
        ax.text(2,8,str(sigmaY)+" MPa")


        ax.annotate("", xy=(2.5, 0), xytext=(2.5, -3),
                arrowprops=dict(arrowstyle="->")) ##sigma_negative_y 
        ax.text(2,-3.5,str(sigmaY)+" MPa")

    

    if TauXY > 0 :
        ax.annotate("", xytext=(5.5, 1.0), xy=(5.5, 4),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_1 

        ax.annotate("", xytext=(1.25, 5.5), xy=(4.25, 5.5),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_2

        ax.text(5,5,str(TauXY) + " MPa")

        ax.annotate("", xytext=(-0.5, 4), xy=(-0.5, 1),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_negative_1
        
        ax.annotate("", xytext=(4.25, -0.5), xy=(1.25, -0.5),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_negative_2

        ax.text(-1.5,-0.5,str(TauXY) + " MPa")
    
    else:
       
        ax.annotate("", xy=(5.5, 1.0), xytext=(5.5, 4),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_1 

        ax.annotate("", xy=(1.25, 5.5), xytext=(4.25, 5.5),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_2

        ax.text(5,5,str(TauXY) + " MPa")

        ax.annotate("", xy=(-0.5, 4), xytext=(-0.5, 1),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_negative_1
        
        ax.annotate("", xy=(4.25, -0.5), xytext=(1.25, -0.5),
                arrowprops=dict(arrowstyle="->")) ##tau_xy_negative_2

        ax.text(-1.5,-0.5,str(TauXY) + " MPa")

    plt.show()
    return 


plane_transformation_no_rotating()



def plane_transformation_rotating(theta):

    fig2 = plt.figure(figsize=[5,5])
    ax2 = fig2.add_subplot()


    rectangle = patches.Rectangle((0,0), 5, 5, color="black",  alpha=1)

    ax2.add_patch(rectangle)

    plt.xlim((-15, 15))
    plt.xticks(rotation = theta)
    plt.ylim(-15, 15)

    plt.axhline(color = 'k')

    plt.grid(True)


    transformation = mpl.transforms.Affine2D().rotate_deg(theta) + ax2.transData
    rectangle.set_transform(transformation)
    if sigma_x_prime > 0:
        ax2.annotate("", xytext=(5.59 * cos(r2d(26.56 + theta)), 5.59 * sin(r2d(26.56 + theta ))), xy=(7.43* cos(r2d(19.70 + theta)), 7.43 * sin(r2d(19.70 + theta))), arrowprops=dict(arrowstyle="->"))
        ax2.text(7.43* cos(r2d(19.70 + theta)), 7.43 * sin(r2d(19.70 + theta)), str(round(sigma_x_prime,2)) + " MPa")
    else: 
        ax2.annotate("", xytext=(7.43* cos(r2d(19.70 + theta)), 7.43 * sin(r2d(19.70 + theta))), xy=(5.59 * cos(r2d(26.56 + theta)), 5.59 * sin(r2d(26.56 + theta ))), arrowprops=dict(arrowstyle="->"))
        ax2.text(7.43* cos(r2d(19.70 + theta)), 7.43 * sin(r2d(19.70 + theta)), str(round(sigma_x_prime,2)) + " MPa")


    if sigma_y_prime > 0:
        ax2.annotate("", xytext=(-6.04 * cos(r2d(114.44 - theta)), 6.04 * sin(r2d(114.44 - theta))), xy=(-8.86 * cos(r2d(106.39 - theta)), 8.86 * sin(r2d(106.39 - theta ))), arrowprops= dict(arrowstyle= "->"))
        ax2.text(-1-8.86 * cos(r2d(106.39 - theta)), 8.86 * sin(r2d(106.39 - theta )), str(round(sigma_y_prime,2)) + " MPa")
    else:
        ax2.annotate("", xy=(-6.04 * cos(r2d(114.44 - theta)), 6.04 * sin(r2d(114.44 - theta))), xytext=(-8.86 * cos(r2d(106.39 - theta)), 8.86 * sin(r2d(106.39 - theta ))), arrowprops= dict(arrowstyle= "->"))
        ax2.text(-0.5-8.86 * cos(r2d(106.39 - theta)), 8.86 * sin(r2d(106.39 - theta )), str(round(sigma_y_prime,2)) + " MPa")

    
    if tau_x_prime_y_prime > 0:
        ax2.annotate("", xytext=(5.59 * cos(r2d(10.30 + theta)), 5.59 * sin(r2d(10.30 + theta))), xy=(6.80 * cos(r2d(36.03 + theta)), 6.80 * sin(r2d(36.03 + theta ))), arrowprops= dict(arrowstyle= "->"))
        ax2.annotate("", xytext=(-5.59 * cos(r2d(100.30 - theta)), 5.59 * sin(r2d(100.30 - theta))), xy=(-6.80 * cos(r2d(126.03 - theta)), 6.80 * sin(r2d(126.03 - theta ))), arrowprops= dict(arrowstyle= "->"))
        ax2.text(6.80 * cos(r2d(36.03 + theta)), 6.80 * sin(r2d(36.03 + theta )), str(round(tau_x_prime_y_prime,2)) + " MPa")
    else:
        ax2.annotate("", xy=(5.59 * cos(r2d(10.30 + theta)), 5.59 * sin(r2d(10.30 + theta))), xytext=(6.80 * cos(r2d(36.03 + theta)), 6.80 * sin(r2d(36.03 + theta ))), arrowprops= dict(arrowstyle= "->"))
        ax2.annotate("", xy=(-5.59 * cos(r2d(100.30 - theta)), 5.59 * sin(r2d(100.30 - theta))), xytext=(-6.80 * cos(r2d(126.03 - theta)), 6.80 * sin(r2d(126.03 - theta ))), arrowprops= dict(arrowstyle= "->"))
        ax2.text(6.80 * cos(r2d(36.03 + theta)), 6.80 * sin(r2d(36.03 + theta )), str(round(tau_x_prime_y_prime,2)) + " MPa")
    plt.text(0, -12,"rotating angle " +str(theta), va = 'bottom', ha = 'right', fontsize = 12)
    plt.show()
    return 


plane_transformation_rotating(theta)
