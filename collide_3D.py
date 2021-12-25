import numpy as np

r1 = 1
r2 = 1
m1 = 1
m2 = 1

x1 = 0
y1 = 0
z1 = 0

vx1 = 0
vy1 = 0
vz1 = 1

vx2 = 0
vy2 = 0
vz2 = 0

x2 = 0
y2 = 0
z2 = -5

# initial impulse and kinetic energy
print('\n initial state:')
print('---------------')
print(f'impulse along x: {m1*vx1 + m2*vx2:.3f}')
print(f'impulse along y: {m1*vy1 + m2*vy2:.3f}')
print(f'impulse along z: {m1*vz1 + m2*vz2:.3f}')
ke1i = 0.5*m1*(vx1**2+vy1**2+vz1**2)
ke2i = 0.5*m1*(vx2**2+vy2**2+vz2**2)
print(f'kinetic energy: {ke1i+ke2i:.3f}')
print('')

# initialize some variables ****
pi=np.pi
error=0
r12=r1+r2
m21=m2/m1
x21=x2-x1
y21=y2-y1
z21=z2-z1
vx21=vx2-vx1
vy21=vy2-vy1
vz21=vz2-vz1
       
#       vx_cm = (m1*vx1+m2*vx2)/(m1+m2) ;
#       vy_cm = (m1*vy1+m2*vy2)/(m1+m2) ;
#       vz_cm = (m1*vz1+m2*vz2)/(m1+m2) ;  

	   
# calculate relative distance and relative speed ***
d=np.sqrt(x21*x21 +y21*y21 +z21*z21)
v=np.sqrt(vx21*vx21 +vy21*vy21 +vz21*vz21)
       
# return if distance between balls smaller than sum of radii ****
if (d<r12):
    pass
    # raise ValueError('distance too small')
       
# return if relative speed = 0 ****
if (v==0):
    raise ValueError('relative velocity 0')
       

# shift coordinate system so that ball 1 is at the origin ***
x2=x21
y2=y21
z2=z21
       
# boost coordinate system so that ball 2 is resting ***
vx1=-vx21
vy1=-vy21
vz1=-vz21

# find the polar coordinates of the location of ball 2 ***
theta2=np.arccos(z2/d)
if (x2==0 and y2==0):
    phi2 = 0
else: 
    phi2=np.arctan2(y2,x2)
st=np.sin(theta2)
ct=np.cos(theta2)
sp=np.sin(phi2)
cp=np.cos(phi2)


# express the velocity vector of ball 1 in a rotated coordinate        system where ball 2 lies on the z-axis ******
vx1r=ct*cp*vx1+ct*sp*vy1-st*vz1
vy1r=cp*vy1-sp*vx1
vz1r=st*cp*vx1+st*sp*vy1+ct*vz1
fvz1r = vz1r/v
if (fvz1r>1):
    fvz1r=1
elif (fvz1r<-1): 
    fvz1r=-1 
thetav=np.arccos(fvz1r)
if (vx1r==0 and vy1r==0):
    phiv=0
else:
    phiv=np.arctan2(vy1r,vx1r)

        						
# calculate the normalized impact parameter ***
dr=d*np.sin(thetav)/r12

'''
//     **** return old positions and velocities if balls do not collide ***
       if (thetav>pi/2 || fabs(dr)>1) {
           x2=x2+x1;
           y2=y2+y1;
           z2=z2+z1;
           vx1=vx1+vx2;
           vy1=vy1+vy2;
           vz1=vz1+vz2;
           error=1;
           return;
        }
'''   

# calculate impact angles if balls do collide ***
alpha=np.arcsin(-dr)
beta=phiv
sbeta=np.sin(beta)
cbeta=np.cos(beta)
        
       
# calculate time to collision ***
t=(d*np.cos(thetav) -r12*np.sqrt(1-dr*dr))/v

'''
//     **** update positions and reverse the coordinate shift ***
       x2=x2+vx2*t +x1;
       y2=y2+vy2*t +y1;
       z2=z2+vz2*t +z1;
       x1=(vx1+vx2)*t +x1;
       y1=(vy1+vy2)*t +y1;
       z1=(vz1+vz2)*t +z1;
'''
 
       
#  update velocities ***

a=np.tan(thetav+alpha)

dvz2=2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*(1+m21))

vz2r=dvz2
vx2r=a*cbeta*dvz2
vy2r=a*sbeta*dvz2
vz1r=vz1r-m21*vz2r
vx1r=vx1r-m21*vx2r
vy1r=vy1r-m21*vy2r

       
# rotate the velocity vectors back and add the initial velocity      vector of ball 2 to retrieve the original coordinate system ****
                     
vx1=ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2
vy1=ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2
vz1=ct*vz1r-st*vx1r               +vz2
vx2=ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2
vy2=ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2
vz2=ct*vz2r-st*vx2r               +vz2

print('time of collision: ',t)
print('v1f', vx1, vy1, vz1)
print('v2f', vx2, vy2, vz2)

# final impulse and kinetic energy
print('\n final state:')
print('---------------')
print(f'impulse along x: {m1*vx1 + m2*vx2:.3f}')
print(f'impulse along y: {m1*vy1 + m2*vy2:.3f}')
print(f'impulse along z: {m1*vz1 + m2*vz2:.3f}')
ke1i = 0.5*m1*(vx1**2+vy1**2+vz1**2)
ke2i = 0.5*m1*(vx2**2+vy2**2+vz2**2)
print(f'kinetic energy: {ke1i+ke2i:.3f}')
print('')
        
'''
//     ***  velocity correction for inelastic collisions ***

       vx1=(vx1-vx_cm)*R + vx_cm;
       vy1=(vy1-vy_cm)*R + vy_cm;
       vz1=(vz1-vz_cm)*R + vz_cm;
       vx2=(vx2-vx_cm)*R + vx_cm;
       vy2=(vy2-vy_cm)*R + vy_cm;
       vz2=(vz2-vz_cm)*R + vz_cm;  

       return;
}
'''

'''
Fortran Code
(Back To) Description
'''