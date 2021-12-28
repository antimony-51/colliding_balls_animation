import numpy as np

class Ball():
    ID = 0

    def __init__(self, env, color, m, c, r, v) -> None:
        self.env = env
        self.color = color
        self.m = m
        self.c = c
        self.r = r
        self.v = v
        Ball.ID += 1
        self.id = Ball.ID

    def get_impulse(self):
        return self.m * self.v

    def get_kinetic_energy(self):
        return 0.5 * self.m * np.dot(self.v, self.v)

    def move(self, dt):
        self.c += dt*self.v

    def collide(self, b2):
        # change the velocity vector.
        # b2 is an instance of Ball or a boundary.
        # 1D elastic collision along the normal.
        if (isinstance(b2, Ball)):
            self.collide_with_ball(b2)
            '''
            self.n = (b2.c - self.c)
            self.n = self.n/np.linalg.norm(self.n)
            b2.n = -self.n

            self.vn = np.dot(self.v, self.n)*self.n
            b2.vn = np.dot(b2.v, b2.n)*b2.n

            self.vt = self.v - self.vn
            b2.vt = b2.v - b2.vn

            total_mass = self.m + b2.m

            v_n_self_i = np.linalg.norm(self.vn)
            v_n_b2_i = np.linalg.norm(b2.vn)
            # initial_vn = self.vn
            # self.vn = 2*b2.m/total_mass*b2.vn + (self.m - b2.m)/total_mass*self.vn
            # b2.vn = 2*self.m/total_mass*initial_vn + (b2.m - self.m)/total_mass*b2.vn
            v_n_self_f = 2*b2.m/total_mass*v_n_b2_i + (self.m - b2.m)/total_mass*v_n_self_i
            v_n_b2_f = 2*self.m/total_mass*v_n_self_i + (b2.m - self.m)/total_mass*v_n_b2_i

            self.vn = self.vn/np.linalg.norm(self.vn)*v_n_self_f
            b2.vn = b2.vn/np.linalg.norm(b2.vn)*v_n_b2_f

            self.v = self.vn + self.vt
            b2.v = b2.vn + b2.vt
            print('done')
            '''
        else:
            axis = b2//2
            self.v[axis] = -self.v[axis]

    def collide_with_ball(self, b2):
        r1 = self.r
        r2 = b2.r
        m1 = self.m
        m2 = b2.m

        x1 = self.c[0]
        y1 = self.c[1]
        z1 = self.c[2]

        vx1 = self.v[0]
        vy1 = self.v[1]
        vz1 = self.v[2]

        vx2 = b2.v[0]
        vy2 = b2.v[1]
        vz2 = b2.v[2]

        x2 = b2.c[0]
        y2 = b2.c[1]
        z2 = b2.c[2]

        # initial impulse and kinetic energy of balls b1 and b2.
        impulse_i = self.get_impulse() + b2.get_impulse()
        ke_i = self.get_kinetic_energy() + b2.get_kinetic_energy()
        '''
        print('\n initial state:')
        print('---------------')
        print(f'impulse along x: {m1*vx1 + m2*vx2:.3f}')
        print(f'impulse along y: {m1*vy1 + m2*vy2:.3f}')
        print(f'impulse along z: {m1*vz1 + m2*vz2:.3f}')
        ke1i = 0.5*m1*(vx1**2+vy1**2+vz1**2)
        ke2i = 0.5*m1*(vx2**2+vy2**2+vz2**2)
        print(f'kinetic energy: {ke1i+ke2i:.3f}')
        print('')
        '''

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
            
        # return without modifying the velocity vectors if relative speed = 0 ****
        if (v==0):
            return
            # raise ValueError('relative velocity 0')
            

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

        '''
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
        # adjust the velocity vectors of balls b1 and b2.
        self.v[0], self.v[1], self.v[2] = vx1, vy1, vz1
        b2.v[0], b2.v[1], b2.v[2] = vx2, vy2, vz2

        # final impulse and kinetic energy of balls b1 and b2.
        impulse_f = self.get_impulse() + b2.get_impulse()
        ke_f = self.get_kinetic_energy() + b2.get_kinetic_energy()

        # validation and checksums
        if (1e-13 < np.linalg.norm(impulse_f-impulse_i)):
            raise ValueError('impulse not preserved for elastic collision of 2 balls.')
        elif (1e-13 < abs(ke_f-ke_i)):
            raise ValueError('kinetic energy not preserved for elastic collision of 2 balls.')
        elif (1e-13 < abs(t)):
            print('WARNING: collision time instant not accurately predicted.')

        return

    def get_collision_time_with(self, b2) -> np.array:
        # b2 can be an instance of Ball or an array of 4 elements that represents walls.
        # solve (r1 + r2)**2 = (c1 + t*v1)**2 + (c2 + t*v2)**2 for t
        # select the smallest t.
        # None if no collision
        Dc = self.c - b2.c
        Dv = self.v - b2.v
        Sr = self.r + b2.r

        a = np.square(np.linalg.norm(Dv))
        b = 2*np.dot(Dc, Dv)
        c = np.square(np.linalg.norm(Dc)) - np.square(Sr)

        # print(np.roots([a, b, c]))
        return np.roots([a, b, c])

    def get_collision_time_with_boundaries(self, boundaries):
        # exclude boundary of last event?
        collision_times = np.zeros(6)
        collision_times[0] = (boundaries[0] + self.r - self.c[0])/self.v[0]
        collision_times[1] = (boundaries[1] - self.r - self.c[0])/self.v[0]
        collision_times[2] = (boundaries[2] + self.r - self.c[1])/self.v[1]
        collision_times[3] = (boundaries[3] - self.r - self.c[1])/self.v[1]
        collision_times[4] = (boundaries[4] + self.r - self.c[2])/self.v[2]
        collision_times[5] = (boundaries[5] - self.r - self.c[2])/self.v[2]
        collision_times = np.where(collision_times > 1e-12, collision_times, np.inf) # 2 events; sufficiently spaced in time.
        return (collision_times.min(), self, collision_times.argmin())

    def overlap_with(self, b2) -> bool:
        d = np.linalg.norm(self.c - b2.c)
        return d <= (self.r + b2.r)

    def get_keyframe(self):
        return ('Ball', self.id, True, self.c[0], self.c[1], self.c[2], 0, 0, 0, self.r, self.r, self.r)

    def __str__(self):
        return f'   * Ball {self.id:d} with mass {self.m:.1f} kg and radius {self.r:.1f} m at position {self.c} m with velocity vector {self.v} m/s'

if __name__=='__main__':
    # 2 equal balls in opposite directions.
    b1 = Ball('None', 'blue', 1, np.array([-5,0,0]), 1, np.array([+1,0,0]))
    b2 = Ball('None', 'blue', 1, np.array([+5,0,0]), 1, np.array([-1,0,0]))
    print(b1.get_collision_time_with(b2))

    # 1 ball from the left, and one from the top, cover equal distances at equal velocities.
    b3 = Ball('None', 'blue', 1, np.array([-5,0,0]), 1, np.array([+1,0,0]))
    b4 = Ball('None', 'blue', 1, np.array([0,+5,0]), 1, np.array([0,-1,0]))
    print(b3.get_collision_time_with(b4))
