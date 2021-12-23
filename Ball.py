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

    def move(self, dt):
        self.c += dt*self.v

    def collide(self, b2):
        # change the velocity vector.
        # b2 is an instance of Ball or a boundary.
        # 1D elastic collision along the normal.
        if (isinstance(b2, Ball)):
            self.n = (b2.c - self.c)
            self.n = self.n/np.linalg.norm(self.n)
            b2.n = -self.n

            self.vn = np.dot(self.v, self.n)
            b2.vn = np.dot(b2.v, b2.n)

            self.vt = self.v - self.vn
            b2.vt = b2.v - b2.vn

            total_mass = self.m + b2.m
            self.vn = 2*b2.m/total_mass*b2.vn + (self.m - b2.m)/total_mass*self.vn
            b2.vn = 2*self.m/total_mass*self.vn + (b2.m - self.m)/total_mass*b2.vn

            self.v = self.vn + self.vt
            b2.v = b2.vn + b2.vt
        else:
            axis = b2//2
            self.v[axis] = -self.v[axis]

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
        return f'Ball {self.id:d} with mass {self.m:.1f} and radius {self.r:.1f} at position {self.c} with velocity vector {self.v}'

if __name__=='__main__':
    # 2 equal balls in opposite directions.
    b1 = Ball('None', 'blue', 1, np.array([-5,0,0]), 1, np.array([+1,0,0]))
    b2 = Ball('None', 'blue', 1, np.array([+5,0,0]), 1, np.array([-1,0,0]))
    print(b1.get_collision_time_with(b2))

    # 1 ball from the left, and one from the top, cover equal distances at equal velocities.
    b3 = Ball('None', 'blue', 1, np.array([-5,0,0]), 1, np.array([+1,0,0]))
    b4 = Ball('None', 'blue', 1, np.array([0,+5,0]), 1, np.array([0,-1,0]))
    print(b3.get_collision_time_with(b4))
