import numpy as np

class Ball():
    def __init__(self, env, color, m, c, r, v) -> None:
        self.env = env
        self.color = color
        self.m = m
        self.c = c
        self.r = r
        self.v = v

    def move(self, dt):
        self.c += dt*self.v

    def collide(self, b2):
        # change the velocity vector and the mass
        # b2 is an instance of Ball.
        # 1D elastic collision along the normal
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

    def get_collision_time_with(self, b2) -> np.array:
        # b2 can be an instance of Ball or an array of 4 elements that represents walls.
        # solve r1 + r2 = (c1 + t*v1)**2 + (c2 + t*v2)**2 for t
        # select the smallest t.
        # None if no collision
        Dc = self.c - b2.c
        Dv = self.v - b2.v
        Dr = self.r - b2.r

        a = np.square(np.linalg.norm(Dv))
        b = 2*np.dot(Dc, Dv)
        c = np.square(np.linalg.norm(Dc)) - np.square(Dr)

        # print(np.roots([a, b, c]))
        return np.roots([a, b, c])

    def get_collision_time_with_boundaries(self, boundaries):
        collision_times = np.zeros(6)
        collision_times[0] = (boundaries[0] - self.c[0])/self.v[0]
        collision_times[1] = (boundaries[1] - self.c[0])/self.v[0]
        collision_times[2] = (boundaries[2] - self.c[1])/self.v[1]
        collision_times[3] = (boundaries[3] - self.c[1])/self.v[1]
        collision_times[5] = (boundaries[5] - self.c[2])/self.v[2]
        collision_times[6] = (boundaries[6] - self.c[2])/self.v[2]
        return np.where(collision_times > 0, collision_times, np.inf).min()

    def overlap_with(self, b2) -> bool:
        d = np.linalg.norm(self.c - b2.c)
        return d <= (self.r + b2.r)

    def __str__(self):
        return f'Ball with mass {self.m:.1f} and radius {self.r:.1f} at position {self.c} with velocity vector {self.v}'
