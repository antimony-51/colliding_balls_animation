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