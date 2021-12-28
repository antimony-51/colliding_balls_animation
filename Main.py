import numpy as np
np.random.seed(42)

from Environment import Environment

B = 25

boundaries = [-B, B, -B, B, -B, B] # -x, +x, -y, +y, -z, +z
Menv = Environment(boundaries)

ball_count = 40

for i in range(ball_count):
    Menv.add_ball()

Menv.run(t_end=60.0)
Menv.write_csv('keyframes.csv')