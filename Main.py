import numpy as np
np.random.seed(42)

from Environment import Environment
from Ball import Ball

B = 10

boundaries = [-B, B, -B, B, -B, B]
Menv = Environment(boundaries)

ball_count = 5

for i in range(ball_count):
    Menv.add_ball()

for b in Menv.balls:
    print(str(b))

Menv.run()