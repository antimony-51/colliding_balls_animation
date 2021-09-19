import numpy as np
np.random.seed(42)

from Environment import Environment
from Ball import Ball

B = 10

boundaries = [-B, B, -B, B]
Menv = Environment(boundaries)

ball_count = 25

for i in range(ball_count):
    Menv.add_ball(
        Ball(
            Menv,                           # environment
            'blue',                         # color
            np.random.uniform(1, 10),       # mass
            np.random.uniform(1, B, 3),     # center
            np.random.uniform(1, 5),        # radius
            np.random.uniform(1, B, 3)      # velocity
            )
        )

Menv.run()