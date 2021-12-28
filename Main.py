import numpy as np
np.random.seed(42)

from Environment import Environment

B = 25

Menv = Environment([-B, B, -B, B, -B, B]) # -x, +x, -y, +y, -z, +z
Menv.add_n_balls(40)
Menv.run(t_end=60.0)
Menv.write_csv('keyframes.csv')