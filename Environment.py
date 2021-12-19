import csv
import pandas as pd

import numpy as np
from numpy.lib.type_check import isreal

from Ball import Ball

class Environment():

    def __init__(self, boundaries) -> None:
        self.time = 0.0
        self.boundaries = boundaries
        self.balls = []
        self.keyframes = []
        self.next_event = None

    def add_ball(self) -> None:
        valid_ball_FLAG = False
        while not valid_ball_FLAG:
            # Construct a ball.
            b = Ball(
                self,                           # environment
                'blue',                         # color
                np.random.uniform(1, 10),       # mass
                np.random.uniform(1, self.boundaries[1], 3),     # center
                np.random.uniform(1, 5),        # radius
                np.random.uniform(1, self.boundaries[1], 3)/25.0      # velocity
                )

            # Add the constructed ball to this environment if it does not overlap with any ball in this environment.
            if not self.overlap_with(b) and not self.outside_boundaries(b):
                self.balls.append(b)
                valid_ball_FLAG = True

        return

    def run(self, t_end=100) -> None:
        # Initialize the run.
        print('start run - time = 0.0 sec')
        for b in self.balls:
            self.keyframes.append((self.time, *b.get_keyframe()))

        # Run.
        while (self.time < t_end):
            self.proceed_to_next_event()

        # Terminate the run.
        print(f'end run - time = {t_end:.1f} sec')
        for b in self.balls:
            self.keyframes.append((self.time, *b.get_keyframe()))

        # Output.
        print('run terminated')
        print(f'{len(self.keyframes):d} keyframes registered.')
        print(self.keyframes)

    def proceed_to_next_event(self) -> None:
        self.next_collision()
        Dt, b1, b2 = self.next_event
        self.time += Dt

        # All balls move in a straight line for Dt seconds.
        for b in self.balls:
            b.move(Dt)

        print(f'we are now {self.time:.2f} sec and proceed {Dt:.2f} sec')

        # Ball 1 and ball 2 collide if b2 is a Ball.
        # Ball 1 collides with a boundary if b2 is an integer.
        b1.collide(b2)
        self.keyframes.append((self.time, *b1.get_keyframe()))
        if (isinstance(b2, Ball)):
            self.keyframes.append((self.time, *b2.get_keyframe()))

    def next_collision(self) -> tuple:
        last_event = self.next_event
        last_colliding_balls = [last_event[1], last_event[2]] if last_event else []
        next_collision_event = None
        for b1 in self.balls:
            for b2 in self.balls:
                if b1 in last_colliding_balls and b2 in last_colliding_balls:
                    # if 3 balls collide simultaneously, the we end up in an infinite loop!
                    continue
                collision_times = b1.get_collision_time_with(b2)
                if (0 < len(collision_times)) and (np.isreal(collision_times[0])):
                    collision_times = collision_times[0 <= collision_times]
                    if (0 < len(collision_times)):
                        tmp_collision_time = np.min(collision_times)
                        if not next_collision_event or (tmp_collision_time < next_collision_event[0]):
                            next_collision_event = (tmp_collision_time, b1, b2)
        # check collisions with boundaries 
        # check collisions with cuboid boundaries.
        for b1 in self.balls:
            tmp_collision_event = b1.get_collision_time_with_boundaries(self.boundaries)
            if not next_collision_event or (tmp_collision_event[0] < next_collision_event[0]):
                next_collision_event = tmp_collision_event

        self.next_event = next_collision_event

    def overlap_with(self, b):
        overlap_FLAG = False
        for b2 in self.balls:
            if b.overlap_with(b2):
                overlap_FLAG = True
                continue
        return overlap_FLAG

    def outside_boundaries(self, b):
        return False if (self.boundaries[0] < (b.c[0] - b.r)) and ((b.c[0] + b.r) < self.boundaries[1]) \
        and (self.boundaries[2] < (b.c[1] - b.r)) and ((b.c[1] + b.r) < self.boundaries[3]) \
        and (self.boundaries[4] < (b.c[2] - b.r)) and ((b.c[2] + b.r) < self.boundaries[5]) else True

    def write_csv(self, filename):
        '''
        with open(filename, 'w') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(('time', 'entity', 'id', 'visibility', 'x', 'y', 'z', 'alpha', 'beta', 'gamma', 'scalex', 'scaley', 'scalez'))
            for row in self.keyframes:
                csv_writer.writerow(row)
                '''
        pd.DataFrame(self.keyframes, columns=('Time', 'entity', 'id', 'visibility', 'x', 'y', 'z', 'alpha', 'beta', 'gamma', 'scalex', 'scaley', 'scalez')).to_csv(filename, index=False)

                