# -*- coding: utf-8 -*-
from .point import Point 
from .velocity import Velocity  


class State():
    """Represents a configuration in state space.  State space is the space of 
        all 3D positions and associated velocity vectors.  

    position <utils.Point>
    velocity <utils.Velocity>
        velocity has an orientation that points towards the trajectory arc 

    TODO: orientation (the way the fragment is facing)
    """

    def __init__(self, position, velocity):

        if not isinstance(position, Point):
            raise TypeError('Expected input position to be of type ' 
                            'kinematics.utils.Point')

        if not isinstance(velocity, Velocity):
            raise TypeError('Expected input velocity to be of type ' 
                            'kinematics.utils.Velocity')


        if position.frame is not velocity.frame: 
            raise ValueError('Expected both position and velocity to be '
                             'described in reference to the same '
                             'CartesianFrame')

        self.position = position 
        self.velocity = velocity


    def __str__(self):
        return "{0}({1}, {2})".format(self.__class__.__name__, 
                                      self.position, self.velocity)


    def __repr__(self):
        return "{0}({1}, {2})".format(self.__class__.__name__, 
                                      self.position, self.velocity)
