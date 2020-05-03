# -*- coding: utf-8 -*-

from .vector3 import Vector3
import quaternion
from astropy import units 
from measures.api import Angle 


class CoordinateFrame:
    """
    origin (Vector3) 
    orientation (quarternion)
    """        

    def __init__(self, origin=Vector3(), 
                 orientation=3 * (Angle(0, units.rad), ),
                 name='world'):
        self.origin = origin  # * dimension_unit

        if not all(isinstance(i, Angle) for i in orientation):
            raise TypeError('Expected input transformation to be a list/tuple '
                            'of length 3 with all elements of type Angle.')        
        euler_angles = [i.to(units.rad).value for i in orientation]
        self.orientation = quaternion.from_euler_angles(euler_angles)
        self.name = name



class BaseFrame(CoordinateFrame):
    """
    origin (Vector3) 
    orientation (quaternion?)
    """

    def __init__(self, origin=Vector3(), 
                 orientation=3 * (Angle(0, units.rad), ),
                 *args):
        super().__init__(*args)
