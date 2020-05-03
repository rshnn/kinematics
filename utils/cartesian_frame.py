
# -*- coding: utf-8 -*-

import quaternion
from astropy import units 
from .vector3 import Vector3
from .point import Point
from measures.api import Angle 
from .coordinate_frame import CoordinateFrame, BaseFrame 


class CartesianFrame(CoordinateFrame):
    """x, y, z

    origin (Point) 
    orientation (tuple/list(Angle) : represent euler angles 
    base_frame (BaseFrame)
    """

    def __init__(self, base_frame=BaseFrame(), 
                 translation=Vector3(),
                 orientation=3 * (Angle(0, units.rad), ),
                 dimension_unit=units.m,
                 axes_convention='', 
                 name='my_frame'):
        """
        axes_convention: Select a common frame transformation 
            ENU - east, north, up 
            NED - north, east, down 

        orientation (list/tuple(Angle) x3: 
        Euler angles that describe the rotation of the coordinate frame with 
        respect to the base_frame.  That is (alpha, beta, gamma). 
        The rotation is given by an initial rotation through gamma about the 
        z axis, followed by a rotation through beta about the new y axis, 
        followed by a final rotation through alpha about the new z axis.
        """

        # Handle axes_convention inputs 
        if axes_convention and isinstance(axes_convention, str):
            if axes_convention.upper() not in COMMON_FRAMES: 
                raise ValueError('Input axes_convention did not match internal ' 
                                 'dictionary of common axes conventions.')
            translation, orientation = COMMON_FRAMES[axes_convention.upper()]
            if name is 'my_frame': 
                name = axes_convention


        super().__init__(orientation=orientation, name=name)

        self.origin = Point(translation, base_frame, dimension_unit)
        if not isinstance(base_frame, BaseFrame):
            raise TypeError('Expected input base_frame to be of type BaseFrame.')
        self.base_frame = base_frame


    def __repr__(self):
        return "{0}({1}, origin: {2}, orientation: {3})".format(self.__class__.__name__, 
                                                                self.name, 
                                                                self.origin, 
                                                                quaternion.as_euler_angles(self.orientation))




COMMON_FRAMES = {
    'ENU': [Vector3(), 
            (Angle(0, units.deg), Angle(0, units.deg), Angle(0, units.deg))], 
    'NED': [Vector3(), 
            (Angle(0, units.deg), Angle(180, units.deg), Angle(90, units.deg))], 

    # EUN: changes the order of how to describe the point 
    # 'PRISIM': [Vector3(), 
    #         (not a simpple rotation...switches points y and z)],
}
