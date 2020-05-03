# -*- coding: utf-8 -*-

import math 
import numbers 
import numpy as np
import quaternion
from astropy import units 
from measures.api import Length, Angle, UnitsError
from .vector3 import Vector3
from .coordinate_frame import CoordinateFrame


class Point:
    """Represents a point in 3D vector space in reference to a coordinate frame. 

        Attributes: 
            coords <Vector3>
            frame <CoordinateFrame>
            dimension_unit <astropy.units of length> 
    """

    def __init__(self, coords, frame, dimension_unit=units.m):

        if isinstance(coords, (list, tuple, np.ndarray)):
            if all(isinstance(i, numbers.Real) for i in coords):
                self.coords = Vector3(coords)
            else:
                raise TypeError('Expected all elements of input list coords '
                                'to be of type int or float. ')

        elif isinstance(coords, Vector3):
            self.coords = coords 
        else: 
            raise TypeError('Expected input coords to be of type Vector3 or a '
                            'list of three elements of numeric type.')


        if isinstance(frame, CoordinateFrame):
            self._frame = frame


        if not dimension_unit.is_equivalent(units.m):
            raise UnitsError('Expected input dimension_unit to be a '
                             'astroy.unit representing length.')
        # TODO make units a single astropy unit.  It will be illegal to 
        #   have mixed units in a Point instance.    
        self.units = {0: dimension_unit, 1: dimension_unit, 2: dimension_unit}

        self._x = coords[0]
        self._y = coords[1]
        self._z = coords[2]


    @property
    def frame(self):
        return self._frame


    @property
    def x(self):
        return self._x 

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z


    @x.setter 
    def x(self, value):
        self._x = value 
        self.coords[0] = value 


    @y.setter 
    def y(self, value):
        self._y = value 
        self.coords[1] = value 

    @z.setter 
    def z(self, value):
        self._z = value 
        self.coords[2] = value 




    def __str__(self):
        if isinstance(self._frame, CoordinateFrame):
            a, b, c = self.to_units_array()
            # a, b, c = [i.to(self.units[0]) for i in [a, b, c]]

        else: 
            raise TypeError('Type of coordinate frame is invalid.')

        frame_name = str(self._frame.__class__.__name__)
        short_name = frame_name.replace('Frame', '')

        return "{0}[{5}({4})]({1:.3f}, {2:.3f}, {3:.3f})".format(self.__class__.__name__,
                                                                 a, b, c, short_name, 
                                                                 self._frame.name,)



    def __repr__(self):
        if isinstance(self._frame, CoordinateFrame):
            a, b, c = self.to_units_array()
            # a, b, c = [i.to(self.units[0]) for i in [a, b, c]]            
        else: 
            raise TypeError('Type of coordinate frame is invalid.')

        frame_name = str(self._frame.__class__.__name__)
        short_name = frame_name.replace('Frame', '')
        return "{0}[{5}({4})]({1:.9f}, {2:.9f}, {3:.9f})".format(self.__class__.__name__,
                                                                 a, b, c, short_name, 
                                                                 self._frame.name,)

    def __add__(self, other):
        """Returns new point taht is the sum of self and other. 
            Input puts are not mutated. 
        """
        if self._frame is other._frame:
            coords = [sum(x).to(units.m).value 
                      for x in zip(self.to_units_array(), 
                                   other.to_units_array())]
            return Point(Vector3(coords), 
                         frame=self._frame, 
                         dimension_unit=units.m)
        else:
            raise ValueError('Expected {0} and {1} to refer to the same '
                             'coordinate frame instance.'
                             .format(self, other)) 



    def __sub__(self, other):
        """Returns new point that is the difference of self and other.
            Input points are not mutated.
        """
        if self._frame is other._frame:
            coords = [(x[0] - x[1]).to(units.m).value
                      for x in zip(self.to_units_array(), 
                                   other.to_units_array())]

            return Point(Vector3(coords), 
                         frame=self._frame, 
                         dimension_unit=units.m)
        else:
            raise ValueError('Expected {0} and {1} to refer to the same '
                             'coordinate frame instance.'
                             .format(self, other)) 




    def distance(self, other):
        """Returns distance from self to other.  
        """
        return (self - other).magnitude()



    def magnitude(self):
        """Get magnitude of point (norm of the point).  
            Distance from the origin of the point's reference frame.  
        """ 
        return np.linalg.norm(self.coords)


    def to_units_array(self):
        """Returns list of coordinates with associated astropy.units  
        """
        return [Length(self.coords[i], self.units[i]) for i in range(3)]           


    def as_spherical_coords(self):
        """Returns spherical coordinates 
        """

        # Convert all coords to meters
        x, y, z = self.to_units_array()
        x_m, y_m, z_m = [i.to(units.m).value for i in [x, y, z]]

        r = (self.magnitude() * self.units[0]).to(units.m).value 
        theta = math.acos(z_m / r)
        phi = math.atan2(y_m, x_m)

        return (Length(r, self.units[0]), 
                Angle(theta, units.rad), 
                Angle(phi, units.rad))




    def to_frame(self, new_frame):
        """Returns a new Point in reference to the new coordinate frame.  
            Input point is not mutated.    
        """
        if self._frame is new_frame: 
            return self

        world_coords = quaternion.rotate_vectors(self._frame.orientation, 
                                                 self.coords)
        x = quaternion.rotate_vectors(new_frame.orientation, world_coords)
        x = np.around(x, decimals=12)
        new_coords = Vector3(new_frame.origin.coords + x)

        return Point(coords=new_coords, frame=new_frame, 
                     dimension_unit=self.units[0])


    # TODO allow increment to have units.  (i.e. be a list of Lengths) 
    def translate(self, increment):
        """Returns a point that is translated by the input array values.  
            Assumes that the units of increment values are the same as the 
            current point's units.  
            Input point is not mutated.  
        """
        coords = self.coords + Vector3(increment)
        return Point(coords=coords, frame=self._frame, 
                     dimension_unit=self.units[0])
