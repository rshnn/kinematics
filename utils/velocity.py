# -*- coding: utf-8 -*-

import math 
import numbers 
import numpy as np 
import quaternion
from astropy import units 
from measures.api import Speed, Angle, UnitsError 
from .vector3 import Vector3 
from .cartesian_frame import CoordinateFrame  


class Velocity:
    """Represents a velocity vector in 3D vector space in reference to a 
        coordinate frame.coordinate

        Attributes: 
            velocity <utils.Vector3>
            frame <utils.CartesianFrame> 
            dimension_unit <something about m/s or speed equivalent>

    """

    def __init__(self, vel, frame, dimension_unit=units.m / units.s):

        if isinstance(vel, (list, tuple, np.ndarray)):
            if all(isinstance(i, numbers.Real) for i in vel):
                self.vector = Vector3(vel)
            else: 
                raise TypeError('Expected all elements of input list vel ' 
                                'to be of type int or float')

        elif isinstance(vel, Vector3):
            self.vector = vel

        else: 
            raise TypeError('Expected input vel to be of type Vector3 or a ' 
                            'list of three elements of numeric type.')

        if isinstance(frame, CoordinateFrame):
            self._frame = frame 

        if not dimension_unit.is_equivalent(units.m / units.s):
            raise UnitsError('Expected input dimension unit to be a ' 
                             'astropy.unit combination representing speed.')  

        self.units = {0: dimension_unit, 1: dimension_unit, 2: dimension_unit}

        self._x = self.vector[0]
        self._y = self.vector[1]
        self._z = self.vector[2]


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
        self.vector[0] = value 


    @y.setter 
    def y(self, value):
        self._y = value 
        self.vector[1] = value 

    @z.setter 
    def z(self, value):
        self._z = value 
        self.vector[2] = value 



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
            vect = [sum(x).to(units.m / units.s).value 
                    for x in zip(self.to_units_array(), 
                                 other.to_units_array())]
            return Velocity(Vector3(vect), 
                            frame=self._frame, 
                            dimension_unit=units.m / units.s)
        else:
            raise ValueError('Expected {0} and {1} to refer to the same '
                             'coordinate frame instance.'
                             .format(self, other)) 


    def __sub__(self, other):
        """Returns new point taht is the sum of self and other. 
            Input puts are not mutated. 
        """
        if self._frame is other._frame:
            vect = [(x[0] - x[1]).to(units.m / units.s).value 
                    for x in zip(self.to_units_array(), 
                                 other.to_units_array())]
            return Velocity(Vector3(vect), 
                            frame=self._frame, 
                            dimension_unit=units.m / units.s)
        else:
            raise ValueError('Expected {0} and {1} to refer to the same '
                             'coordinate frame instance.'
                             .format(self, other)) 

    def magnitude(self):
        """Get magnitude of point (norm of the point).  
            Distance from the origin of the point's reference frame.  
        """ 
        return np.linalg.norm(self.vector)

    def get_azimuth_elevation(self):
        """Returns azimuth and elevation angle 

            Azimuth is the angle measurment of a vector from the x=0 plane 
            Elevation is the angle measurement of a vector along the z=0 plane

            Azimuth and elevation angles coorespond to phi and pi/2-theta 
            respectively from the standard sphereical coordinate system 
            (r, theta, phi).  
            https://www.mathworks.com/help/phased/ug/spherical-coordinates.html
        """

        # Convert all coords to meters/second (SI units)
        x, y, z = self.to_units_array()
        x_si, y_si, z_si = [i.to(units.m / units.s).value for i in [x, y, z]]

        # get spherical coords 
        r = (self.magnitude() * self.units[0]).to(units.m / units.s).value
        theta = math.acos(z_si / r)
        phi = math.atan2(y_si, x_si)

        az = phi 
        el = np.pi / 2 - theta

        return (Angle(az, units.rad).to(units.deg), 
                Angle(el, units.rad).to(units.deg))



    def to_frame(self, new_frame):
        """Returns a new Point in reference to the new coordinate frame.  
            Input point is not mutated.    
        """
        if self._frame is new_frame: 
            return self

        world_coords = quaternion.rotate_vectors(self._frame.orientation, 
                                                 self.vector)
        x = quaternion.rotate_vectors(new_frame.orientation, world_coords)
        x = np.around(x, decimals=12)
        new_vel = Vector3(new_frame.origin.coords + x)

        return Velocity(vel=new_vel, frame=new_frame, 
                        dimension_unit=self.units[0])


    def to_units_array(self):
        """Returns list of coordinates with associated astropy.units  
            Default return unit is meters per second (units.m/units.s)
        """
        return [Speed(self.vector[i], self.units[i]) for i in range(3)]           
