# -*- coding: utf-8 -*-

import numpy as np 


class Vector3(np.ndarray):
    """Defines a Vector of shape (3,) that inherits from numpy.ndarray 
    """

    def __new__(cls, *args):
        if len(args) == 1: 
            if isinstance(args[0], Vector3):
                return args[0].copy()

            if isinstance(args[0], np.matrix):
                return Vector3(args[0].flatten().tolist()[0])
        data = _massage_args(args)
        array = np.array(data, dtype=np.float, copy=True)
        return np.ndarray.__new__(cls, shape=(3,), buffer=array)


    def __repr__(self):
        return '{0}{1}'.format(self.__class__.__name__, repr(tuple(self)))


    def magnitude(self):
        return np.linalg.norm(self)



class UnitVector3(Vector3):
    """Defines a Vector3 of shape (3,) that is a Vector3 of magnitude 1.
    """

    def __new__(cls, *args):
        vec = super().__new__(cls, *args)
        if len(args) is 0:
            return vec
        return vec / vec.magnitude()



def _massage_args(args):
    """Handles various argument inputs to Vector3 
    Returns a tuple/list that numpy.array can parse appropriately  
    """

    if len(args) == 0:
        data = 3 * (0, )

    elif len(args) == 1:
        data = args[0]
        if len(data) != 3: 
            raise TypeError('Vector3 expects a sequence with 3 elements and 1 '
                            'argument was given with {} elements.'
                            .format(len(data)))
    elif len(args) == 3: 
        data = args

    else: 
        raise TypeError('Vector3 expects a sequence with 3 elements.')


    assert len(data) == 3

    try:
        return tuple(map(float, data))

    except (TypeError, ValueError): 
        raise TypeError('Vector3 generation could not convert input elements '
                        'to floats.')
