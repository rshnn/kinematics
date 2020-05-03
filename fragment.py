# -*- coding: utf-8 -*-
import pandas as pd 
from astropy import units
from math import sqrt
from measures.api import Angle, Mass, Measure, Speed 
from kinematics.utils import Point, Velocity, State, CartesianFrame, BaseFrame 
from kinematics.three_dof import Traj3DOF


class Fragment:    
    """TODO fill out this docstring with a description of Fragmnet 

        Attributes:
            diameter (measures.Length)

            mass (measures.Mass)

            inital_state <kinematics.utils.State>
                Contains all initial conditions (angles and velocity)  
                i.e.  velocity magnitude, az/el angles will be available through
                    this object 

            trajectory (list(kinematics.utils.state))
                List of states that the fragment goes through 
                Populated via 3dof.  Store position at each timestep 
                Might make for fun plotting later.
                We can set up for large computation of all trajectories of a Munition
                 and then plot all the trajectories using vtk 

            polar_zone (detonation.PolarZone)
                reference back to the polarzone that this fragment's mass and 
                fire_solution were sampled from 

            3dof (kinematics.3dof)
                Might be helpful to include a boolean in 3dof that states 
                whether the 3dof was already executed or not 
                Check if trajectory is empty?  Does .move return?
    """

    def __init__(self, posX=0, posY=0, posZ=0, initVelocity=0, azimuth=0, 
                 elevation=0, mass=0, presentedArea=0, dragFile='',
                 debugMode=False):

        # Save local copies of inputs
        # TODO store all the initial state information using the State class. 
        self.init_x = posX
        self.init_y = posY
        self.init_z = posZ
        self.initial_velocity = initVelocity
        self.init_azimuth = azimuth
        self.init_elevation = elevation

        self.mass = mass

        # Can this belong to the threeDOF object? 
        self.drag_file = dragFile
        self.threeDOF = Traj3DOF(debugMode)

        # Read in the drag file
        self.threeDOF.readFullFile(self.drag_file)

        # Set: mass properties, sim wind, sim wind errors, launch errors, MET errors
        self.threeDOF.setMassProperties(mass=Mass(self.mass, units.kg), 
                                        diameter=Measure(self.diameter, units.m))

        # Set the launch conditions for the fragment
        self.threeDOF.setLaunchConditions(v=Speed(self.initial_velocity, units.m / units.s),
                                          azi=Measure(self.init_azimuth, units.radian), 
                                          elv=Measure(self.init_elevation, units.radian))

        self.threeDOF.setSimWind(windRange=Speed(0, units.m / units.s),
                                 windCross=Speed(0, units.m / units.s),
                                 windVertical=Speed(0, units.m / units.s))

        self.threeDOF.setLaunchErrors(initAzi_1StdDev=Measure(0, units.radian),
                                      initElv_1StdDev=Measure(0, units.radian),
                                      initVel_1StdDev=Speed(0, units.m / units.s))

        self.threeDOF.setMETErrors(seaLvlTemp_perctErr=0 * units.dimensionless_unscaled,
                                   airDensity_perctErr=0 * units.dimensionless_unscaled,
                                   windCross_1StdDev=Speed(0, units.m / units.s),
                                   windRange_1StdDev=Speed(0, units.m / units.s),
                                   windVert_1StdDev=Speed(0, units.m / units.s))

        # Bookeeping for the fragment state along the trajectory
        self.colNames = ['t',
                         'x',
                         'y',
                         'z',
                         'vx',
                         'vy',
                         'vz',
                         'azi',
                         'elv']

        # Create a blank data frame to hold a decompresed version of self.polar_zones
        self.trajDF = pd.DataFrame(columns=self.colNames)



    def run_3dof(self, dt=0.001, lowerKineticLimit=100, lowerVelLimit=0):        
        # Fire the munition, given our initial conditions and error/MET data
        self.threeDOF.fire(initVel=Speed(self.initial_velocity, units.m / units.s), 
                           posX=Measure(self.init_x, units.m), 
                           posY=Measure(self.init_y, units.m), 
                           posZ=Measure(self.init_z, units.m), 
                           oriAzi=Measure(self.init_azimuth, units.radian), 
                           oriElv=Measure(self.init_elevation, units.radian),
                           lowerKineticLimit=Measure(lowerKineticLimit, units.J),
                           lowerVelLimit=Speed(lowerVelLimit, units.m / units.s), 
                           dt=Measure(dt, units.s))

        # Initialize the trajectory bookeeping with the initial launch conditions
        self.update_track(0,
                          self.init_x,
                          self.init_y,
                          self.init_z,
                          self.threeDOF.prev_velX,
                          self.threeDOF.prev_velY,
                          self.threeDOF.prev_velZ,
                          self.init_azimuth,
                          self.init_elevation)

        # Add the current fragment state to our track record
        self.update_track(self.threeDOF.simTime,
                          self.threeDOF.curr_posX,
                          self.threeDOF.curr_posY,
                          self.threeDOF.curr_posZ,
                          self.threeDOF.curr_velX,
                          self.threeDOF.curr_velY,
                          self.threeDOF.curr_velZ,
                          self.threeDOF.curr_azi,
                          self.threeDOF.curr_elv)

        # Run the trajectory to the ground or to the lower velocity limit
        while self.threeDOF.move():

            # Add the current fragment state to our track record
            self.update_track(self.threeDOF.simTime,
                              self.threeDOF.curr_posX,
                              self.threeDOF.curr_posY,
                              self.threeDOF.curr_posZ,
                              self.threeDOF.curr_velX,
                              self.threeDOF.curr_velY,
                              self.threeDOF.curr_velZ,
                              self.threeDOF.curr_azi,
                              self.threeDOF.curr_elv)

    def update_track(self, t, x, y, z, vx, vy, vz, azi, elv):
        # Build the initial row values
        vals = {'t': t,
                'x': x,
                'y': y,
                'z': z,
                'vx': vx,
                'vy': vy,
                'vz': vz,
                'azi': azi,
                'elv': elv}
        # Create a temporary dataframe representing a row to be concatenated
        tmp = pd.DataFrame(data=vals, columns=self.colNames, index=[0])

        # Append the single row dataframe into the uncompressed container
        self.trajDF = pd.concat([self.trajDF, tmp], ignore_index=True)

