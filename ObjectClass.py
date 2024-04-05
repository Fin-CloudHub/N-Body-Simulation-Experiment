#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 10:49:05 2024

@author: finlaysime

"""

import math as m
import numpy as np

class Object():
        
    """
    Class to calculate forces, erengry, etc for simulation.
    
    """
    
    def __init__(self, name, mass, orbit, colour, G, dt, RY):
        """
        Magic method, creates an object.  Assigns all variables ready to be used in the code.

        """
        
        self.name = str(name)
        self.mass = float(mass)
        self.orbit = float(orbit)
        self.colour = str(colour)
        self.G = G
        self.dt = dt
        self.year = [0.0]
        self.t = 0.0
        self.angle = 0.0
        self.angleQ = 'TR'
        self.RY = RY
        
        
    def Begin(self, p):
        """
        A method to create the inital vectors.  Defines the position & velocity vectors as
        numpy arrays based of positions.  Generates the accelerations with a nunmpy array
        or the NewAcceleration() function.

        """
                
        self.t = 0.0
        #Position vectors
        self.r = np.array([self.orbit, 0])
        self.rE = np.array([self.orbit, 0])
        p.r = np.array([p.orbit, 0])
        p.rE = np.array([p.orbit, 0])
            
        if self.orbit == 0: #Velocity vectors
            self.v1 = np.array([0, 0])
            self.v2 = np.array([0, 0])
        else:
            speed = m.sqrt(self.G*p.mass/self.orbit)
            self.v1 = np.array([0, speed])
            self.v2 = np.array([0, speed])
            
        if self.orbit == 0: #Acceleration vectors
            self.a = np.array([0, 0])
        else:
            self.a = self.NewAcceleration(p)
        self.old_a = self.a
            
                                          
    def NewPosition1(self):
        """
        Updates the position through the Beeman method.
        
        """
        
        self.old_r = self.r
                    
        self.r = self.r+self.v1*self.dt+(4*self.a-self.old_a)*self.dt*self.dt/6

    
    def NewVelocity1(self, p):
        """
        Calculates the new velocity through the Beeman method.

        """
                
        new_a = self.NewAcceleration(p)
        self.v1 = self.v1+(2*new_a+5*self.a-self.old_a)*self.dt/6
        self.old_a = self.a
        self.a = new_a
            
            
    def NewAcceleration(self, p):
        """
        Calculates the acceleration of the next time step.

        """
        
        pos = self.r - p.r #Position of planets relative to each other
        a = -self.G*p.mass*pos/m.pow(np.linalg.norm(pos),3) #Calculates acceleration
        return a
    
    
    def CheckYear(self, p):
        """
        Checks if it is a new year yet by comparing the old position
        and the current position against x-components of the sun.
        
        """
             
        if self.old_r[1] < p.old_r[1] and self.r[1] >= p.r[1]: #Checks if the planet has crossed the
            self.t += 1.0                            #x-componenent of the sun.
            self.year.append(self.t)                          
            if self.name == 'Earth':
                print(f"Current year in Earth years: {self.year[-1]}") #Current year in Earth years
            return True
            
        else:
          return False
      
            
    def NewPosition2(self):
        """
        Calculates the new position of the object by the Direct Euler method.
         
        """
        
        self.rE = self.rE+self.v2*self.dt
        
    def NewVelocity2(self, p):
        """
        Calculates the new velocity of the object by the Direct Euler method.

        """
             
        self.a = self.NewAcceleration(p)
        
        self.v2 = self.v2+self.a*self.dt
            
    def KineticEnergyBeeman(self):
        """
        Calculates the total kinetic energy of the system with the Beeman method.
        
        """
        
        KE1 = 0.5*self.mass*(np.dot(self.v1, self.v1))
        
        return KE1
    
        
    def KineticEnergyEuler(self):
        """
        Calculates teh total kinetic energy of the system with the direct euler method.
        
        """
        
        KE2 = 0.5*self.mass*(np.dot(self.v2, self.v2))
        
        return KE2
    
        
def main():
    
     S = Object()
     S.Begin()
     S.GravForce()
     
#main()
     
     
     
     
     
     
