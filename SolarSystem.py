#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:06:38 2024

@author: finlaysime

s2212677

"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from ObjectClass import Object
import numpy as np
import math as m

class SolarSystem(Object):
    
    #Class variable used to determine which experiment is going to be run over the simulation.
    ExpNo = 3
    
    def __init__(self):
        """
        Magic Method, opens a file and processes the data to then create objects of the
        Object class.  First removes anything in file starting with a hashtag then appends
        a list with the information.  The information is stored line by line so is split
        at ',' then with each paramter now an individual item in a list, using the object class
        the objects of the solar system are created in a list.

        """
        
        #Open and reads the file
        self.data = [line.rstrip() for line in open('Test.txt')] #Appends a list with each line.
        open('Test.txt').close()
        
        #Removes any strings that begin with a #
        for i in range(0, len(self.data)-1):
            if '#' in self.data[i]:
                del self.data[i]
                
        
        #Creates a list of lists for each object to be scripted
        self.data_array = []
        for i in range(0, len(self.data)):
            self.data_array.append(self.data[i].split(', ')) #Split each line in the list by ','
            
        self.objects = []
        
        
        for j in range(0, len(self.data_array)): #Creates the objects of the solar system using the object class
            x = Object(str(self.data_array[j][0]), float(self.data_array[j][1]),
                       float(self.data_array[j][2]), str(self.data_array[j][3]),
                       float(self.data_array[j][4]), float(self.data_array[j][5]),
                       float(self.data_array[j][6]))
            self.objects.append(x)
            

        #Simulation parameters
        self.max = 0
        self.niter = 100000
        
        #Initaliser functions
        for i in range(0, len(self.objects)):
            self.objects[i].Begin(self.objects[0])
        
        #Work out the lagest orbit
        for i in range(0, len(self.objects)):
            if self.objects[i].orbit > self.max:
                self.max = self.objects[i].orbit #Determine the maximum size of the simulation
                
        self.y_data_1_E = []
        self.y_data_2_E = []
        self.time_years = [0]
        
                
    def __str__(self):
        """
        Prints the parameters of the simulation along with each object of the simulation
        and their unique parameters.
        
        """
        
        print("Simulation parameters:\n") #Parameters G and dt stated.
        print(f"Gravitaional constant: {self.objects[0].G}Nm^2/kg^2")
        print(f"Time-step (dt): {self.objects[0].dt}s\n")
        
        print("Here is a list of the objects and their paramters.\n")
        
        #Print the name of each term with it's value
        for j in range(0, len(self.objects)):
            print(f"Name: {self.objects[j].name}") #Name
            print(f"Mass: {self.objects[j].mass}kg") #Mass
            print(f"Orbital Radius: {self.objects[j].orbit}m") #Orbit
            print(f"Colour: {self.objects[j].colour}\n") #Colour code
             
        return ""
    
    @classmethod
    def ExpNoFunc(cls):
        """
        Choose if code runs experiment 2: Checking for energy conservation
        with Beeman method and Euler method, or experiment 4: Finding out
        when the planets align with each other.
        
        """
                
        while True:
            print("This experiment can run one of two experiments at a time.")
            print("For experiment 2 (energy conservation of Beeman and Eueler), enter 2")
            print("For experiment 4 (planetry alignment), enter 4")
            c = input("") #Input value to determine chosen experiment
        
            if c == str(2):
                SolarSystem.ExpNo = 2
                break
            elif c == str(4):
                SolarSystem.ExpNo = 4
                break
            else:
                print("That was not a recognised input.")  #Loops back if input was not 2 or 4
        

        
    def energyB(self):
        """
        Calculates the total energy of the system using the Beeman method.
        First the kinetic energy is claculated in a loop, then the potential energy of every body
        is calculated in a loop inside the first loop.  The total energy is then determined and
        returned from the funciton.

        """
        
        ke = 0.0 #Kinetic energy
        pe = 0.0 #Potential energy

        for j in range(0, len(self.objects)):
            ke += self.objects[j].KineticEnergyBeeman() #Calcualte the kinetic energy using beeman method.
            for k in range(0, len(self.objects)):
                if k != j:
                    r = np.linalg.norm(self.objects[k].r - self.objects[j].r) #Distance between 2 objects
                    pe -= self.objects[j].G*self.objects[j].mass*self.objects[k].mass / r #Calculate the potential
                                                                                #energy using beeman method.
        pe = pe / 2 #Half the potential energy so each object is not counted twice (position vector flips).
        totEnergyB = ke + pe #Calcualte total energy
        return totEnergyB

    def energyE(self):
        """
        Calculates the total energy of the system using the Euler method.
        First the kinetic energy is claculated in a loop, then the potential energy of every body
        is calculated in a loop inside the first loop.  The total energy is then determined and
        returned from the funciton.
        
        """
        
        ke = 0.0 #Kinetic energy
        pe = 0.0 #Potential energy

        for j in range(0, len(self.objects)):
            ke += self.objects[j].KineticEnergyEuler() #Calcualte the kinetic energy using euler method.
            for k in range(0, len(self.objects)):
                if k != j:
                    r = np.linalg.norm(self.objects[k].rE - self.objects[j].rE) #distance between 2 objects
                    pe -= self.objects[j].G*self.objects[j].mass*self.objects[k].mass / r #Calculate the potential
                                                                                #energy using euler method.
        pe = pe / 2 #Half the potential energy so each object is not counted twice (position vector flips).
        totEnergyE = ke + pe #Calculate total energy
        return totEnergyE   
        
                
    def init(self):
        """
        Initaliser for the simulation.
        
        """
        return self.patches
    
    def Alignment(self, p):
        """
        This function checks if the planets (excluding the sun) are aligned with one another within 5 degrees.
        This is done by calculating the arctan of the absolute values of the x and y co-ordinates.  Next a 
        correction must be done, as with the absolute values the same angle could be returned even when the
        objects are in different quadrants.  The function then returns true when all are aligned.
        
        """
        
        #Varaible to determine if all are aligned.
        a = 0
        
        for i in range(1, len(self.objects)-1):   
            self.objects[i].angle = np.arctan2((self.objects[i].r[1]-self.objects[0].r[1]), self.objects[i].r[0])
            
        for j in range(1, len(self.objects)-1):
            for k in range(1, len(self.objects)-1):
                if j != k:
                    if self.objects[j].angle >= self.objects[k].angle-0.087 and self.objects[j].angle <= self.objects[k].angle+0.087:
                         a += 1 #Add one each time another planet is aligned
                         
        if a == 20: #When all the planets are aligned.
            return True #This is fullfilled when they are aligned.
    
    def animate(self, i):
        """
        Animation Function.  Is looped in the run function.  With each iteration & time-step
        the new position, velocity and hence acceleration is calculated.  Only the vectors
        calculated using the Beeman method are returned to the patches displayed on the simulation.
        The total energy is also calculated in this function using both the Beeman and Euler methods.
        If the planetary alignment experiment was chosen then the comparison of the total energys
        is skipped and instead the funciton waits for planetary alignment conditions to be met.
        
        """
        
        time = (i+1)*self.objects[0].dt #Calcualtes the time in the simulation.
        self.time_years.append(time/60/60/24/365.256) #Converts simulation time into years.
        
        for j in range(0, len(self.objects)):
            self.objects[j].NewPosition1() #Determine the position, Beeman method
            self.patches[j].center = self.objects[j].r #Updates the position of the patch using the beeman method
            
            self.objects[j].NewPosition2() #Calculates the position using the euler method
                          
        for k in range(0, len(self.objects)):
            for j in range(0, len(self.objects)):
                if k != j: #Determine the new veocity, Beeman method
                    self.objects[k].NewVelocity1(self.objects[j])
                    self.objects[k].NewVelocity2(self.objects[j]) #Calculates the veocity using the euler method
        
        for j in range(0, len(self.objects)):
            if self.objects[j].CheckYear(self.objects[0]): #Checks the year of each object.
                if self.objects[j].year[-1] < 2:
                    print(f"{self.objects[j].name}") #Prints the period of the orbit in Earth years once calcualted.
                    print(f"Orbital Period In Earth years: {self.time_years[-1]}")
                    x = self.objects[j].RY/self.time_years[-1]
                    print(f"Compared to reality: {x}")
                
                #Once an Earth year, calculates the total energy of the system using both the
                #Beeman method and Euler method.
                if self.objects[j].name == 'Earth':
                    energyB = self.energyB()
                    energyE = self.energyE()
                    
                    #Appends both values for the total energy to seperate lists to be plotted ona  graph.
                    self.y_data_1_E.append(energyB)
                    self.y_data_2_E.append(energyE)
                    f = open("E.txt", 'a') #Creates a text file to store energy information in.
                    print(f"Total energy is: {energyB: .3e}J")
                    f.write(f"Current year in Earth years: {self.objects[3].year[-1]}.\n")
                    f.write(f"Total energy from Beeman is: {energyB: .3e}J\n")
                    f.write(f"Total energy from Euler is: {energyE: .3e}J\n")
        
        #If experiment 2 was chosen, runs code to plot energy information of a graph.
        if SolarSystem.ExpNo == 2:
            for j in range(0, len(self.objects)):
                if self.objects[j].name == 'Earth':
                    if self.objects[j].year[-1] == 300.0: #Creates graph after 300 Earth years have passed.             
                        if len(self.objects[j].year) != len(self.y_data_1_E) and len(self.y_data_2_E):
                            self.objects[j].year.pop(0) #Remove first item of year incase it is 1 extra.
                        fig, ax = plt.subplots(1, 1)
                        plt.plot(self.objects[j].year, self.y_data_1_E, self.objects[j].year, self.y_data_2_E)
                        plt.xlabel("Time (years)")
                        plt.ylabel("Energy (Joules)")
                        plt.show()
                        
        #If experiment 4 was chosen, runs code to check for alignment as simulation runs.
        if SolarSystem.ExpNo == 4:
            if self.Alignment(self.objects[0]):
                print("Planetary alignment has occured!!!")
                print(f"It took the simulation {self.objects[3].year[-1]} Earth years.")
                if self.time_years[-1] > 2:
                    d = open("Alignment.txt", 'a')
                    d.write(f"It took the simulation {self.time_years[-1]} Earth years.\n")
        
        return self.patches
    
    
    def run(self):
        """
        Runs the simulation.  Creates a plot and a list of patches as circles to appear on the simulation
        and sets the limits of the simulation plot.  Calls the self.animate funciton and loops it.

        """
        
        fig = plt.figure()
        ax = plt.axes()
        
        self.patches = [] #Create a list of patches
        
        #Add a shape or object for each body in the simulation
        for i in range(0, len(self.objects)):
            self.patches.append(plt.Circle(self.objects[i].r, self.max*0.01,
                                           color = self.objects[i].colour, animated=True))
            
        for i in range(0, len(self.patches)):
            ax.add_patch(self.patches[i])
        
        ax.axis("scaled")
        ax.set_xlim(-self.max*3, self.max*3) #Set upper and lower limits of simulation
        ax.set_ylim(-self.max*3, self.max*3)
        ax.set_xlabel("x-axis") #Label the axis
        ax.set_ylabel("y-axis")
        
        #Run the animation
        self.anim = FuncAnimation(fig, self.animate, init_func=self.init, frames=self.niter,
                                  repeat=True, interval=1, blit=True)

        
        plt.show()
        

def main():
    
    SS = SolarSystem()
    print(SS)
    SolarSystem.ExpNoFunc()
    SS.run()
    
main()
        