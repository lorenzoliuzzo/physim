
# author:          Lorenzo Liuzzo
# email:           lorenzoliuzzo@outlook.com
# description:     Position(class) keeps track of coordinates and velocity of an object in a 3D system.
# last updated:    24/06/2022

import numpy as np


class Position : 
    
    # class members 
    coord = np.empty(3)
    vel = np.empty(3)
    pos = np.empty(shape = (3, 2), dtype = float) 

    def __init__(self, coord, vel) :
        self.coord = coord
        self.vel = vel
        self.pos[:, 0] = self.coord
        self.pos[:, 1] = self.vel
        

    # set coordinates
    def set_coordinates(self, coord) : 
        self.coord = coord

    def set_coord_x(self, x) : 
        self.coord[0] = x   

    def set_coord_y(self, y) : 
        self.coord[1] = y   

    def set_coord_z(self, z) : 
        self.coord[2] = z   


    # set velocity
    def set_velocity(self, vel) : 
        self.vel = vel

    def set_vel_x(self, x) : 
        self.vel[0] = x   

    def set_vel_y(self, y) : 
        self.vel[1] = y   

    def set_vel_z(self, z) : 
        self.vel[2] = z   


    # set position
    def set_position(self, coord, vel) : 
        self.coord = coord
        self.vel = vel
        self.update_position()
    
    def update_position(self) : 
        self.pos[:, 0] = self.coord
        self.pos[:, 1] = self.vel

    
    # get coordinates
    def get_coordinates(self) : 
        return self.coord 

    def get_coord_x(self) : 
        return self.coord[0]    
        
    def get_coord_y(self) : 
        return self.coord[1]    
    
    def get_coord_z(self) : 
        return self.coord[2]    
    

    # get velocity
    def get_velocity(self) : 
        return self.vel 
    
    def get_vel_x(self) : 
        return self.vel[0]    
        
    def get_vel_y(self) : 
        return self.vel[1]    
    
    def get_vel_z(self) : 
        return self.vel[2]   


    # get position
    def get_position(self) -> np.array((3, 2)):
        return self.pos


    # get magnitude
    def get_coord_magnitude(self) : 
        return np.sqrt(self.get_coord_x() ** 2 +
                        self.get_coord_y() ** 2 +
                        self.get_coord_z() ** 2)

    def get_vel_magnitude(self) : 
        return np.sqrt(self.get_vel_x() ** 2 +
                        self.get_vel_y() ** 2 +
                        self.get_vel_z() ** 2)


    # get distance
    def get_distance(self, pos) : 
        return np.sqrt((pos.get_coord_x() - self.get_coord_x()) ** 2 +
                    (pos.get_coord_y() - self.get_coord_y()) ** 2 +
                    (pos.get_coord_z() - self.get_coord_z()) ** 2)


    # get polar coord
    def get_rho(self) : 
        return np.sqrt(self.get_coord_x() ** 2 +
                       self.get_coord_y() ** 2)
        
    def get_phi(self, pos2 = None) :     
        return np.arctan2(pos2.get_coord_y() - self.y, 
                        pos2.get_coord_x() - self.x)

    def get_theta(self, pos2 = None) : 
        if pos2 is None : 
            if self.get_coord_z() == 0 : 
                return 0
            else : 
                return np.arccos(self.get_coord_z() / self.get_magnitude())
        else : 
            return np.arccos((pos2.get_coord_z() - self.get_coord_z()) / self.get_distance(pos2))


    # get vel direction
    def get_phi_vel(self) :
        return np.arctan2(self.get_vel_y(), self.get_vel_x())

    def get_theta_vel(self) : 
        if self.get_vel_z() == 0 : 
            return 0
        else : 
            return np.arccos(self.get_vel_z() / self.get_magnitude())


    # print
    def print_coordinates(self) : 
        print("Coordinates : [", self.coord[0], "] [", self.coord[1], "] [", self.coord[2], "] \n")

    def print_velocity(self) : 
        print("Velocity : [", self.vel[0], "] [", self.vel[1], "] [", self.vel[2], "] \n")

    def print_position(self) : 
        print("*" * 25)
        print("Position : [ x ] [ y ] [ z ] \n")
        self.print_coordinates()
        self.print_velocity()
        print("*" * 25)
            