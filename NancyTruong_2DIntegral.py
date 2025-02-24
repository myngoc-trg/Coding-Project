# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 14:06:52 2023

@author: leung
"""



from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import norm

class Mesh:
    def __init__(self, *mesh):
        self._coor, self._elem =mesh
        self._corner = self.all_corner()
        self._edge = self.all_edge()
        self._dete = self.all_J()
        
    
    def one_corner(self,i):
        CoorT = self._coor.transpose()
        corner = array([CoorT[int(self._elem[j,i])-1] for j in range(3)])
        #print(corner)
        return corner
    
    def all_corner(self):
        corner = [self.one_corner(i) for i in range(shape(self._elem)[1])]
        return corner    
    
    def one_edge(self,i):
        corner = self._corner[i]
        
        edge = array([corner[j]-corner[j-1] for j in range(3)])
        for ed in edge:
            if norm(ed) ==0:
                raise Exception("Length of edge cant be 0!!!")
        return edge
    
    def all_edge(self):
        edge = [self.one_edge(i) for i in range(shape(self._elem)[1])]
    
        return edge
    
    
    def one_J(self, i):
        corner = self._corner[i]
        r1 = corner[1] - corner[0]
        r2 = corner[2] - corner[0]
        
        M = vstack((r1,r2))
        #print(M)
        return linalg.det(M)
    
    def one_min_angle(self, i):
                
        edge = self._edge[i]
        a = [arccos(-(edge[j] @ edge[j-1] )/(norm(edge[j])*norm(edge[j-1])) ) for j in range(3)]
        ''' 
        for aa in a:
            print(aa)
         '''   
        return min(a)
    
    def all_J(self):
        for i in range(shape(self._elem)[1]):
            if self.one_min_angle(i) <1.e-14:
                raise Exception("The triangle is too thin!!!")
        J = [self.one_J(i) for i in range(shape(self._elem)[1])]
       
        return J 
    
  
    def sum_of_area(self):
        area = [abs(cross(self._edge[i][0],self._edge[i][1]))/2 for i in range(shape(self._elem)[1]) ]
        return sum(area)
   
    def I_i(self,f,i):
        x = [self._coor[0,int(self._elem[j,i])-1] for j in range(3)]
        y = [self._coor[1,int(self._elem[j,i])-1] for j in range(3)]
        I = (abs(self._dete[i])/6)*(f(x[0],y[0])+f(x[1],y[1])+f(x[2],y[2])) 
       # print(I)
        return I
    
    def I_omega(self,f):
        I = [self.I_i(f,i) for i in range(shape(self._elem)[1])]
        return sum(I)
    
    def plot(self):
        color= ['red', 'orange', 'green', 'blue','purple', 'black']
        
        r = 0.1
        fig, ax = subplots()
        for i in range(shape(self._elem)[1]):
           
            
            corner = self._corner[i]
            edge = self._edge[i]
                        
            edge_length = [norm(ed) for ed in edge]
            parameter = sum(edge_length)
            weighted_corner = [corner[j]*edge_length[j-1] for j in range(3)]
            weighted_x = [weighted_corner[j][0] for j in range(3)]
            weighted_y = [weighted_corner[j][1] for j in range(3)]
            incentre = array([sum(weighted_x)/parameter, sum(weighted_y)/parameter])
            #print(incentre)
            direction_vector = [incentre - corner[j] for j in range(3)]
            a = [arccos(-(edge[j] @ edge[j-2] )/(norm(edge[j])*norm(edge[j-2])) ) for j in range(3)]
            #print(a)
            
            move_vector = [r*direction_vector[j] for j in range(3)]
                        
            new_coor = array([corner[j]+move_vector[j] for j in range(3)])
           
            x = [new_coor[j,0] for j in range(-1,3)]
            y = [new_coor[j,1] for j in range(-1,3)]
                       
            ax.plot(x, y, color=color[i%len(color)], linewidth =1)
               
        ax.axis('equal')
    

        
def testF(x,y):
    return sqrt(x**2+y**2)

    
with open ('coord1.txt', 'r') as myfile:
    coord = loadtxt(myfile)

with open ('elementnode1.txt', 'r') as myfile:
    elements = loadtxt(myfile)


#print(coord)
#print (elements)

ex1 = Mesh(coord, elements)

sum_of_ex1 = ex1.sum_of_area()
I_omega_ex1_1 = ex1.I_omega(lambda x,y :1)

print(f"The total area is {sum_of_ex1} and the integral of the constant function one is {I_omega_ex1_1}. ")

ex1.plot()


with open ('coordinates_dolfin_coarse.txt', 'r') as myfile:
    coord_d = loadtxt(myfile)

with open ('nodes_dolfin_coarse.txt', 'r') as myfile:
    elements_d = loadtxt(myfile)



dophine = Mesh(coord_d, elements_d)

dophine.plot()

print(f"The area of dophine is {1-dophine.sum_of_area()}.")

with open ('coordinates_unitcircle_10000.txt', 'r') as myfile:
    coord_c = loadtxt(myfile)

with open ('nodes_unitcircle_10000.txt', 'r') as myfile:
    elements_c = loadtxt(myfile)



circle = Mesh(coord_c, elements_c)

circle.plot()

print(f"The area of circle is {circle.sum_of_area()}.")

v_testF = circle.I_omega(testF)

print(f"The volume is {v_testF}.")

