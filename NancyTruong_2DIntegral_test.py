#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 23:27:21 2023

@author: mexmex
"""

from Final_Project_2D_Integral import Mesh
from numpy import *
import unittest


class TestMesh(unittest.TestCase):
    def setUp(self):
        
        coord = array([[5,5,1,0],
                       [5,1,4,1]])
        
                 
        elements = array([[1,2],
                    [2,3],
                    [3,4]])
        self.mesh = Mesh(coord, elements)
        '''Node 1 = (0,0) (5,5)
        Node 2 = (0,1) (5,1)
        Node 3 = (1,0) (1,4)
        Node 4 = (1,1) (0,1)
        
        Triangle 1 by nodes (1,2,3): (0,0), (0,1), (1,0)
        Triangle 2 by nodes (2,3,4): (0,1), (1,0), (1,1)'''

        

        
        
        coord2 = array([[0,0,1,1],
                        [0,1,0,1]])
        
        self.mesh2 = Mesh(coord2,elements)
  
    def test_one_J(self):
        self.assertAlmostEqual(self.mesh.one_J(0),-16)
        self.assertAlmostEqual(self.mesh.one_J(1),15)
    
    def test_one_edge(self):
        self.assertEqual(self.mesh.one_edge(0)[1, 0], 0)
        self.assertEqual(self.mesh.one_edge(1)[0, 0], 5)
        
    def test_one_edge_zero_length(self):
        
        corner = [
            array([[5, 5], [5, 5], [1, 4]]),
            array([[5, 1], [1, 4], [0, 1]])]
        # Setting the second point same as the first
        self.mesh._corner = corner

        
        with self.assertRaises(Exception) as assert_error:
            self.mesh.one_edge(0)
        
        self.assertEqual(str(assert_error.exception), "Length of edge cant be 0!!!")
        
    def test_one_min_angle(self):
        self.assertAlmostEqual(self.mesh.one_min_angle(0), 0.888, places = 3)
        self.assertAlmostEqual(self.mesh.one_min_angle(1), 0.644, places = 3)
    
    def test_all_J(self):
        self.assertEqual(len(self.mesh.all_J()), 2)
    
    def test_thin_triangle_exception(self):
        #print("TEST THIN HERE")
        
        
        
        with self.assertRaises(Exception) as assert_error:
            thin_triangle_coord = array([[0.0, 1, 2], 
                                         [0.0, 0.0, 0]])
            thin_triangle_elem = array([[1], [2], [3]])
            self.mesh = Mesh(thin_triangle_coord, thin_triangle_elem)
            self.mesh.all_J()

        
        self.assertEqual(str(assert_error.exception), "The triangle is too thin!!!")
    
    def test_I_i(self):
       self.assertAlmostEqual(self.mesh.I_i(lambda x, y: x, 0), 29.333, places = 3)
       self.assertAlmostEqual(self.mesh.I_i(lambda x, y: x+y, 1), 30, places = 3)
        
    def test_I_omega(self):
        self.assertAlmostEqual(self.mesh.I_omega(lambda x, y: x), 44.333, places = 3)  
        self.assertAlmostEqual(self.mesh.I_omega(lambda x, y: x + y), 86, places = 3)
        
        self.assertAlmostEqual(self.mesh2.I_omega(lambda x, y: x), 0.5, places = 3)
        
        self.assertEqual(self.mesh2.I_omega(lambda x, y: 1), 1)
        self.assertEqual(self.mesh2.I_omega(lambda x, y: 2), 2)
        
    def test_all_edge(self):
        self.assertEqual(len(self.mesh.all_edge()), 2)
    
    def test_sum_of_area(self):
        
        self.assertEqual(self.mesh.sum_of_area(), 15.5)
        self.assertEqual(self.mesh2.sum_of_area(), 1)
    
    def test_one_corner(self):
        corner1 = self.mesh.one_corner(0)
        expected1 = array([[5, 5], [5, 1], [1, 4]])
        self.assertTrue(array_equal(corner1, expected1))
    
    def test_all_corner(self):
        
        all_corners = self.mesh.all_corner()
        expected = [
            array([[5, 5], [5, 1], [1, 4]]),
            array([[5, 1], [1, 4], [0, 1]])]
        
        self.assertTrue(array_equal(all_corners, expected))

if __name__ == '__main__':
    unittest.main()