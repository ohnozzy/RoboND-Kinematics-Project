#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
from sympy.solvers import solve

def calculate_angle(a, b, c):
    return acos((b*b+c*c-a*a)/(2*b*c))

def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:

        ### Your FK code here
        # Create symbols
	d = symbols("d0:8")
        a = symbols("a0:7")
        alpha = symbols("alpha0:7") 
        q = symbols("q0:8")
	# Create Modified DH parameters
	s = {
            alpha[0]: 0,    a[0]: 0        , d[1]: 0.75, 
            alpha[1]: -pi/2., a[1]: 0.35     , d[2]: 0,  q[2]: q[2] - pi/2,
            alpha[2]: 0,    a[2]: 1.25     , d[3]: 0,
            alpha[3]: -pi/2., a[3]: -0.054   , d[4]: 1.5,
            alpha[4]: pi/2, a[4]: 0        , d[5]: 0,
            alpha[5]: -pi/2., a[5]: 0        , d[6]: 0,
            alpha[6]: 0   , a[6]: 0        , d[7]: 0.303, q[7]: 0
        }
	
	# Define Modified DH Transformation matrix
        T=[]
	for i in range(1, 8):
          qi = q[i]
          ai = a[i-1]
          alphai = alpha[i-1]
          di = d[i]
          T.append(Matrix([
           [            cos(qi),            -sin(qi),            0,              ai],
           [sin(qi)*cos(alphai), cos(qi)*cos(alphai), -sin(alphai), -sin(alphai)*di],
           [sin(qi)*sin(alphai), cos(qi)*sin(alphai), cos(alphai), cos(alphai)*di],
           [                  0,                   0,            0,               1] 
          ]).subs(s))
	T0=T[0]
        for i in range(1, len(T)):
          T0 = simplify(T0 * T[i])
        
        R_z = Matrix([
           [ cos(pi), -sin(pi), 0, 0],
           [ sin(pi),  cos(pi), 0, 0],
           [       0,        0, 1, 0],
           [       0,        0, 0, 1] 
          ])
        R_y = Matrix([
           [ cos(pi/2),0,-sin(pi/2),0],
           [         0,1,         0,0],
           [ sin(pi/2),0, cos(pi/2),0],
           [         0,0,         0,1]
          ])
        R_cor = simplify(R_z * R_y)
        T_total = simplify(T0 * R_cor)
        ###

        # Initialize service response
        joint_trajectory_list = []

        r, p, y = symbols('r p y')
        ROT_x = Matrix([[1,          0,       0],
                        [0,     cos(r), -sin(r)],
                        [0,     sin(r),  cos(r)]])

        ROT_y = Matrix([[cos(p),     0,  sin(p)],
                         [0,          1,       0],
                         [-sin(p),    0,  cos(p)]])
            
        ROT_z = Matrix([[cos(y), -sin(y),     0],
                        [sin(y),  cos(y),     0],
                        [0,            0,     1]])
            
        ROT_EE = ROT_z * ROT_y * ROT_x

        Rot_Error = ROT_z.subs(y, radians(180)) * ROT_y.subs(p, radians(-90))
        ROT_EE = ROT_EE * Rot_Error
        R0_3_analysis = T[0][0:3, 0:3]*T[1][0:3, 0:3]*T[2][0:3, 0:3]
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            ### Your IK code here
	    # Compensate for rotation discrepancy between DH parameters and Gazebo
            
            ROT_EE_sub_r = ROT_EE.subs({'r': roll, 'p': pitch, 'y': yaw})
            
            EE = Matrix([[px],[py],[pz]])
            WC = EE - 0.303 * ROT_EE_sub_r[:,2]
	    #
	    # Calculate joint angles using Geometric IK method
	    #
	    #
            ###
            theta1 = atan2(WC[1], WC[0])
       
            side_a = 1.501
            translate_xy = hypot(WC[0], WC[1]) - 0.35
            translate_z  = WC[2] - 0.75
            
            side_b = hypot(translate_xy, translate_z)
            side_c = 1.25
            
             
            angle_a = calculate_angle(side_a, side_b, side_c)
            angle_b = calculate_angle(side_b, side_a, side_c)
            angle_c = calculate_angle(side_c, side_a, side_b)

            theta2 = pi / 2 - angle_a - atan2(translate_z, translate_xy)
            theta3 = pi / 2 - angle_b - 0.036
            
            
            R0_3 = R0_3_analysis.evalf(subs={q[1]: theta1, q[2]: theta2, q[3]: theta3})
            R3_6 = R0_3.inv("LU")*ROT_EE_sub_r

            theta4 = atan2(R3_6[2,2], -R3_6[0,2])
            theta5 = atan2(hypot(R3_6[0,2],R3_6[2,2]),R3_6[1,2])
            theta6 = atan2(-R3_6[1,1], R3_6[1,0])
            
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5,theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
