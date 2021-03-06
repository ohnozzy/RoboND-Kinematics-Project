from sympy import *
from time import time
from mpmath import radians
from mpmath import hypot
import tf

'''
Format of test case is [ [[EE position],[EE orientation as quaternions]],[WC location],[joint angles]]
You can generate additional test cases by setting up your kuka project and running `$ roslaunch kuka_arm forward_kinematics.launch`
From here you can adjust the joint angles to find thetas, use the gripper to extract positions and orientation (in quaternion xyzw) and lastly use link 5
to find the position of the wrist center. These newly generated test cases can be added to the test_cases dictionary.
'''

test_cases = {1:[[[2.16135,-1.42635,1.55109],
                  [0.708611,0.186356,-0.157931,0.661967]],
                  [1.89451,-1.44302,1.69366],
                  [-0.65,0.45,-0.36,0.95,0.79,0.49]],
              2:[[[-0.56754,0.93663,3.0038],
                  [0.62073, 0.48318,0.38759,0.480629]],
                  [-0.638,0.64198,2.9988],
                  [-0.79,-0.11,-2.33,1.94,1.14,-3.68]],
              3:[[[-1.3863,0.02074,0.90986],
                  [0.01735,-0.2179,0.9025,0.371016]],
                  [-1.1669,-0.17989,0.85137],
                  [-2.99,-0.12,0.94,4.06,1.29,-4.12]],
              4:[],
              5:[]}

def calculate_angle(a, b, c):
    return acos((b*b+c*c-a*a)/(2*b*c))

def test_code(test_case):
    ## Set up code
    ## Do not modify!
    x = 0
    class Position:
        def __init__(self,EE_pos):
            self.x = EE_pos[0]
            self.y = EE_pos[1]
            self.z = EE_pos[2]
    class Orientation:
        def __init__(self,EE_ori):
            self.x = EE_ori[0]
            self.y = EE_ori[1]
            self.z = EE_ori[2]
            self.w = EE_ori[3]

    position = Position(test_case[0][0])
    orientation = Orientation(test_case[0][1])

    class Combine:
        def __init__(self,position,orientation):
            self.position = position
            self.orientation = orientation

    comb = Combine(position,orientation)

    class Pose:
        def __init__(self,comb):
            self.poses = [comb]

    req = Pose(comb)
    start_time = time()
    
    ########################################################################################
    ## 

    ## Insert IK code here!
    
    roll, pitch, yaw = tf.transformations.euler_from_quaternion(
                   [orientation.x, orientation.y,
                    orientation.z, orientation.w])

    px = req.poses[x].position.x
    py = req.poses[x].position.y
    pz = req.poses[x].position.z
            
    
    

    d = symbols("d0:8")
    a = symbols("a0:7")
    alpha = symbols("alpha0:7") 
    q = symbols("q0:8")
	# Create Modified DH parameters
    s = {
            alpha[0]: 0,    a[0]: 0        , d[1]: 0.75, 
            alpha[1]: -pi/2, a[1]: 0.35     , d[2]: 0,  q[2]: q[2] - pi/2,
            alpha[2]: 0,    a[2]: 1.25     , d[3]: 0,
            alpha[3]: -pi/2, a[3]: -0.054   , d[4]: 1.5,
            alpha[4]: pi/2, a[4]: 0        , d[5]: 0,
            alpha[5]: -pi/2, a[5]: 0        , d[6]: 0,
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
           [ cos(pi/2.),0,-sin(pi/2.),0],
           [         0,1,         0,0],
           [ sin(pi/2.),0, cos(pi/2.),0],
           [         0,0,         0,1]
          ])
    R_cor = simplify(R_z * R_y)
    T_total = simplify(T0 * R_cor)

    print("Start IK")
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
    ROT_EE = ROT_EE.subs({'r': roll, 'p': pitch, 'y': yaw})
            
    EE = Matrix([[px],[py],[pz]])
    WC = EE - 0.303 * ROT_EE[:,2]
	    #
	    # Calculate joint angles using Geometric IK method
	    #
	    #
            ###
    theta1 = atan2(WC[1], WC[0])
    print("theta1")   
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
    print("theta2-3")        
    R0_3_analysis = T[0][0:3, 0:3]*T[1][0:3, 0:3]*T[2][0:3, 0:3]
    R0_3 = R0_3_analysis.evalf(subs={q[1]: theta1, q[2]: theta2, q[3]: theta3})
    R3_6 = R0_3.inv("LU")*ROT_EE

    theta4 = atan2(R3_6[2,2], -R3_6[0,2])
    theta5 = atan2(hypot(R3_6[0,2],R3_6[2,2]),R3_6[1,2])
    theta6 = atan2(-R3_6[1,1], R3_6[1,0])
    print("theta4-6")        
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
    #joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5,theta6]    


    ## 
    ########################################################################################
    ########################################################################################
    ## For additional debugging add your forward kinematics here. Use your previously calculated thetas
    ## as the input and output the position of your end effector as your_ee = [x,y,z]
    theta_map = {q[1]: test_case[2][0], q[2]: test_case[2][1], q[3]: test_case[2][2], q[4]:test_case[2][3] ,q[5]:test_case[2][4], q[6]:test_case[2][5]} 

    ## End your code input for forward kinematics here!
    ########################################################################################

    ## For error analysis please set the following variables of your WC location and EE location in the format of [x,y,z]
    your_wc = [1,1,1] # <--- Load your calculated WC values in this array
    your_ee = T_total.subs(theta_map).col(-1) # <--- Load your calculated end effector value from your forward kinematics
    ########################################################################################

    ## Error analysis
    print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time))

    # Find WC error
    if not(sum(your_wc)==3):
        wc_x_e = abs(your_wc[0]-test_case[1][0])
        wc_y_e = abs(your_wc[1]-test_case[1][1])
        wc_z_e = abs(your_wc[2]-test_case[1][2])
        wc_offset = sqrt(wc_x_e**2 + wc_y_e**2 + wc_z_e**2)
        print ("\nWrist error for x position is: %04.8f" % wc_x_e)
        print ("Wrist error for y position is: %04.8f" % wc_y_e)
        print ("Wrist error for z position is: %04.8f" % wc_z_e)
        print ("Overall wrist offset is: %04.8f units" % wc_offset)

    # Find theta errors
    t_1_e = abs(theta1-test_case[2][0])
    t_2_e = abs(theta2-test_case[2][1])
    t_3_e = abs(theta3-test_case[2][2])
    t_4_e = abs(theta4-test_case[2][3])
    t_5_e = abs(theta5-test_case[2][4])
    t_6_e = abs(theta6-test_case[2][5])
    print ("\nTheta 1 error is: %04.8f" % t_1_e)
    print ("Theta 2 error is: %04.8f" % t_2_e)
    print ("Theta 3 error is: %04.8f" % t_3_e)
    print ("Theta 4 error is: %04.8f" % t_4_e)
    print ("Theta 5 error is: %04.8f" % t_5_e)
    print ("Theta 6 error is: %04.8f" % t_6_e)
    print ("\n**These theta errors may not be a correct representation of your code, due to the fact \
           \nthat the arm can have muliple positions. It is best to add your forward kinmeatics to \
           \nconfirm whether your code is working or not**")
    print (" ")

    # Find FK EE error
    if not(sum(your_ee)==3):
        ee_x_e = abs(your_ee[0]-test_case[0][0][0])
        ee_y_e = abs(your_ee[1]-test_case[0][0][1])
        ee_z_e = abs(your_ee[2]-test_case[0][0][2])
        ee_offset = sqrt(ee_x_e**2 + ee_y_e**2 + ee_z_e**2)
        print ("\nEnd effector error for x position is: %04.8f" % ee_x_e)
        print ("End effector error for y position is: %04.8f" % ee_y_e)
        print ("End effector error for z position is: %04.8f" % ee_z_e)
        print ("Overall end effector offset is: %04.8f units \n" % ee_offset)




if __name__ == "__main__":
    # Change test case number for different scenarios
    test_case_number = 1

    test_code(test_cases[test_case_number])
