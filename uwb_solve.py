#!/usr/bin/env python
import rospy
import numpy as np
from nlink_parser.msg import LinktrackTagframe0
from gps_common.msg import GPSFix
from math import sin, cos, sqrt, atan2, pow, radians, degrees

# Constants
M_PI = 3.14159265358979323846
deg = M_PI / 180.0
f = 1.0 / 298.257
e = sqrt(f * 2 - pow(f, 2))
Re = 6378137
Rp = (1 - f) * Re
ep = sqrt(pow(Re, 2) - pow(Rp, 2)) / Rp

# Variables
d0, d1, d2, d3 = 0, 0, 0, 0

# GPS data
Lon = np.array([118.70418395375, 118.70423659438, 118.70415988838, 118.70421629093, 118.70421313409])
Lat = np.array([32.15520190344, 32.15520817160, 32.15547203086, 32.15547707852, 32.1551839796])
Alt = np.array([16.5119 + 0.8, 16.4566 + 0.8, 16.4996 + 0.8, 16.4159 + 0.8, 15.6761])

N = Re / np.sqrt(1 - pow(e, 2) * np.sin(Lat * deg) * np.sin(Lat * deg))
X = (N + Alt) * np.cos(Lat * deg) * np.cos(Lon * deg)
Y = (N + Alt) * np.cos(Lat * deg) * np.sin(Lon * deg)
Z = (N * (1 - pow(e, 2)) + Alt) * np.sin(Lat * deg)

Cet = np.array([
    [-sin(Lon[1] * deg), cos(Lon[1] * deg), 0],
    [-sin(Lat[1] * deg) * cos(Lon[1] * deg), -sin(Lat[1] * deg) * sin(Lon[1] * deg), cos(Lat[1] * deg)],
    [cos(Lat[1] * deg) * cos(Lon[1] * deg), cos(Lat[1] * deg) * sin(Lon[1] * deg), sin(Lat[1] * deg)]
])

C2 = np.array([X[0], Y[0], Z[0]])
p = np.zeros((5, 3))
for i in range(5):
    C1 = np.array([X[i], Y[i], Z[i]])
    p[i] = np.dot(Cet, (C1 - C2)).T

Localx = p[:, 0]
Localy = p[:, 1]
Localz = p[:, 2]

ekf_uwb = p[4, :].reshape(3, 1)
A = np.eye(3)
P = np.array([[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.2]])
Q = np.array([[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]])
R = np.array([[5.0, 0.0, 0.0, 0.0], [0.0, 5.0, 0.0, 0.0], [0.0, 0.0, 5.0, 0.0], [0.0, 0.0, 0.0, 5.0]])

# ROS Callback
def tagframe0Callback(msg):
    global d0, d1, d2, d3
    d0 = msg.dis_arr[1]
    d1 = msg.dis_arr[2]
    d2 = msg.dis_arr[4]
    d3 = msg.dis_arr[5]
    rospy.loginfo("Enter success!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

# Main function
def main():
    rospy.init_node('uwb_dis2llt', anonymous=True)
    rospy.Subscriber("nlink_linktrack_tagframe0", LinktrackTagframe0, tagframe0Callback)
    pub = rospy.Publisher("uwb_ekf", GPSFix, queue_size=1000)
    rate = rospy.Rate(50)

    global ekf_uwb, P, Q, R, A
    position = GPSFix()

    while not rospy.is_shutdown():
        rospy.spinOnce()

        # Prediction Step
        ekf_uwb = np.dot(A, ekf_uwb)

        # Calculate H matrix
        D0 = sqrt(pow(ekf_uwb[0] - Localx[0], 2) + pow(ekf_uwb[1] - Localy[0], 2) + pow(ekf_uwb[2] - Localz[0], 2))
        D1 = sqrt(pow(ekf_uwb[0] - Localx[1], 2) + pow(ekf_uwb[1] - Localy[1], 2) + pow(ekf_uwb[2] - Localz[1], 2))
        D2 = sqrt(pow(ekf_uwb[0] - Localx[2], 2) + pow(ekf_uwb[1] - Localy[2], 2) + pow(ekf_uwb[2] - Localz[2], 2))
        D3 = sqrt(pow(ekf_uwb[0] - Localx[3], 2) + pow(ekf_uwb[1] - Localy[3], 2) + pow(ekf_uwb[2] - Localz[3], 2))

        H = np.array([
            [(ekf_uwb[0] - Localx[0]) / D0, (ekf_uwb[1] - Localy[0]) / D0, (ekf_uwb[2] - Localz[0]) / D0],
            [(ekf_uwb[0] - Localx[1]) / D1, (ekf_uwb[1] - Localy[1]) / D1, (ekf_uwb[2] - Localz[1]) / D1],
            [(ekf_uwb[0] - Localx[2]) / D2, (ekf_uwb[1] - Localy[2]) / D2, (ekf_uwb[2] - Localz[2]) / D2],
            [(ekf_uwb[0] - Localx[3]) / D3, (ekf_uwb[1] - Localy[3]) / D3, (ekf_uwb[2] - Localz[3]) / D3]
        ])

        # Update P matrix
        P = np.dot(np.dot(A, P), A.T) + Q

        # Calculate Kalman Gain
        K = np.dot(np.dot(P, H.T), np.linalg.inv(np.dot(np.dot(H, P), H.T) + R))

        # Measurements
        Z_UWB = np.array([d0, d1, d2, d3]).reshape(4, 1)
        UZ = np.array([D0, D1, D2, D3]).reshape(4, 1)

        # Update step
        ekf_uwb = ekf_uwb + np.dot(K, (Z_UWB - UZ))

        # Update P matrix
        P = np.dot((np.eye(3) - np.dot(K, H)), P)

        # Local to ECEF
        ecef_p = np.dot(Cet.T, ekf_uwb).flatten() + C2

        # ECEF to LLA
        v = sqrt(pow(ecef_p[0], 2) + pow(ecef_p[1], 2))
        theta = atan2(ecef_p[2] * Re, v * Rp)
        lon = degrees(atan2(ecef_p[1], ecef_p[0]))
        lat = degrees(atan2(ecef_p[2] + pow(ep, 2) * Rp * pow(sin(theta), 3), v - pow(e, 2) * Re * pow(cos(theta), 3)))
        NN = Re / sqrt(1 - pow(e, 2) * sin(radians(lat)) * sin(radians(lat)))
        alt = v / cos(radians(lat)) - NN

        # Publish position
        position.header.stamp = rospy.Time.now()
        position.longitude = lon
        position.latitude = lat
        position.altitude = alt
        rospy.loginfo("lon: %s lat: %s", lon, lat)
        pub.publish(position)

        rate.sleep()

if __name__ == '__main__':
    try:
        main()
    except rospy.ROSInterruptException:
        pass
