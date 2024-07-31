
#https://blog.csdn.net/qq_34213260/article/details/109133847

#wgs84坐标系转换为ECEF坐标系
import pyproj
def lla2ecef(lon, lat, alt):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    x, y, z = pyproj.transform(lla,ecef, lon, lat, alt)
    return x, y, z

#ECEF坐标系转换为wgs84坐标系
def ecef2lla(x,y,z):
    #ecef转换为经纬高
    ecef = pyproj.Proj(proj='geocent',ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong',ellps='WGS84', datum='WGS84')
    lon, lat, alt = pyproj.transform(ecef,lla,x,y,z,radians=False) #用弧度返回值
    return lat,lon,alt

######################手动计算
from math import radians, cos, sin

def lla_to_xyz(lat, lon, alt):
  '''
  将"lla"坐标系的经纬度转换为"xyz"坐标系的笛卡尔坐标
  '''
    a = 6378137.0  # 地球长半径
    f = 1 / 298.257223563  # 地球扁率
    e2 = 2*f - f**2  # 第一偏心率的平方

    # 将经纬度转换为弧度
    lat_rad = radians(lat)
    lon_rad = radians(lon)

    # 计算地球的第一卯线曲率半径
    N = a / (1 - e2*sin(lat_rad)**2)**0.5

    # 计算笛卡尔坐标系的坐标
    x = (N + alt) * cos(lat_rad) * cos(lon_rad)
    y = (N + alt) * cos(lat_rad) * sin(lon_rad)
    z = ((1 - e2) * N + alt) * sin(lat_rad)

    return x, y, z


def xyz_to_enu(x, y, z, lat_ref, lon_ref):
'''
将"xyz"坐标系转换为"enu"坐标系
'''
    lat_ref_rad = radians(lat_ref)
    lon_ref_rad = radians(lon_ref)

    # 计算转换矩阵
    T = [[-sin(lon_ref_rad), cos(lon_ref_rad), 0],
         [-sin(lat_ref_rad)*cos(lon_ref_rad), -sin(lat_ref_rad)*sin(lon_ref_rad), cos(lat_ref_rad)],
         [cos(lat_ref_rad)*cos(lon_ref_rad), cos(lat_ref_rad)*sin(lon_ref_rad), sin(lat_ref_rad)]]

    # 计算ENU坐标
    e = T[0][0]*x + T[0][1]*y + T[0][2]*z
    n = T[1][0]*x + T[1][1]*y + T[1][2]*z
    u = T[2][0]*x + T[2][1]*y + T[2][2]*z

    return e, n, u



