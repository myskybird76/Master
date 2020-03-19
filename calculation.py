#x = 1 + 3 * ((4.5 - 4.0) / (4.5 - 1.0)) # 1.86 1
#y = ((4.5 - 1) * ((2.05 - 1) / 3)) - 4.5
# 3.24 / 4.30
# 3.49 / 4.50
#print(x)
#print(y)

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy import optimize

data = np.genfromtxt('D:/Master thesis/Final_refill/result_final1.csv', dtype = None, delimiter=',', skip_header = 1, encoding = 'UTF8', names = ('Name', 'c', 'cp', 'cm', 'w', 'wp', 'wm', 'G', 'Gp', 'Gm', 'M', 'Mp', 'Mm', 'A1', 'A1p', 'A1m', 'A2', 'A2p', 'A2m'
                                                                                                                                                 , 'P2', 'P2p', 'P2m', 'P3', 'P3p', 'P3m', 'P4', 'P4p', 'P4m'))


for i in range(len(data)) :

    Name = data['Name'][i]

    w = data['w'][i]
    w_p = data['wp'][i]
    w_m = data['wm'][i]
    w1 = '{:f}'.format((w_p + w_m) / 2)


    c = data['c'][i]
    c_p = data['cp'][i]
    c_m = data['cm'][i]
    c1 = '{:f}'.format((c_p + c_m) / 2)


    G = data['G'][i]
    G_p = data['Gp'][i]
    G_m = data['Gm'][i]
    G1 = '{:f}'.format((G_p + G_m) / 2)

    M = data['M'][i]
    M_p = data['Mp'][i]
    M_m = data['Mm'][i]
    M1 = '{:f}'.format((M_p + M_m) / 2)

    A180 = data['A1'][i]
    A180_p = data['A1p'][i]
    A180_m = data['A1m'][i]
    A1801 = '{:f}'.format((A180_p + A180_m) / 2)

    A90 = data['A2'][i]
    A90_p = data['A2p'][i]
    A90_m = data['A2m'][i]
    A901 = '{:f}'.format((A90_p + A90_m) / 2)

    P2 = data['P2'][i]
    P2_p = data['P2p'][i]
    P2_m = data['P2m'][i]
    P21 = (P2_p + P2_m) / 2


    P3 = data['P3'][i]
    P3_p = data['P3p'][i]
    P3_m = data['P3m'][i]
    P31 = (P3_p + P3_m) / 2

    P4 = data['P4'][i]
    P4_p = data['P4p'][i]
    P4_m = data['P4m'][i]
    P41 = (P4_p + P4_m) / 2
    print(
        r'{} & \num[round-precision=3]{} & \(\pm\) & \num[round-precision=3]{} & \num[round-precision=3]{} & \(\pm\) & \num[round-precision=3]{} & \num[round-precision=3]{} & \(\pm\) & \num[round-precision=3]{} & \num[round-precision=3]{} & \(\pm\) & \num[round-precision=3]{} & \num[round-precision=3]{} & \(\pm\) & \num[round-precision=3]{} & \num[round-precision=3]{} & \(\pm\) & \num[round-precision=3]{} & \num[round-precision=2]{} & \(\pm\) & \num[round-precision=2]{} & \num[round-precision=2]{} & \(\pm\) & \num[round-precision=2]{} & \num[round-precision=2]{} & \(\pm\) & \num[round-precision=2]{} &  &  \\'.format(
            Name, {c}, {c1}, {w}, {w1}, {G}, {G1}, {M}, {M1}, {A180}, {A1801}, {A90}, {A901}, {P2}, {P21}, {P3}, {P31}, {P4}, {P41}))



