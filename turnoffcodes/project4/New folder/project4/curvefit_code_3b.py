import math


import numpy as np
import pandas as pd
from scipy import optimize
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def func_cgd_vdg_low(vdg, k1a, m1a, k2a, k3a):
    return k1a / (((1 + (vdg / k2a)) ** m1a) + k3a)

def func_cgd_vdg_high(vdg, k1b, m1b, k2b):
    return k1b / ((((1 + ((vdg - vtd) / k2b))) ** m1b))


def func_coss_vds(vds, k4, k5, m2):
    return k4 / (1 + vds / k5) ** (m2)



# to find threshold using np.piecewise
def func_id_vgs(vgs, vth_th, kp):
    return np.piecewise(vgs, [vgs < vth_th, vgs >= vth_th],
                        [lambda vgs: 0,
                         lambda vgs: (kp * ((vgs - vth_th) ** 2)) / (2 * (1 + (theta * (vgs - vth_th))))])

def func(sec, pvf, y,theta):
    return (0.5 * kp * y * pvf / ((y - 1) * (1 + theta * (sec[:,0] - v_th)))) * (
                (sec[:,0] - v_th) * sec[:, 1] - (pvf ** (y - 1) / y) * ((sec[:,0] - v_th) ** (2 - y)) * (sec[:, 1] ** y))  # sec[:,0]=Vgs, sec[:,1]=Vds

def reverse_transfer_capacitance(file_rtfr):
    vds = pd.read_csv(file_rtfr, header=None, usecols=[0])
    cgd_vds = pd.read_csv(file_rtfr, header=None, usecols=[1])
    vds = vds.squeeze()
    cgd_vds = cgd_vds.squeeze()

    vds = vds.to_numpy()
    cgd_vds = cgd_vds.to_numpy()
    vds = pd.Series(vds)
    cgd_vds = (pd.Series(cgd_vds)) / 1000
    vds_original = vds.copy()
    cgd_vds_original = cgd_vds.copy()


    vgd = np.array(vds)
    cgd = np.array(cgd_vds)
    vds1 = pd.Series(vds)
    cgd_vds1 = (pd.Series(cgd_vds))
    count_value = 0
    vds2 = []
    vds3 = []
    cgd_vds2 = []
    cgd_vds3 = []
    while (count_value <len(vds)):
        # count_value += 1
        if (float(vds1.get(count_value)) >= 1.0):
            # when vdg value is not in the range 0 to 1, the value is stored to a new array vds2
            # corresponding capacitance cgd value is stored into a new array cgd_vds2
            # Take natural logarithm of the values in the new arrays and create two new arrays vds3,cgd_vds3
            vds2.append(vds1[count_value])
            vds3.append(math.log(vds1[count_value]))
            cgd_vds2.append(cgd_vds1[count_value])
            cgd_vds3.append(math.log(cgd_vds1[count_value]))
        count_value += 1


    res = []
    global vtd

    # calculating the slope dcgd by dvdg
    diff_list_vds = np.diff(vds3)
    diff_list_cgd_vds = np.diff(cgd_vds3)
    # res = list(map(truediv, diff_list_cgd_vds , diff_list_vds))
    j = 0
    for j in range(len(vds2) - 1):
        res.append(diff_list_cgd_vds[j] / diff_list_vds[j])


    # finding the index at which the slope is minimum,vdg at that point gives the voltage vtd
    for j in range(len(vds2)):
        if j == res.index(min(res)):
            vtd = vds2[j]


    k = 0
    vdg1 = []
    cgd1 = []
    cgd2 = []
    vdg2 = []
    # vdg1 contains vdg values less than vtd
    # vdg2 contains vdg values greater than or equal to vtd
    for k in range(len(vgd)):
        if vgd[k] <= vtd:
            vdg1.append(vgd[k])
            cgd1.append(cgd[k])
        else:
            vdg2.append(vgd[k])
            cgd2.append(cgd[k])
        # #When vdg values are less than  vtd using eq cgd_vds = k1a / ((1 + (vdg1 / k2a)) ** m1a +k3a)
    params, _ = curve_fit(func_cgd_vdg_low, np.transpose(vdg1), np.transpose(cgd1), maxfev=20000)
    k1a, m1a, k2a, k3a = params[0], params[1], params[2], params[3]

    cgd_vds_fit = k1a / ((1 + (vdg1 / k2a)) ** m1a + k3a)

    # generates figure 3
    plt.subplot(2, 2, 2)
    # plt.subplot2grid((3, 2), (0, 1))
    plt.plot(vdg1, cgd1, marker='o', ms=3, color='black', label="Cgd_vds data")
    plt.plot(vdg1, cgd_vds_fit, label="Cgd_vds = k1a / ((1 + (vdg1 / k2a)) ** m1a +k3a)")
    plt.xlabel('Vds(V)')
    plt.ylabel('Cgd_vds(nF)')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.yscale('log')
    plt.grid(True)

    # #When vdg values are greater than vtd using eq cgd_vds = k1b / (((1 + ((vdg2 - vtd) / k2b)))
    params, _ = curve_fit(func_cgd_vdg_high, np.transpose(vdg2), np.transpose(cgd2), maxfev=100000)
    k1b = params[0]
    m1b = params[1]
    k2b = params[2]
    cgd_vds_fit = k1b / (((1 + ((vdg2 - vtd) / k2b))) ** m1b)

    plt.subplot(2, 2, 2)
    # plt.subplot2grid((3, 2), (0, 1))
    plt.plot(vdg2, cgd2, marker='o', ms=3, color='black')
    plt.plot(vdg2, cgd_vds_fit, label="Cgd_vds = k1b / (((1+((vdg1-vtd)/k2b)))**m1b)")
    plt.xlabel('Vds(V)')
    plt.ylabel('Cgd_vds(nF)')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.yscale('log')
    plt.grid(True)

    return k1b, k2b, m1b, k1a, k2a, m1a, k3a, vtd




def Outputcapacitance(file_dscap):
    vds = pd.read_csv(file_dscap, header=None, usecols=[0])
    coss = pd.read_csv(file_dscap, header=None, usecols=[1])
    vds = vds.squeeze()
    coss = coss.squeeze()
    vds = vds.to_numpy()
    coss = coss.to_numpy()
    vds = pd.Series(vds)
    coss_vds = (pd.Series(coss)) / 1000
    params_coss_vds, _ = curve_fit(func_coss_vds, np.transpose(vds), np.transpose(coss_vds))
    k4, k5, m2 = params_coss_vds[0], params_coss_vds[1], params_coss_vds[2]
    coss_vds_fit = k4 / (1 + vds / k5) ** (m2)

    plt.subplot(2, 2, 1)
    plt.plot(vds, coss_vds, marker='o', ms=3, color='black', label="Coss_vds data")
    plt.plot(vds, coss_vds_fit, label="Coss_vds = k4/(1+vds/k5)**(m2)")
    plt.xlabel('Vds(V)')
    plt.ylabel('Coss_vds(nF)')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.yscale('log')
    plt.grid(True)

    return k4, k5, m2


def transfer_characteristics(vgs_th, id_current_th,theta1):
    vgs_th1 = pd.Series(vgs_th)
    id_current_th1 = pd.Series(id_current_th)
     # Pandas  Series is a     one - dimensional    labeled    array    capable    of    holding    data    of    any    type
    #numpy.squeeze() function is used when we want to remove single-dimensional entries from the shape of an array.
    vgs_th = vgs_th1.squeeze()
    id_current = id_current_th1.squeeze()

    id_current_th_new = np.array(id_current_th)
    vgs_th_new = np.array(vgs_th)
    global theta
    theta = theta1
    params, cov = optimize.curve_fit(func_id_vgs, np.transpose(vgs_th_new), np.transpose(id_current_th_new),
                                     bounds=((0, 0), (np.inf, np.inf)), method='trf', maxfev=100000)

    vth_th, kp = params[0], params[1]
    return vth_th, kp



def output_characteristics(vds3, id3, array1,vth, k_p,theta_old):
    global v_th, kp
    v_th = vth
    kp = k_p
    vds33 = np.transpose(vds3)
    #obtaining combined array of vgs,vds
    a_3d_list1 = zip(array1, vds3)
    a_3d_list = np.array(list(a_3d_list1))
    #calling the curve_fit fuction
        # print("rd before curvefit function is called=",rd[iy])
    if (theta_old == 0):
       thetaprime = np.inf
    else:
       thetaprime = theta_old

    params, pcov = optimize.curve_fit(func, a_3d_list, np.transpose(id3), bounds=((0, 1, 0), (np.inf, 2, thetaprime)), method='trf')
    pvf = params[0]
    y = params[1]
    theta = params[2]
    kf = (y * pvf * 0.5) / (y - 1)

    return pvf, y, theta, kf, a_3d_list



#--------------------data for transfer_characteristics----------------------------------------------

def data_truncation_transfer_characteristics(file_tfr,i_0_max):
    vgs12 = pd.read_csv(file_tfr, header=None, usecols=[0])
    id_current12 = pd.read_csv(file_tfr, header=None, usecols=[1])

    # performing squeeze operation
    vgs111 = vgs12.squeeze()
    id_current111 = id_current12.squeeze()

    # converting to numpy arrays
    vgs11 = vgs111.to_numpy()
    id_current11 = id_current111.to_numpy()
    id_current1 = []
    vgs1=[]
    for y1 in range(0, len(id_current11)):
        if id_current11[y1] < (1.5 * i_0_max):
            id_current1.append(id_current11[y1])
            vgs1.append(vgs11[y1])
        else:
            break

    return vgs1,id_current1

#--------------------data for output_characteristics----------------------------------------------

def data_preparation_output_characteristics(v_gs1,id_vds,i_0_max):
    i = 0
    vds3 = []
    id3 = []
    array1 = []
    array2 = []
    l2 = []
    sum2 = 0
    vds22 = []
    id22 = []
    array2 = []
    while(i < len(v_gs1)):
    # creating arrays for plotting graphs
        vds22.append([])
        id22.append([])
        array2.append([])
        vds2n = []
        id2n = []
        id11 = []
        array2n = []
    # Reading the values from the files
        vds2 = pd.read_csv(id_vds[i], header=None, usecols=[0])
        id2 = pd.read_csv(id_vds[i], header=None, usecols=[1])
        vds = vds2.squeeze()
        id = id2.squeeze()
        vds11 = vds.to_numpy()
        id111 = id.to_numpy()
        id11 = np.array(id111)
        for it in range(0, len(id11)):
            if id11[it] < (i_0_max * 1.5):
                #vds3 has the truncated vds values from all the id_vds files
                vds3.append(vds11[it])
                vds2n.append((vds11[it]))
                # id3 has the truncated id values from all the id_vds files
                id3.append(id11[it])
                id2n.append(id11[it])
            else:
                break
        #l2 is an array with the lengths of vds2n
        l2.append(len(vds2n))
        #vds22[0] ,vds22[1] etc are arrays with the truncated vds values from the individual files
        vds22[i].append(vds2n)
        # id22[0] ,id22[1] etc are arrays with the truncated id values from the individual files
        id22[i].append(id2n)
    # l2 is the array which contains the length of id array with values within the limit

    # forming the array with vgs values
        for pk in range(0, l2[i]):
            array1.append(v_gs1[i])
        #     array2n.append(v_gs1[i])
        # array2[i].append(array2n)
    # reading the vds and id values for creating the graph
        sum2 = sum2 + l2[i]
        i = i + 1
    return vds3, id3, array1, id22, vds22



 # function for finding roots

def equationroots(a, b, c):
    # calculating discriminant using formula
    dis = b * b - 4 * a * c
    sqrt_val = math.sqrt(abs(dis))

    # checking condition for discriminant
    if dis > 0:
        print(" real and different roots ")
        p1 = ((-b + sqrt_val) / (2 * a))
        p2 = ((-b - sqrt_val) / (2 * a))

    elif dis == 0:
        print(" real and same roots")
        p1 = (-b / (2 * a))
        p2 = (-b / (2 * a))

        # when discriminant is less than 0
    else:
        print("Complex Roots")
        p1 = (- b / (2 * a), " + i", sqrt_val)
        p2 = (- b / (2 * a), " - i", sqrt_val)
    return p1, p2
