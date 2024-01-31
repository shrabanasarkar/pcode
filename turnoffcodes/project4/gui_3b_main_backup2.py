import PySimpleGUI as sg
from scipy.integrate import quad

from curvefit_code_3b import *
from gui_second import *
from scipy import integrate
from scipy.integrate import quad


# layout = [[sg.Text('Turn off Switching Loss Calculator for SiC MOSFET in Half Bridge Configuration', auto_size_text=False,
#                    justification='center', font=("Arial", 20), size=(80, 1))],
#           [sg.Text('OPERATING CONDITION', font=("Arial", 16), justification='center')],
#           [sg.Text('DC bus voltage (Vdc)', font=(30), size=(20, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key = 'Vdc',justification='left', size=(5, 1)),
#            sg.Text('(V)', font=(30), size=(3, 1)), sg.Text('', size=(2, 1)),
#            sg.Text('Minimum Load current (Iomin)', font=(30), size=(30, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key='Iomin',justification='left', size=(5, 1)),
#            sg.Text('(A)', font=(30), size=(3, 1)), sg.Text('', size=(2, 1)),
#            sg.Text('Maximum Load current (Iomax)', font=(30), size=(30, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key='Iomax' ,justification='left', size=(5, 1)),
#            sg.Text('(A)', font=(30), size=(3, 1))],
#             [sg.Text('External capacitance', font=(30), size=(20, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key = 'Cext',justification='left', size=(5, 1)),
#            sg.Text('(pF)', font=(30), size=(3, 1)), sg.Text('', size=(2, 1)),],
#           # [sg.Text('DC bus voltage', font=(30), size=(13, 1)), sg.InputText(font=(30), do_not_clear=False, justification='left', size=(5,1)), sg.Text('(V)', font=(30), size=(3, 1)), sg.Text('', size=(8,1)), sg.Text('Load current', font=(30), size=(12, 1)), sg.InputText(font=(30), do_not_clear=False, justification='left', size=(5,1)), sg.Text('(A)', font=(30), size=(3, 1))],
#           [sg.Text('GATE DRIVER PARAMETERS', font=("Arial", 16), justification='center')],
#           [sg.Text('Vgg', font=(30), size=(5, 1)), sg.InputText(font=(30), do_not_clear=True,key ='Vgg', justification='left', size=(5,1)), sg.Text('(V)', font=(30),
#            size=(5, 1)), sg.Text('', size=(2,1)), sg.Text('Vee', font=(30), size=(5, 1)), sg.InputText(font=(30), do_not_clear=True,key ='Vee', justification='left', size=(5,1)),
#            sg.Text('(V)', font=(30), size=(3, 1)), sg.Text('', size=(2,1)), sg.Text('Rgext', font=(30), size=(7, 1)),
#
#            sg.InputText(font=(30), do_not_clear=True,key ='Rgext', justification='left', size=(5, 1)),
#            sg.Text('(Ω)', font=(30), size=(5, 1))
#            ],
#           [sg.Text('EXTERNAL CIRCUIT PARASITICS', font=("Arial", 16), justification='center')],
#           [sg.Text('Ld', font=(30), size=(12, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key='Ld', justification='left', size=(5, 1)),
#            sg.Text('(nH)', font=(30), size=(5, 1)), sg.Text('', size=(5, 1)), sg.Text('Ls', font=(30), size=(12, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key='Ls', justification='left', size=(5, 1)),
#            sg.Text('(nH)', font=(30), size=(8, 1)), sg.Text('Cgd(ext)', font=(30), size=(12, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key ='Cgd(ext)', justification='left', size=(5, 1)),
#            sg.Text('(pF)', font=(30), size=(4, 1)), sg.Text('', size=(2, 1)),
#             ],
#           [sg.Text('SiC MOSFET PARAMETERS', font=("Arial", 16), justification='center')],
#           [sg.Text('Ron', font=(30), size=(6, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key ='Ron', justification='left', size=(5, 1)),
#            sg.Text('(Ω)', font=(30), size=(5, 1)), sg.Text('', size=(2, 1)),
#            sg.Text('Rgint', font=(30), size=(6, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key = 'Rgint', justification='left', size=(5, 1)),
#            sg.Text('(Ω)', font=(30), size=(5, 1)), sg.Text('', size=(2, 1)),
#            sg.Text('Ciss', font=(30), size=(6, 1)),
#            sg.InputText(font=(30), do_not_clear=True,key = 'Ciss', justification='left', size=(5, 1)),
#            sg.Text('(nF)', font=(30), size=(3, 1)), sg.Text('', size=(2, 1)),],
#           [sg.Text('Dataset for output characteristics', font=(30), justification='center'),sg.Button('Insert',font=(30))],
#           [sg.Text('Dataset for Transfer Characteristics (Id(A)-Vgs(V)) (.csv):', size=(57, 1), font=(30)),
#            sg.Input(key='file for  transfer characteristics',do_not_clear=True), sg.FileBrowse(font=(30))],
#           [sg.Text('Dataset for Reverse Transfer capacitance (Crss(pF)-Vds(V)) (.csv):', size=(57, 1), font=(30)),
#            sg.Input(key='file for reverse transfer characteristics'), sg.FileBrowse(font=(30))],
#           [sg.Text('Dataset for Output Capacitance (Coss(pF)-Vds(V)) (.csv):', size=(57, 1), font=(30)),
#            sg.Input(key='file for Output Capacitance',do_not_clear=True), sg.FileBrowse(font=(30))],
#           [sg.Text('', font=(30), size=(20, 1)), sg.Submit(font=(30)), sg.Text('', font=(30), size=(10, 1)),
#            sg.Cancel(font=(30))]]
# window_first = sg.Window('Calculator', layout)
#
#
#
# #---Start of the main function----
# count = 0
# while True :
#     #count is set to 1 , when user enters the data in the main window fields for the first time
#     count = count + 1
#     # when all field values are correct b is set to 2
#     b =2
#     a = 0
#     event, values = window_first.Read()
#     if event == "Cancel" or event == sg.WIN_CLOSED:
#         break
#     global v_gs1
#     global id1_vds
#     global flagb
#     #---if the user is entering the values in the main window for the first time , the dataset values are empty,so flagb is set to 1
#     #---Then, if user tries to submit with those empty values, the popup will come showing that outdataset characteristics should not be empty
#     if count ==1:
#         v_gs1 = []
#         id1_vds = []
#         flagb =1
#
#     if  event == 'Insert':
#         v_gs12,id1_vds2 = second_window()
#         v_gs1 = []
#         id1_vds = []
#         for t1 in range(len(v_gs12)):
#             if v_gs12[t1] != '':
#                 v_gs1.append(float(v_gs12[t1]))
#         for t12 in range(len(id1_vds2)):
#             if id1_vds2[t12] != '':
#                 id1_vds.append(id1_vds2[t12])
#         # when the v_gs and id1_vds values are empty ,flagb value is set to 1
#         if v_gs1 == [] or id1_vds == []:
#             flagb = 1
#         else:
#             flagb = 0
#     if  event=='Submit':
#         b1 = ''
#         # c has the count of values that are not provided by the user
#         c = 0
#         for v, l in values.items():
#             if values[v] == '':
#                 c = c + 1
#                 b1 = b1 + v + ','
#         if b1 != '':
#             te = 0
#             while (b1.rfind('Browse') != -1):
#                 #te has the count of the number of times the string 'Browse' appears in a string
#                 te = te + 1
#                 p = b1.rfind('Browse')
#                 b2 = str(b1[p + 8:])
#
#                 if (te == 4):
#                     b1 = str(b1[:p]) + 'f' + b2
#                 else:
#                     b1 = str(b1[:p]) + b2
#             tr = b1.rfind(',')
#             #b is set to 0 when one or more fields are empty
#             b=0
#             if b==0 and flagb == 1:
#                 sg.popup(b1[:tr] + " and dataset for output characteristics  are missing.")
#             if b == 0 and c >1 and flagb!=1:
#                 sg.popup(b1[:tr] + " are missing.")
#             if b == 0 and c ==1 and flagb!=1:
#                 sg.popup(b1[:tr] + " is missing.")
#             # if some values are missing in the gui , then b is set to 0
#
#         else:
#             if float(values['Vdc']) <= 0:
#                 sg.popup("Enter correct value of Vdc")
#                 b = 1
#                 # if some values are incorrect in the gui , then b is set to 1
#
#             if float(values['Iomin']) <= 0:
#                 sg.popup("Enter correct value of Iomin")
#                 # if value of Iomin is  incorrect in the gui , then a is set to 1
#                 a = 1
#                 b = 1
#
#             if float(values['Iomax']) <= 0:
#                 sg.popup("Enter correct value of Iomax")
#                 # if value of Iomax is  incorrect in the gui , then a is set to 1
#                 b = 1
#                 a = 1
#
#             if float(values['Iomin']) >= float(values['Iomax']) and a!=1:
#                 sg.popup("Enter correct values of Iomin & Iomax")
#                 b = 1
#
#             if float(values['Vgg']) <= 0:
#                 sg.popup("Enter correct value of Vgg")
#                 # if value of Vgg is  incorrect in the gui , then q is set to 1
#                 b = 1
#                 q = 1
#
#             if float(values['Vee']) > 0:
#                 sg.popup("Enter correct value of Vee")
#                 # if value of Vee is  incorrect in the gui , then q is set to 1
#                 b = 1
#                 q = 1
#
#             if float(values['Vee']) >= float(values['Vgg']) and q!=1:
#                 sg.popup("Enter correct values of Vee & Vgg")
#                 b = 1
#
#             if float(values['Rgext']) <= 0:
#                 sg.popup("Enter correct value of Rgext")
#                 b = 1
#             if float(values['Cext']) < 0:
#                 sg.popup("Enter correct value of Cext")
#                 b = 1
#             if float(values['Ld']) <= 0:
#                 sg.popup("Enter correct value of Ld")
#                 b = 1
#
#             if float(values['Ls']) <= 0:
#                 sg.popup("Enter correct value of Ls")
#                 b = 1
#
#             if float(values['Cgd(ext)']) <= 0:
#                 sg.popup("Enter correct value of Cgd(ext)")
#                 b = 1
#
#             if float(values['Ron']) < 0:
#                 sg.popup("Enter correct value of Ron")
#                 b = 1
#
#             if float(values['Rgint']) <= 0:
#                 sg.popup("Enter correct value of Rgint")
#                 b = 1
#
#             if float(values['Ciss']) <= 0:
#                 sg.popup("Enter correct value of Ciss")
#                 b = 1
#             if v_gs1 == [] or id1_vds == []:
#                 #b is set to 3 , when all other data is present, but the dataset for output characteristics is not present.
#                 b = 3
#
#             if b == 3:
#                 sg.popup("Enter dataset for output characteristics")
#
#
#
#     # if all values are set properly ,b value is set to 2
#             if b == 2:
#                 r_gint = float(values['Rgint'])
#                 r_gext = float(values['Rgext'])
#                 c_iss = float(values['Ciss'])
#                 l_s = float(values['Ls'])
#
#                 if (r_gext + r_gint) <= (2 * l_s / c_iss) ** .5:
#                     sg.popup("gate circuit is underdamped, increase Rgext")
#                     b = 4
#                 else:
#                     plt.figure(1)

                    # -------------------- curvefit and data extraction section--------------------------------

i_0_max = 25
v_gs1 = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
id_vds = ['C:/proj/RediffMail.245151617893178/Id_Vds_C2M0080120D_Vgs_10V.csv',
          'C:/proj/RediffMail.245151617893178/Id_Vds_C2M0080120D_Vgs_12V.csv',
          'C:/proj/RediffMail.245151617893178/Id_Vds_C2M0080120D_Vgs_14V.csv',
          'C:/proj/RediffMail.245151617893178/Id_Vds_C2M0080120D_Vgs_16V.csv',
          'C:/proj/RediffMail.245151617893178/Id_Vds_C2M0080120D_Vgs_18V.csv',
          'C:/proj/RediffMail.245151617893178/Id_Vds_C2M0080120D_Vgs_20V.csv']
value1 = 'C:/proj/RediffMail.245151617893178/Transfer_characteristics.csv'
value2 = 'C:/proj/RediffMail.245151617893178/Reverse_transfer_characteristics.csv'
value3 = 'C:/proj/RediffMail.245151617893178/Output_capacitance.csv'
value4 = 'C:/proj/RediffMail.245151617893178/Diode_capacitance.csv'
k1b, k2b, m1b, k1a, k2a, m1a, k3a, vtd = reverse_transfer_capacitance(value2)
k4, k5, m2 = Outputcapacitance(value3)
vgs1, id_current1 = data_truncation_transfer_characteristics(value1, i_0_max)
vds3, id3, array1,id22,vds22 = data_preparation_output_characteristics(v_gs1, id_vds, i_0_max)
vth = 0
k_p = 0
theta_old = 0
theta_new = 0
for ig in range(0, 100):
    # Evaluate vth and k_p from the transfer characteristics
    vth, k_p = transfer_characteristics(
        vgs1, id_current1, theta_old)
    # Evaluate pvf,y,theta,kf from the output characteristics
    pvf, y, theta_new, kf, a_3d_list = output_characteristics(vds3, id3, array1, vth, k_p, theta_old)
    # Get difference between theta obtained from output ch and the theta that was fed into transfer ch
    del_theta = theta_new - theta_old

    if (abs(del_theta) < 1e-4):
        break
    else:
        theta_old = theta_new
thetaa = theta_new
#SSE calculation
# idnew1 = (0.5 * k_p * y * pvf / ((y - 1) * (1 + thetaa * (a_3d_list[:, 0] - vth)))) * (
#         (a_3d_list[:, 0] - vth) * a_3d_list[:, 1] - (pvf ** (y - 1) / y) * (
                    #             (a_3d_list[:, 0] - vth) ** (2 - y)) * (a_3d_list[:, 1] ** y))
                    # error = np.subtract(idnew1, id3)
                    # sq = np.power(error, 2)
                    # # aq = np.sum(sq)
                    # sum3 = 0
                    # for p in range(0, len(sq)):
                    #     sum3 = sum3 + sq[p]
                    #idnew is the current obtained using the parameters derived after 3D curvefitting
idnew = []
for re in range(len(vds22)):
    a = v_gs1[re] - vth
    b = a ** (2 - y)
    c = np.power(vds22[re], y)
    d = b * c
    e = ((pvf ** (y - 1)) / y) * d
    f = ((k_p * 0.5 * y * pvf) / (y - 1))
    g = 1 / (1 + thetaa * a)
    g1 = np.dot(a, vds22[re])
    h = np.subtract(g1, e)
    idnew.append(np.dot((g * f), h))
plt.subplot(2, 2, 3)
plt.plot(vgs1, id_current1, marker='o', ms=3, color='black', label="Id_Vgs data")
plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
             label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
plt.xlabel('Vgs(V)')
plt.ylabel('Id(A)')
plt.legend(loc='best', fancybox=True, shadow=True)
plt.grid(True)

plt.subplot(2, 2, 4)
for i in range(len(id22)):
    plt.plot(vds22[i], id22[i], marker='o', ms=3, color='black')
    plt.plot(vds22[i], idnew[i], marker='o', ms=3, color='red')
plt.xlabel('vds(V)')
plt.ylabel('Id(A)')
plt.grid(True)
# -------------------------------------Input parameters----------------------

v_dc = 800
i_0_min = 20
r_gint = 4.6
c_iss = .95
c_gd_ext = .015
c_ext = .47
# c_ext = .01

v_gg = 20
v_ee = -5
r_gext = 3
# r_gext = 15
# r_gext = 20
l_d = 6
l_s = 9
r_on = .08
ldc = 45


# -------------------------------------Analytical model section--------------
r_d = 0
i = 0
j = 0
i_0_step = float((i_0_max - i_0_min) / 4)
i_0 = []
for i in np.arange(i_0_min, i_0_max + 1, i_0_step):
    i_0.append(i)
i_0_len = len(i_0)
r_g = r_gint + r_gext

# t_off = []
# e_off = []
# dv_by_dt = []
# di_by_dt = []
# vds_ov = []


# -----------------------------Mode I-------------------------------

#-----------------------------Submode A-------------------------------

coxd = c_iss
c_gs = c_iss

n = 10**6
tau = r_g * c_iss * .001
a = l_s * (c_gs + coxd)
b = r_g * (c_gs + coxd)
c = 1
for ii in np.arange(20, 21, i_0_step):
    Vds1A = (r_on - r_d) * ii
    E1A = 0
    T1A = 0
    vgs = []
    vds = []
    vdpsp = []
    id = []
    ich = []
    vgs.append(v_gg)
    vds.append(Vds1A)
    vdpsp.append((r_on - r_d) * ii)
    vds.append(Vds1A)
    id.append(ii)
    Flag1 = 0
    i = 0
    print("vth=", vth, "pvf=", pvf,"thetaa=",thetaa,"y=",y ,"kf=",kf)
    while (i < n):
        if vgs[i] > vth and ((vgs[i] - vds[i]) > 0):
            ich.append(k_p * kf * (((vgs[i] - vth) * vds[i]) - (((pvf ** (y - 1)) / y) * (vgs[i] - vth) ** (2 - y) *
                                   vds[i] ** y)) / (1 + thetaa * (vgs[i] - vth)))
            cgd = coxd
            # print("ich",ich)
        if ((vgs[i] - vds[i]) <= 0) and (vds[i] <= (vgs[i] - vth) / pvf):
            break
        if ((vgs[i] - vds[i]) > 0) and (vds[i] > (vgs[i] - vth) / pvf):
            Flag1 = 1
            break
        if vds[i] < vtd:
            cds = (k4 / (1 + vds[i] / k5) ** (m2)) - (k1a / (((1 + (vds[i] / k2a)) ** m1a) + k3a))
        else:
            cds = (k4 / (1 + vds[i] / k5) ** (m2)) - (k1b / ((1 + ((vds[i] - vtd) / k2b))) ** m1b)

        vdpsp.append(vdpsp[i] + (tau * ((ii - id[i]) / c_ext)))
        id.append(id[i] + ((tau / (l_d + l_s)) * (vdpsp[i] - vds[i])))

        vgs.append(vgs[i] + tau * ((
                    v_ee + ((coxd * r_g / (coxd + cds))  * (ii - ich[i])) - vgs[i] - ((vdpsp[i] - vds[i]) * l_s / (
                        l_d + l_s))) / ((r_g * (c_gs + coxd)) - (r_g * (coxd ** 2) / (coxd + cds)))))
        vds.append(
            vds[i] + tau * ((ii - ich[i]) + (coxd / (r_g * (c_gs + coxd))) * (
                        v_ee - vgs[i] - ((vdpsp[i] - vds[i]) * (l_s / (l_d + l_s)))) / (
                                    (coxd + cds) - ((coxd ** 2) / (c_gs + coxd)))))

        E1A = E1A + vds[i] * ich[i] * tau
        i = i + 1
    # print("i=",i)
    # print("vds=",vds)
    # print("ich=", ich)
    Vgs1A = vgs[i]
    Vds1A = vds[i]
    Id1A = id[i]
    Vdpsp1A = vdpsp[i]
    T1A = (i - 1) * tau

    print("T1A=",T1A,"Vgs1A=", Vgs1A,"Vds1A=",Vds1A,"Vdpsp1A=",Vdpsp1A,"ich=",ich[i-1],"Id1A=",Id1A,"Io=",ii,"vd=",v_dc - Id1A* r_on,"E1A=",E1A)
    t = np.linspace(0, T1A,len(vgs))

    plt.figure(2)
    plt.subplot(2, 3, 1)
    # plt.subplot(3, 3, 1)
    plt.plot(t, vgs, marker='o', ms=1, color='black', label="Vgs vs t")
    plt.xlabel('t')
    plt.ylabel('vgs')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)

    t = np.linspace(0, T1A, len(vds))
    plt.subplot(2, 3, 2)
    plt.plot(t, vds, marker='o', ms=1, color='black', label="Vds vs t")
    plt.xlabel('t')
    plt.ylabel('vds')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)

    t = np.linspace(0, T1A, len(vdpsp))
    plt.subplot(2, 3, 3)
    plt.plot(t, vdpsp, marker='o', ms=1, color='black', label="vdpsp vs t")
    plt.xlabel('t')
    plt.ylabel('vdpsp')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)
    t = np.linspace(0, T1A, len(id))

    plt.subplot(2, 3, 4)
    plt.plot(t, id, marker='o', ms=1, color='black', label="id vs t")
    plt.xlabel('t')
    plt.ylabel('id')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)
    t = np.linspace(0, T1A, len(ich))

    plt.subplot(2, 3, 5)
    plt.plot(t, ich, marker='o', ms=1, color='black', label="ich vs t")
    plt.xlabel('t')
    plt.ylabel('ich')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)


    # ----------------------------------------------submode B-------------------
    ich = []
    vdpsp = []
    id =[]
    vds = []
    vgs = []
    vgs.append(Vgs1A)
    vds.append(Vds1A)
    vdpsp.append(Vdpsp1A)
    id.append(Id1A)
    vds1B = 0
    vdpsp1B =0
    Id1B = 0
    E1B = 0
    T1B =0
    j = 0
    # Flag1 = 0
    while (j < n):
        if (Flag1 == 0):
            ich.append(k_p * kf * (
                    (vgs[j] - vth) * vds[j] - (pvf ** (y - 1) / y) * (vgs[j] - vth) ** (2 - y) *
                    vds[j] ** y) / (1 + thetaa * (vgs[j] - vth)))
            cgd = k1a / ((1 + ((vds[j] - vgs[j]) / k2a)) ** m1a + k3a)

        if (Flag1 == 1):
            cgd = coxd
            if (vgs[j] > vth):
                ich.append(0.5 * k_p * (vgs[j] - vth) ** 2 / (1 + thetaa * (vgs[j] - vth)))
            else:
                ich.append(0)
        if ((vgs[j] - vds[j]) < 0 and (Flag1 == 1)):
            break
        if ((vds[j] > (vgs[j] - vth) / pvf) and (Flag1 == 0)):
            break

        if (vds[j] < vtd):
            cds = (k4 / (1 + vds[j] / k5) ** (m2)) - (k1a / (((1 + (vds[j] / k2a)) ** m1a) + k3a))
        else:
            cds = (k4 / (1 + vds[j] / k5) ** (m2)) - (k1b / ((1 + ((vds[j] - vtd) / k2b))) ** m1b)

        vdpsp.append(vdpsp[j] + tau * (ii - id[j]) / c_ext)
        id.append(id[j] + (tau / (l_d + l_s)) * (vdpsp[j] - vds[j]))


        vgs.append(vgs[j] + tau * ((v_ee + ((cgd * r_g / (cgd + cds)) * (ii - ich[j])) - vgs[j] - (
                                               (vdpsp[j] - vds[j]) * l_s / (
                                               l_d + l_s))) / (
                                               (r_g * (c_gs + cgd)) - (r_g * (cgd ** 2) / (cgd + cds)))))
        vds.append(
            vds[j] + tau * ((ii - ich[j]) + (cgd / (r_g * (c_gs + cgd))) * (
                    v_ee - vgs[j] - ((vdpsp[j] - vds[j]) * (l_s / (l_d + l_s)))) / (
                                    (cgd + cds) - ((cgd ** 2) / (c_gs + cgd)))))
        s= vds[j] * ich[j] * tau
        E1B = E1B + s
        j = j + 1
    Vgs1B = vgs[j]
    Vds1B = vds[j]
    Vdpsp1B =vdpsp[j]
    Id1B = id[j]
    T1B = (j - 1) * tau
    TI =  T1A + T1B
    EI = (E1A +  E1B) * 1e-3
    print(j, len(vgs), len(vds), len(id), len(vdpsp),len(ich))
    # print(T1B, Vgs1B, Vds1B, Vdpsp1B, ich[j-1], Id1B, ii,  E1B)
    print("T1B=", T1B, "Vgs1B=", Vgs1B, "Vds1B=", Vds1B, "Vdpsp1B=", Vdpsp1B, "ich=", ich[j], "Id1B=", Id1B, "Io=",
          ii, "vd=", v_dc - Id1B * r_on, "E1B=", E1B)
    print('T1B', T1B)
    print("E1B", E1B * 1e-3)
    print("TI", TI)
    print("EI", EI)
    t = np.linspace(0, T1B, len(vgs))
    print(j,len(vgs), len(vds), len(id), len(vdpsp))

    plt.figure(3)
    plt.subplot(2, 3, 1)
    plt.plot(t, vgs, marker='o', ms=3, color='black', label="Vgs vs t")

    plt.xlabel('t')
    plt.ylabel('vgs')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)

    t = np.linspace(0, T1B, len(vds))
    plt.subplot(2, 3, 2)
    plt.plot(t, vds, marker='o', ms=3, color='black', label="Vds vs t")

    plt.xlabel('t')
    plt.ylabel('vds')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)
    t = np.linspace(0, T1B, len(vdpsp))
    # print(len(vdpsp))
    plt.subplot(2, 3, 3)
    plt.plot(t, vdpsp, marker='o', ms=3, color='black', label="vdpsp vs t")

    plt.xlabel('t')
    plt.ylabel('vdpsp')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)
    t = np.linspace(0, T1B, len(id))
    plt.subplot(2, 3, 4)
    plt.plot(t, id, marker='o', ms=3, color='black', label="id vs t")

    plt.xlabel('t')
    plt.ylabel('id')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)
    t = np.linspace(0, T1B, len(ich))
    plt.subplot(2, 3, 5)
    plt.plot(t, ich, marker='o', ms=1, color='black', label="ich vs t")
    plt.xlabel('t')
    plt.ylabel('ich')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)

                    #     # -----------------------------Mode II-------------------------------
                    #
                    #     # ----------------------------------submode C--------------------------------------------

    vgs = []
    vds = []
    idc = []
    vd = []
    ich = []
    vgs.append(Vgs1B)
    vds.append(Vds1B)
    idc.append(ii)
    vd.append(v_dc)
    E2C = 0
    Flag2 = 0
    i = 0
    while i < n:
        ich.append((0.5 * k_p * (vgs[i] - vth) ** 2) / (1 + thetaa * (vgs[i] - vth)))

        cgd = k1a / ((1 + ((vds[i] - vgs[i]) / k2a)) ** m1a + k3a)
        coss = (k4 / (1 + vds[i] / k5) ** (m2))
        cd = (k4 / (1 + vd[i] / k5) ** (m2)) +c_ext + c_gd_ext

        time_1 = r_g * (c_gs + cgd) + (r_gext * c_gd_ext) + k_p*l_s*((vgs[i]-vth)- 3*thetaa*.5*((vgs[i] -vth)**2))
        time_2 = (r_g * cgd) + (r_gext * c_gd_ext)

        idc.append(idc[i] + tau * ((v_dc - vds[i] - vd[i]) / ldc))
        vd.append(vd[i] + tau * ((idc[i] - ii) / cd ))
        vds.append(vds[i] + tau * ((idc[i] - ich[i]) / (coss + c_gd_ext + c_ext)))
        vgs.append(vgs[i] + (tau/time_1) * (v_ee -vgs[i] + time_2*(idc[i] - ich[i])/(coss + c_gd_ext + c_ext)))
        E2C = E2C + (vds[i] * ich[i] * tau)
        if (vds[i] - vgs[i]) > vtd:
            break
        elif  ich[i] < (0.001*ii):
            Flag2=1
            break

        i = i + 1
    Vgs2C = vgs[i]
    Vds2C = vds[i]
    Vd2C = vd[i]
    Id2C = idc[i]
    T2C = (i - 1) * tau
    print("E2C=",E2C/1000,"T2C=",T2C,"Flag2=",Flag2)
    t = np.linspace(0, T2C, len(vgs))
    print(j, len(vgs), len(vds), len(idc), len(vdpsp))
    # print(T2C, Vgs2C, Vds2C, Vds2C, ich[i], Id1B, ii, Vd2C, E1B)
    print("T2C=", T2C, "Vgs2C=", Vgs2C, "Vds2C=", Vds2C, "Vdpsp2C=", Vds2C, "ich=", ich[i], "Id2C =", Id2C , "Io=",
          ii, "vd=", Vd2C, "E2C=", E2C)
    plt.figure(4)
    plt.subplot(2, 3, 1)
    plt.plot(t, vgs, marker='o', ms=3, color='black', label="Vgs vs t")
    # plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
    #          label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
    plt.xlabel('t')
    plt.ylabel('vgs')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)

    t = np.linspace(0, T2C, len(vds))
    plt.subplot(2, 3, 2)
    plt.plot(t, vds, marker='o', ms=3, color='black', label="Vds vs t")
    # plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
    #          label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
    plt.xlabel('t')
    plt.ylabel('vds')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)

    t = np.linspace(0, T2C, len(idc))
    plt.subplot(2, 3, 4)
    plt.plot(t, idc, marker='o', ms=3, color='black', label="idc vs t")
    # plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
    #          label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
    plt.xlabel('t')
    plt.ylabel('idc')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)
    t = np.linspace(0, T2C, len(ich))
    plt.subplot(2, 3, 5)
    plt.plot(t, ich, marker='o', ms=1, color='black', label="ich vs t")
    plt.xlabel('t')
    plt.ylabel('ich')
    plt.legend(loc='best', fancybox=True, shadow=True)
    plt.grid(True)
    # plt.show()


                    #
                    #     #  --------------------------------------Submode D-------------------------------------------------
    vgs = []
    vds = []
    idc = []
    vd = []
    ich = []
    E2D = 0
    Flag3 = 0
    vgs.append(Vgs2C)
    vds.append(Vds2C)
    idc.append(Id2C)
    vd.append(Vd2C)
    if (Flag2 == 1):
        Vgs2D = Vgs2C
        Vds2D = Vds2C
        VD2D = Vd2C
        Id2D = Id2C
        T2D = 0
        break
    else:
        i = 0
        while i <n:
            if (vgs[i] > vth):
                ich.append((0.5 * k_p * ((vgs[i] - vth) ** 2)) / (1 + thetaa * (vgs[i] - vth)))
                cgd = k1b / ((1 + ((vds[i] - vgs[i] - vtd) / k2b)) ** m1b)
                coss = k4 / ((1 + (vds[i] / k5)) ** m2)
                cd = (k4 / (1 + (vd[i] / k5)) ** m2) + c_ext + c_gd_ext
                time_1 = r_g * (c_gs + cgd) + (r_gext * c_gd_ext) + k_p * l_s * (
                            (vgs[i] - vth) - 3 * thetaa * .5 * ((vgs[i] - vth) ** 2))
                time_2 = (r_g * cgd) + (r_gext * c_gd_ext)
                vd.append(vd[i] + (tau * ((idc[i] - ii) / cd)))
                idc.append(idc[i] + (tau * ((v_dc - vds[i] - vd[i]) / ldc)))

                vds.append(vds[i] + (tau * ((idc[i] - ich[i]) / (coss + c_gd_ext + c_ext))))
                vgs.append(vgs[i] + (tau/time_1) * (v_ee -vgs[i] + time_2*(idc[i] - ich[i])/(coss + c_gd_ext + c_ext)))

                E2D = E2D + (vds[i] * ich[i] * tau)
                if ich[i] < (0.001 * ii):
                    break
                elif vd[i]<0:
                    Flag3 = 1
                    break

            i = i + 1
        Vgs2D = vgs[i]
        print("ich[i]=", ich[i], 'Flag3=', Flag3)
        Vds2D = vds[i]
        Vd2D = vd[i]
        Id2D = idc[i]
        T2D = (i - 1) * tau
        TII = T2C + T2D
        EII = (E2C + E2D) * .001
        print('EII = ', EII)
        print('Ich2A2 = ', ich[i - 1])
        print('TII = ', TII)
        print('EII = ', EII)
        t = np.linspace(0, T2D, len(vgs))
        print(i, len(vgs), len(vds), len(idc), len(vdpsp))
        # print(T2D, Vgs2D, Vds2D, Vds2D, ich[i], Id2D, ii, Vd2D, E2D)
        print("T2D=", T2D, "Vgs2D=", Vgs2D, "Vds2D=", Vds2D, "Vdpsp2D=", Vds2D, "ich=", ich[i], "Id2D =", Id2D, "Io=",
              ii, "vd=", Vd2D, "E2D=", E2D)
        plt.figure(5)
        plt.subplot(2, 3, 1)
        plt.plot(t, vgs, marker='o', ms=3, color='black', label="Vgs vs t")
        # plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
        #          label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
        plt.xlabel('t')
        plt.ylabel('vgs')
        plt.legend(loc='best', fancybox=True, shadow=True)
        plt.grid(True)

        t = np.linspace(0, T2D, len(vds))
        plt.subplot(2, 3, 2)
        plt.plot(t, vds, marker='o', ms=3, color='black', label="Vds vs t")
        # plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
        #          label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
        plt.xlabel('t')
        plt.ylabel('vds')
        plt.legend(loc='best', fancybox=True, shadow=True)
        plt.grid(True)
        # t = np.linspace(0, T1B, len(vdpsp))
        # # print(len(vdpsp))
        # plt.subplot(2, 3, 3)
        # plt.plot(t, vdpsp, marker='o', ms=3, color='black', label="vdpsp vs t")
        # # plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
        # #          label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
        # plt.xlabel('t')
        # plt.ylabel('vdpsp')
        # plt.legend(loc='best', fancybox=True, shadow=True)
        # plt.grid(True)
        t = np.linspace(0, T2D, len(idc))
        plt.subplot(2, 3, 4)
        plt.plot(t, idc, marker='o', ms=3, color='black', label="idc vs t")
        # plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
        #          label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
        plt.xlabel('t')
        plt.ylabel('idc')
        plt.legend(loc='best', fancybox=True, shadow=True)
        plt.grid(True)
        t = np.linspace(0,T2D, len(ich))
        plt.subplot(2, 3, 5)
        plt.plot(t, ich, marker='o', ms=1, color='black', label="ich vs t")
        plt.xlabel('t')
        plt.ylabel('ich')
        plt.legend(loc='best', fancybox=True, shadow=True)
        plt.grid(True)
        plt.show()
        # -----------------Mode 3-------------------------
        def eqn_cq1(vds):
            vd = v_dc - vds
            coss = k4 / ((1 + (vds / k5)) ** m2)
            cd = (k4 / (1 + (vd/ k5)) ** m2)
            cteq = cd + c_ext + c_gd_ext
            cbeq = coss * (c_ext + c_gd_ext)
            return ((cteq * cbeq)/(cteq + cbeq))


        def eqn_cq2(vds):
            vd = v_dc - vds
            coss = k4 / ((1 + (vds / k5)) ** m2)
            cd = (k4 / (1 + (vd / k5)) ** m2)
            cteq = cd + c_ext + c_gd_ext
            cbeq = coss * (c_ext + c_gd_ext)
            return (cteq + cbeq)

        #
        # for i in range(i_0_len):

        global cbeqconj, vdpsp3
        int_I = []
        int_II = []
        if Flag3==0:
            sum1 = 0
            # vds = np.linspace(0, v_dc, v_dc*10)
            # vds = np.linspace(0, v_dc, (int)(v_dc /.1))

            # vd = v_dc -  vds

            # coss = k4 / ((1 + (vds[i] / k5)) ** m2)
            # cd = (k4 / (1 + (vd[i] / k5)) ** m2)

            # int1 = (cteq * cbeq)/((cteq + cbeq)*v_dc)
            vds = list(np.arange(0, v_dc, 0.1).tolist())

            int1, err1 = quad(eqn_cq1, 0, vds[len(vds) - 1])
            int_I.append(int1)

            int2, err2 = quad(eqn_cq2, 0, vds[len(vds) - 1])
            int_II.append(int2)

            cq1 = int1/v_dc
            cq2 = int2 / v_dc


            # cq1,er1 = integrate.quad(int1, 0, v_dc)
            # int2 = (cteq + cbeq)/v_dc
            # int2 =np.dot(1/v_dc,(np.add(cteq,cbeq)))
            # print("vds = ",vds)
            # print("vd = ", vd)
            # print("coss = ", coss)
            # print("cd = ", cd)
            # print("cteq = ", cteq)
            # print("cbeq = ", cbeq)
            # print("int1 = ", int1)
            # print("int2 = ", int2)
            # print("sum1 ", sum1)
            # cq2 ,er2= integrate.quad(int2, 0, v_dc)
            wo = 1/((ldc * cq1)**.5)
            Ao= Vds2D
            A1 = ii/cq2
            A2 = (Vds2D - A1)/wo
            cbeqconj = (k4 / ((1 + (vds[i] / k5)) ** m2)) + c_ext + c_gd_ext

            A3 = ldc* cbeqconj* (wo**2) - 1
            TIII = (v_dc -Ao)/A1
            vdpsp3 = Ao + A1*TIII + A2*math.sin(wo * TIII)
            Idc3 = cbeqconj* (A1 + A2*wo*math.cos(wo * TIII))
        else:
            TIII = 0
            vdpsp3 = 0
            Idc3 = 0
            break
        print("int1=",int1,"int2=",int2)
        #t3 (ns)	Vd's'3(V)	I'd3(A)	dv/dt(V/ns)

        print("t3=",TIII,"vdpsp3=", vdpsp3,"Idc3=",Idc3,)

        #------------------------ModeIv---------------------
        if Flag3 == 0:
            ldceq = ldc + l_d + l_s
            w1 = 1 / ((ldceq * cbeqconj) ** .5)
            A4 = (((vdpsp3 - v_dc) ** 2) + ((ldceq / cbeqconj) * (Idc3 ** 2))) ** .5
            phi = math.atan((((vdpsp3 - v_dc) / Idc3) * ((cbeqconj / ldceq) ** .5)))
            TIV = (1.57 - phi) / w1
        else:
            TIV = 0
            Vdsmax = 0
            break
        #----------Final output----

        Toff = TI + TII + TII + TIV
        EOFF = EI + EII
        dvdt = (vdpsp3 - Vds2D) / TIII
        didt = id3 / TIV
        vdsmax = v_dc + A4
        print("dvdt=",dvdt)
# window_first.close()


