import PySimpleGUI as sg
from curvefit_code_3c import *
from gui_second import *


layout = [[sg.Text('Turn off Switching Loss Calculator for SiC MOSFET and Schottky Diode Pair in a Half Bridge Configuration', auto_size_text=False,
                   justification='center', font=("Arial", 20), size=(100, 1))],
          [sg.Text('OPERATING CONDITION', font=("Arial", 16), justification='center')],
          [sg.Text('DC bus voltage (Vdc)', font=(30), size=(20, 1)),
           sg.InputText(font=(30), do_not_clear=True,key = 'Vdc',justification='left', size=(5, 1)),
           sg.Text('(V)', font=(30), size=(3, 1)), sg.Text('', size=(2, 1)),
           sg.Text('Minimum Load current (Iomin)', font=(30), size=(30, 1)),
           sg.InputText(font=(30), do_not_clear=True,key='Iomin',justification='left', size=(5, 1)),
           sg.Text('(A)', font=(30), size=(3, 1)), sg.Text('', size=(2, 1)),
           sg.Text('Maximum Load current (Iomax)', font=(30), size=(30, 1)),
           sg.InputText(font=(30), do_not_clear=True,key='Iomax' ,justification='left', size=(5, 1)),
           sg.Text('(A)', font=(30), size=(3, 1))],
          # [sg.Text('DC bus voltage', font=(30), size=(13, 1)), sg.InputText(font=(30), do_not_clear=False, justification='left', size=(5,1)), sg.Text('(V)', font=(30), size=(3, 1)), sg.Text('', size=(8,1)), sg.Text('Load current', font=(30), size=(12, 1)), sg.InputText(font=(30), do_not_clear=False, justification='left', size=(5,1)), sg.Text('(A)', font=(30), size=(3, 1))],
          [sg.Text('GATE DRIVER PARAMETERS', font=("Arial", 16), justification='center')],
          [sg.Text('Vgg', font=(30), size=(5, 1)), sg.InputText(font=(30), do_not_clear=True,key ='Vgg', justification='left', size=(5,1)), sg.Text('(V)', font=(30),
           size=(5, 1)), sg.Text('', size=(2,1)), sg.Text('Vee', font=(30), size=(5, 1)), sg.InputText(font=(30), do_not_clear=True,key ='Vee', justification='left', size=(5,1)),
           sg.Text('(V)', font=(30), size=(3, 1)), sg.Text('', size=(2,1)), sg.Text('Rgext', font=(30), size=(7, 1)),

           sg.InputText(font=(30), do_not_clear=True,key ='Rgext', justification='left', size=(5, 1)),
           sg.Text('(Ω)', font=(30), size=(5, 1)),
           sg.Text('Tf', font=(30), size=(5, 1)),
           sg.InputText(font=(30), do_not_clear=True,key ='Tf' ,justification='left', size=(5, 1)),
           sg.Text('(ns)', font=(30), size=(5, 1)), sg.Text('', size=(2, 1))
           ],
          [sg.Text('EXTERNAL CIRCUIT PARASITICS', font=("Arial", 16), justification='center')],
          [sg.Text('Ld', font=(30), size=(12, 1)),
           sg.InputText(font=(30), do_not_clear=True,key='Ld', justification='left', size=(5, 1)),
           sg.Text('(nH)', font=(30), size=(5, 1)), sg.Text('', size=(5, 1)), sg.Text('Ls', font=(30), size=(12, 1)),
           sg.InputText(font=(30), do_not_clear=True,key='Ls', justification='left', size=(5, 1)),
           sg.Text('(nH)', font=(30), size=(8, 1)), sg.Text('Cgd(ext)', font=(30), size=(12, 1)),
           sg.InputText(font=(30), do_not_clear=True,key ='Cgd(ext)', justification='left', size=(5, 1)),
           sg.Text('(pF)', font=(30), size=(4, 1)), sg.Text('', size=(2, 1)),
           sg.Text('Cak(ext)', font=(30), size=(12, 1)),
           sg.InputText(font=(30), do_not_clear=True,key='Cak(ext)', justification='left', size=(5, 1)),
           sg.Text('(pF)', font=(30), size=(4, 1)), sg.Text('', size=(2, 1)), ],
          [sg.Text('SiC MOSFET PARAMETERS', font=("Arial", 16), justification='center')],
          [sg.Text('Ron', font=(30), size=(6, 1)),
           sg.InputText(font=(30), do_not_clear=True,key ='Ron', justification='left', size=(5, 1)),
           sg.Text('(Ω)', font=(30), size=(5, 1)), sg.Text('', size=(2, 1)),
           sg.Text('Rgint', font=(30), size=(6, 1)),
           sg.InputText(font=(30), do_not_clear=True,key = 'Rgint', justification='left', size=(5, 1)),
           sg.Text('(Ω)', font=(30), size=(5, 1)), sg.Text('', size=(2, 1)),
           sg.Text('Ciss', font=(30), size=(6, 1)),
           sg.InputText(font=(30), do_not_clear=True,key = 'Ciss', justification='left', size=(5, 1)),
           sg.Text('(nF)', font=(30), size=(3, 1)), sg.Text('', size=(2, 1)),],
          [sg.Text('Dataset for output characteristics', font=(30), justification='center'),sg.Button('Insert',font=(30))],
          [sg.Text('Dataset for Transfer Characteristics (Id(A)-Vgs(V)) (.csv):', size=(57, 1), font=(30)),
           sg.Input(key='file for  transfer characteristics',do_not_clear=True), sg.FileBrowse(font=(30))],
          [sg.Text('Dataset for Reverse Transfer capacitance (Crss(pF)-Vds(V)) (.csv):', size=(57, 1), font=(30)),
           sg.Input(key='file for reverse transfer characteristics'), sg.FileBrowse(font=(30))],
          [sg.Text('Dataset for Output Capacitance (Coss(pF)-Vds(V)) (.csv):', size=(57, 1), font=(30)),
           sg.Input(key='file for Output Capacitance',do_not_clear=True), sg.FileBrowse(font=(30))],
            [sg.Text('SiC SCHOTTKY DIODE PARAMETERS', font=("Arial", 16), justification='center')],
          [sg.Text('Dataset for Diode Capacitance (CD(pF)-VD(V)) (.csv):', size=(57, 1), font=(30)),
           sg.Input(key='file for Diode Capacitance',do_not_clear=True), sg.FileBrowse(font=(30))],
          [sg.Text('', font=(30), size=(20, 1)), sg.Submit(font=(30)), sg.Text('', font=(30), size=(10, 1)),
           sg.Cancel(font=(30))]]
window_first = sg.Window('Calculator', layout)



#---Start of the main function----
count = 0
while True :
    #count is set to 1 , when user enters the data in the main window fields for the first time
    count = count + 1
    # when all field values are correct b is set to 2
    b =2
    a = 0
    event, values = window_first.Read()
    if event == "Cancel" or event == sg.WIN_CLOSED:
        break
    global v_gs1
    global id1_vds
    global flagb
    #---if the user is entering the values in the main window for the first time , the dataset values are empty,so flagb is set to 1
    #---Then, if user tries to submit with those empty values, the popup will come showing that outdataset characteristics should not be empty
    if count ==1:
        v_gs1 = []
        id1_vds = []
        flagb =1

    if  event == 'Insert':
        v_gs12,id1_vds2 = second_window()
        v_gs1 = []
        id1_vds = []
        for t1 in range(len(v_gs12)):
            if v_gs12[t1] != '':
                v_gs1.append(float(v_gs12[t1]))
        for t12 in range(len(id1_vds2)):
            if id1_vds2[t12] != '':
                id1_vds.append(id1_vds2[t12])
        # when the v_gs and id1_vds values are empty ,flagb value is set to 1
        if v_gs1 == [] or id1_vds == []:
            flagb = 1
        else:
            flagb = 0
    if  event=='Submit':
        b1 = ''
        # c has the count of values that are not provided by the user
        c = 0
        for v, l in values.items():
            if values[v] == '':
                c = c + 1
                b1 = b1 + v + ','
        if b1 != '':
            te = 0
            while (b1.rfind('Browse') != -1):
                #te has the count of the number of times the string 'Browse' appears in a string
                te = te + 1
                p = b1.rfind('Browse')
                b2 = str(b1[p + 8:])

                if (te == 4):
                    b1 = str(b1[:p]) + 'f' + b2
                else:
                    b1 = str(b1[:p]) + b2
            tr = b1.rfind(',')
            #b is set to 0 when one or more fields are empty
            b=0
            if b==0 and flagb == 1:
                sg.popup(b1[:tr] + " and dataset for output characteristics  are missing.")
            if b == 0 and c >1 and flagb!=1:
                sg.popup(b1[:tr] + " are missing.")
            if b == 0 and c ==1 and flagb!=1:
                sg.popup(b1[:tr] + " is missing.")
            # if some values are missing in the gui , then b is set to 0

        else:
            if float(values['Vdc']) <= 0:
                sg.popup("Enter correct value of Vdc")
                b = 1
                # if some values are incorrect in the gui , then b is set to 1

            if float(values['Iomin']) <= 0:
                sg.popup("Enter correct value of Iomin")
                # if value of Iomin is  incorrect in the gui , then a is set to 1
                a = 1
                b = 1

            if float(values['Iomax']) <= 0:
                sg.popup("Enter correct value of Iomax")
                # if value of Iomax is  incorrect in the gui , then a is set to 1
                b = 1
                a = 1

            if float(values['Iomin']) >= float(values['Iomax']) and a!=1:
                sg.popup("Enter correct values of Iomin & Iomax")
                b = 1

            if float(values['Vgg']) <= 0:
                sg.popup("Enter correct value of Vgg")
                # if value of Vgg is  incorrect in the gui , then q is set to 1
                b = 1
                q = 1

            if float(values['Vee']) > 0:
                sg.popup("Enter correct value of Vee")
                # if value of Vee is  incorrect in the gui , then q is set to 1
                b = 1
                q = 1

            if float(values['Vee']) >= float(values['Vgg']) and q!=1:
                sg.popup("Enter correct values of Vee & Vgg")
                b = 1

            if float(values['Rgext']) <= 0:
                sg.popup("Enter correct value of Rgext")
                b = 1
            if float(values['Tf']) < 0:
                sg.popup("Enter correct value of Tf")
                b = 1
            if float(values['Ld']) <= 0:
                sg.popup("Enter correct value of Ld")
                b = 1

            if float(values['Ls']) <= 0:
                sg.popup("Enter correct value of Ls")
                b = 1

            if float(values['Cgd(ext)']) <= 0:
                sg.popup("Enter correct value of Cgd(ext)")
                b = 1

            if float(values['Cak(ext)']) <= 0:
                sg.popup("Enter correct value of Cak(ext)")
                b = 1

            if float(values['Ron']) < 0:
                sg.popup("Enter correct value of Ron")
                b = 1

            if float(values['Rgint']) <= 0:
                sg.popup("Enter correct value of Rgint")
                b = 1

            if float(values['Ciss']) <= 0:
                sg.popup("Enter correct value of Ciss")
                b = 1
            if v_gs1 == [] or id1_vds == []:
                #b is set to 3 , when all other data is present, but the dataset for output characteristics is not present.
                b = 3

            if b == 3:
                sg.popup("Enter dataset for output characteristics")



    # if all values are set properly ,b value is set to 2
            if b == 2:
                r_gint = float(values['Rgint'])
                r_gext = float(values['Rgext'])
                c_iss = float(values['Ciss'])
                l_s = float(values['Ls'])

                if (r_gext + r_gint) <= (2 * l_s / c_iss) ** .5:
                    sg.popup("gate circuit is underdamped, increase Rgext")
                    b = 4
                else:
                    plt.figure(1)

                    # -------------------- Curve fit and data extraction section--------------------------------

                    i_0_max = float(values['Iomax'])
                    k1b, k2b, m1b, k1a, k2a, m1a, k3a, vtd = reverse_transfer_capacitance(
                        values['file for reverse transfer characteristics'])
                    k4, k5, m2 = Outputcapacitance(values['file for Output Capacitance'])
                    k6, k7, m3 = diode_capacitance(values['file for Diode Capacitance'])

                    vgs1, id_current1 = data_truncation_transfer_characteristics(
                        values['file for  transfer characteristics'], i_0_max)
                    vds3, id3, array1,id22,vds22 = data_preparation_output_characteristics(v_gs1, id1_vds, i_0_max)

                    vth = 0
                    k_p = 0
                    theta_old = 0
                    theta_new = 0
                    for ig in range(0, 100):
                        # Evaluate vth and k_p from the transfer characteristics
                        vth, k_p = transfer_characteristics(
                            vgs1, id_current1, theta_old)
                        # Evaluate pvf,y,theta,kf from the output characteristics
                        pvf, y, theta_new, kf ,a_3d_list= output_characteristics(vds3, id3,array1,vth, k_p,theta_old)
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

                    plt.subplot(2, 3, 4)
                    plt.plot(vgs1, id_current1, marker='o', ms=3, color='black', label="Id_Vgs data")
                    plt.plot(vgs1, func_id_vgs(np.transpose(vgs1), vth, k_p),
                             label="id=(kp*((vgs -vth_th)**2))/(2 * (1 + theta*(vgs -vth_th)))")
                    plt.xlabel('Vgs(V)')
                    plt.ylabel('Id(A)')
                    plt.legend(loc='best', fancybox=True, shadow=True)
                    plt.grid(True)

                    plt.subplot(2, 3, 5)
                    for i in range(len(id22)):
                        plt.plot(vds22[i],id22[i], marker='o', ms=3, color='black')
                        plt.plot(vds22[i], idnew[i], marker='o', ms=3, color='red')
                    plt.xlabel('vds(V)')
                    plt.ylabel('Id(A)')
                    plt.grid(True)

                    # -------------------- Input parameters--------------------------------

                    v_dc = float(values['Vdc'])
                    i_0_min = float(values['Iomin'])
                    i_0_max = float(values['Iomax'])
                    # r_gint = float(values['Rgint'])
                    # c_iss = float(values['Ciss'])
                    c_gd_ext = (float(values['Cgd(ext)'])) * 0.001
                    c_ak_ext = (float(values['Cak(ext)'])) * 0.001
                    v_gg = float(values['Vgg'])
                    v_ee = float(values['Vee'])
                    # r_gext = float(values['Rgext'])
                    tf = float(values['Tf'])
                    l_d = float(values['Ld'])
                    # l_s = float(values['Ls'])
                    r_on = float(values['Ron'])


                    # ------------------------Analytical model section--------------

                    r_d = 0
                    i = 0
                    j = 0
                    i_0_step = float((i_0_max - i_0_min) / 4)
                    i_0 = []
                    for i in np.arange(i_0_min, i_0_max + 1, i_0_step):
                        i_0.append(i)
                    i_0_len = len(i_0)
                    r_g = r_gint + r_gext

                    t_off = []
                    e_off = []
                    dv_by_dt = []
                    di_by_dt = []
                    vds_ov = []

                    # # -----------------------------Mode I-------------------------------

                    # # -----------------------------Submode A-------------------------------
                    #
                    coxd = c_iss
                    c_gs = c_iss
                    n = r_g * c_iss * 1000
                    tau = r_g * c_iss * .001
                    sampl = (tf - 0) / tau
                    t = np.linspace(0, tf, (int)(sampl - 1))
                    a = l_s * (c_gs + coxd)
                    b = r_g * (c_gs + coxd)
                    c = 1

                    # If a is 0, then incorrect equation
                    if a == 0:
                        print("Input correct quadratic equation")
                    else:
                        p1, p2 = equationroots(a, b, c)

                    sum = 0
                    vgs = []

                    # Doing matrix operations step by step
                    subterm1a_2 = np.divide(np.power(2.718, np.dot(p1, t)), (p1 ** 2))
                    subterm2a_2 = np.divide(np.power(2.718, np.dot(p2, t)), (p2 ** 2))
                    prod = (p1 * p2) / (p2 - p1)
                    a11 = p1 + p2
                    a22 = p1 * p2
                    v_gs = np.add(v_gg, np.dot(((v_ee - v_gg) / tf), np.add(
                        np.subtract(t, np.dot(prod, np.subtract(subterm1a_2, subterm2a_2))), (a11 / a22))))

                    Vgs1A = v_gs[len(v_gs) - 1]

                    for ii in np.arange(i_0_min, i_0_max + 1, i_0_step):
                        Vds1A = (r_on - r_d) * ii

                        ich = []
                        sum = 0
                        for ir in range(0, len(v_gs)):
                            ich.append(k_p * kf * ((v_gs[ir] - vth) * Vds1A - ((pvf ** (y - 1)) / y) * (
                                        (v_gs[ir] - vth) ** (2 - y)) * (
                                                           Vds1A ** y)) / (1 + thetaa * (v_gs[ir] - vth)))
                            sum = sum + ich[ir]

                        E1A = (sum * Vds1A * tau)

                        # -----------------------------------submode B1---------------

                        vgs = []
                        vds = []
                        ich = []
                        vds1 = []
                        vgs1 = []
                        Vgs1B1 = 0
                        Vds1B1 = 0
                        E1B1 = 0
                        vgs.append(Vgs1A)

                        vds.append(Vds1A)
                        Flag1 = 0
                        xflag = 0
                        ichflag = 0

                        i = 0

                        while (i < n):
                            if ((vgs[i] - vds[i]) <= 0):
                                xflag = 1
                                break
                            if (vds[i] > (vgs[i] - vth) / pvf):
                                Flag1 = 1
                                break

                            ich.append(k_p * kf * (
                                    (vgs[i] - vth) * vds[i] - (pvf ** (y - 1) / y) * (vgs[i] - vth) ** (2 - y) * vds[
                                i] ** y) / (
                                               1 + thetaa * (vgs[i] - vth)))
                            cgd = coxd

                            if vds[i] < vtd:
                                cds = (k4 / (1 + vds[i] / k5) ** (m2)) + (k6/((1 + (vds[i]/k7)) ** m3)) - (k1a / (((1 + (vds[i] / k2a)) ** m1a) + k3a))
                            else:
                                cds = (k4 / (1 + vds[i] / k5) ** (m2)) + (k6/((1 + (vds[i]/k7)) ** m3)) - (k1b / ((1 + ((vds[i] - vtd) / k2b))) ** m1b)

                            vgs.append(vgs[i] + tau * (v_ee + (cgd / (cgd + cds)) * r_g * (ii - ich[i]) - vgs[i]) / (
                                    r_g * (c_iss + cgd) - r_g * (cgd ** 2) / (cgd + cds)))
                            vds.append(
                                vds[i] + tau * ((ii - ich[i]) + ((cgd / (c_iss + cgd)) * ((v_ee - vgs[i]) / r_g))) / (
                                        (cgd + cds) - (cgd ** 2) / (c_iss + cgd)))

                            E1B1 = E1B1 + vds[i] * ich[i] * tau
                            i = i + 1

                        Vgs1B1 = vgs[i]
                        Vds1B1 = vds[i]
                        T1B1 = (i - 1) * tau

                        # ----------------------------------------------submode B2-------------------

                        vgs = []
                        vds = []
                        ich = []
                        vgs.append(Vgs1B1)
                        vds.append(Vds1B1)
                        vgs1B2 = 0
                        vds1B2 = 0
                        E1B2 = 0
                        i = 0
                        while (i < n):
                            if ((vgs[i] - vds[i]) < 0 and (Flag1 == 1)):
                                break
                            if ((vds[i] > (vgs[i] - vth) / pvf) and (Flag1 == 0)):
                                break

                            if (Flag1 == 0):
                                ich.append(k_p * kf * (
                                        (vgs[i] - vth) * vds[i] - (pvf ** (y - 1) / y) * (vgs[i] - vth) ** (2 - y) *
                                        vds[
                                            i] ** y) / (1 + thetaa * (vgs[i] - vth)))
                                cgd = k1a / ((1 + ((vds[i] - vgs[i]) / k2a)) ** m1a + k3a)

                            if (Flag1 == 1):
                                if (vgs[i] > vth):
                                    ich.append(0.5 * k_p * (vgs[i] - vth) ** 2 / (1 + thetaa * (vgs[i] - vth)))
                                else:
                                    ich.append(0)
                                cgd = coxd

                            if (vds[i] < vtd):
                                cds = (k4 / (1 + vds[i] / k5) ** (m2)) + (k6/((1 + (vds[i]/k7)) ** m3)) - (k1a / (((1 + (vds[i] / k2a)) ** m1a) + k3a))
                            else:
                                cds = (k4 / (1 + vds[i] / k5) ** (m2)) + (k6/((1 + (vds[i]/k7)) ** m3)) - (k1b / ((1 + ((vds[i] - vtd) / k2b))) ** m1b)

                            vgs.append(vgs[i] + tau * (v_ee + ((cgd / (cgd + cds)) * r_g * (ii - ich[i])) - vgs[i]) / (
                                    r_g * (c_iss + cgd) - (r_g * (cgd ** 2) / (cgd + cds))))
                            vds.append(
                                vds[i] + tau * ((ii - ich[i]) + ((cgd / (c_iss + cgd)) * ((v_ee - vgs[i]) / r_g))) / (
                                        (cgd + cds) - ((cgd ** 2) / (c_iss + cgd))))

                            E1B2 = E1B2 + (vds[i] * ich[i] * tau)
                            i = i + 1

                        ti = np.linspace(0, tau * (len(ich) - 1), len(ich))
                        # print("ti =", ti)
                        Vgs1B2 = vgs[i]
                        Vds1B2 = vds[i]
                        T1B2 = (i - 1) * tau
                        TI = tf + T1B1 + T1B2
                        EI = (E1A + E1B1 + E1B2) * 1e-3

                        print("T1B2", T1B2)

                        print("E1B2", E1B2 * 1e-3)
                        print("TI", TI)
                        print("EI", EI)

                        # ----------------------------------Mode II--------------------------------------------

                        # ----------------------------------Submode A1--------------------------------------------

                        vgs = []
                        vds = []
                        id = []
                        vd = []
                        ich = []
                        vgs.append(Vgs1B2)
                        vds.append(Vds1B2)
                        id.append(ii)
                        vd.append(v_dc)
                        vgs2A1 = 0
                        vds2A1 = 0
                        vd2A1 = 0
                        id2A1 = 0
                        t2A1 = 0
                        E2A1 = 0
                        i = 0
                        while i < n:
                            if (vds[i] - vgs[i]) > vtd:
                                break

                            if (vgs[i] > vth):
                                ich.append(0.5 * k_p * (vgs[i] - vth) ** 2 / (1 + thetaa * (vgs[i] - vth)))
                            else:
                                ich.append(0)
                            cgd = k1a / ((1 + ((vds[i] - vgs[i]) / k2a)) ** m1a + k3a)
                            coss = (k4 / (1 + vds[i] / k5) ** (m2)) + (k6/((1 + (vds[i] / k7)) ** m3))
                            cd = (k4 / (1 + (vd[i] / k5)) ** m2) + (k6/((1 + (vd[i] / k7)) ** m3))

                            time_1 = r_g * (c_iss + cgd) + (r_gext * c_gd_ext)
                            time_2 = r_g * cgd + (r_gext * c_gd_ext)

                            A = (v_ee - ((v_dc * l_s) / (l_s + l_d))) / time_1
                            B = (l_s / ((l_s + l_d) * time_1))
                            C = time_2 / (time_1 * (coss + c_gd_ext))
                            D = -1 / time_1

                            id.append(id[i] + tau * ((v_dc - vds[i] - vd[i]) / (l_d + l_s)))
                            vd.append(vd[i] + tau * ((id[i] - ii) / (cd + c_ak_ext)))
                            vds.append(vds[i] + tau * ((id[i] - ich[i]) / (coss + c_gd_ext)))
                            vgs.append(vgs[i] + tau * (A + B * (vds[i] + vd[i]) + C * (id[i] - ich[i]) + D * vgs[i]))
                            E2A1 = E2A1 + (vds[i] * ich[i] * tau)
                            E2A1new = E2A1 / 1000
                            i = i + 1

                        vgs2A1 = vgs[i]
                        vds2A1 = vds[i]
                        vd2A1 = vd[i]
                        id2A1 = id[i]
                        t2A1 = (i - 1) * tau

                        #  --------------------------------------submode A2-------------------------------------------------

                        print("*****----start  of submode A2-----***********")

                        vgs = []
                        vds = []
                        id = []
                        vd = []
                        ich = []
                        vgs.append(vgs2A1)
                        vds.append(vds2A1)
                        id.append(id2A1)
                        vd.append(vd2A1)
                        vgs2A2 = 0
                        vds2A2 = 0
                        vd2A2 = 0
                        id2A2 = 0
                        t2A2 = 0
                        E2A2 = 0
                        i = 0
                        # tau1 = .001
                        while i < n:
                            if (vd[i] < 0):
                                break

                            if (vgs[i] > vth):
                                Flagich = 0
                                ich.append(0.5 * k_p * (vgs[i] - vth) ** 2 / (1 + thetaa * (vgs[i] - vth)))
                            else:
                                Flagich = 1
                                ich.append(0)

                            cgd = k1b / ((1 + ((vds[i] - vgs[i] - vtd) / k2b)) ** m1b)
                            coss = (k4 / (1 + (vds[i] / k5)) ** (m2)) + (k6 / ((1 + (vds[i] / k7)) ** m3))
                            cd = (k4 / (1 + (vd[i] / k5)) ** (m2)) + (k6 / ((1 + (vd[i] / k7)) ** m3))

                            # print("value of vd", vd[i])

                            time_1 = r_g * (c_iss + cgd) + (r_gext * c_gd_ext)
                            time_2 = r_g * cgd + (r_gext * c_gd_ext)

                            A = (v_ee - ((v_dc * l_s) / (l_s + l_d))) / time_1
                            B = (l_s / ((l_s + l_d) * time_1))
                            C = time_2 / (time_1 * (coss + c_gd_ext))
                            D = -1 / time_1

                            id.append(id[i] + tau * ((v_dc - vds[i] - vd[i]) / (l_d + l_s)))
                            vd.append(vd[i] + tau * ((id[i] - ii) / (cd + c_ak_ext)))
                            vds.append(vds[i] + tau * ((id[i] - ich[i]) / (coss + c_gd_ext)))
                            vgs.append(vgs[i] + tau * (A + B * (vds[i] + vd[i]) + C * (id[i] - ich[i]) + D * vgs[i]))
                            E2A2 = E2A2 + (vds[i] * ich[i] * tau)
                            E2A2new = E2A2 / 1000
                            i = i + 1

                        vgs2A2 = vgs[i]
                        vds2A2 = vds[i]
                        vd2A2 = vd[i]
                        Id2A2 = id[i]
                        t2A2 = (i - 1) * tau
                        TII = t2A1 + t2A2
                        EII = E2A1new + E2A2new
                        dv_dt = (vds2A2 - Vds1B2) / TII

                        print('E2A2 = ', E2A2new)
                        print('Ich2A2 = ', ich[i - 1])
                        print('TII = ', TII)
                        print('EII = ', EII)
                        print('dv_dt = ', dv_dt)

                        #     #----------------------Mode 3------------------------
                         #   -----------------------submode A-------------------

                        Ich3A = ich[i - 1]
                        if Flagich == 0:
                            Ich3A = (0.5 * k_p * (vgs2A2 - vth) ** 2 / (1 + thetaa * (vgs2A2 - vth)))
                        else:
                            Ich3A = 0
                        cosseq3 = (k4 / (1 + (v_dc / k5)) ** m2) + c_gd_ext + k6/((1 + (v_dc/k7)) ** m3)
                        w3 = 1 / (((l_d + l_s) * cosseq3) ** 0.5)
                        r1 = (vds2A2 - v_dc) ** 2
                        ra = (Id2A2 - Ich3A) ** 2
                        r2 = ((l_d + l_s) / cosseq3) * ra
                        A4 = (r1 + r2) ** .5

                        # A4 = (((vds2A2 - v_dc) ** 2) + ((l_d + l_s) / cosseq3) * ((id2A2 - Ich3A) ** 2)) ** .5
                        phi = math.atan(((vds2A2 - v_dc) / (Id2A2 - Ich3A)) * ((cosseq3 / (l_d + l_s)) ** .5))
                        t3A = (1.57 - phi) / w3
                        # E3A = VdcIch3AT3A + (Ld + Ls)(Id2B − Ich3A)Ich3A
                        E3A = (v_dc * Ich3A * t3A) + (l_d + l_s) * (Id2A2 - Ich3A) * Ich3A
                        vds3A = v_dc + A4
                        E3Anew = E3A / 1000

                        TII = t2A1 + t2A2
                        EII = E2A1new + E2A2new
                        #
                        #     #   -----------------------submode B-------------------
                        print("----submode B----")
                        if Flagich == 0:
                            Crssq = k1b / ((1 + ((v_dc - vtd) / k2b)) ** m1b)
                            Cossq = (k4 / ((1 + (v_dc / k5)) ** m2)) + c_gd_ext
                            time_1 = r_g * (c_iss + Crssq) + r_gext * c_gd_ext
                            time_2 = r_g * Crssq + r_gext * c_gd_ext
                            T3B = -(time_1 + k_p * l_s * (v_ee - vth)) * (
                                math.log(1 - (vth - vgs2A2) / (v_ee - vgs2A2))) - k_p * l_s * (vth - vgs2A2)

                            d1 = time_1 + k_p * l_s * (v_ee - vth)
                            d2 = k_p * l_s * (v_ee - vth)
                            d3 = (vgs2A2 - vth) / (v_ee - vth)
                            E3B = (0.5 * k_p * vds3A * ((v_ee - vth) ** 2) * (
                                    d1 * (d3 + (d3 ** 2) * .5 + (math.log(1 - d3))) + d2 * (d3 ** 3) / 3))

                        else:
                            T3B = 0
                            E3B = 0
                        TIII = t3A + T3B
                        EIII = (E3A + E3B) / 1000

                        print("T3B = ", T3B)
                        print("E3B = ", E3B)
                        print("EIII =", EIII)
                        print("TIII =", TIII)
                        print("---final output---")

                        Toff = (TI + TII + TIII)
                        Eoff = (EI + EII + EIII)
                        dv_dt = (vds2A2 - Vds1B2) / TII
                        di_dt = Id2A2 / TIII
                        t_off.append(Toff)
                        e_off.append(Eoff)
                        di_by_dt.append(di_dt)
                        dv_by_dt.append(dv_dt)
                        vds_ov.append(vds3A)


                        print("Toff", Toff)
                        print("dv/dt", dv_dt)
                        print("Eoff", Eoff)
                        print("di/dt", di_dt)
                        print("EII", EII)
                        print("EIII", EIII)

                    if event == 'Submit':
                        otpt = pd.DataFrame(list(zip(*[i_0, t_off, e_off, dv_by_dt, di_by_dt, vds_ov]))).add_prefix('Col')
                        otpt.to_csv('output_orig_a.csv', mode='w', index=False,
                                    header=['I_0(A)', 'T_off(ns)', 'E_off(uJ)', 'dv/dt(V/ns)', 'di/dt(A/ns)', 'Vds_ov(V)'])

                    plt.figure(2)
                    plt.subplot(2, 3, 1)
                    plt.plot(i_0, t_off, marker='o', ms=3, color='black')
                    plt.plot(i_0, t_off, label="T_off vs I_0")
                    plt.xlabel('I_0(A)')
                    plt.ylabel('T_off(ns)')
                    plt.legend(loc='best', fancybox=True, shadow=True)
                    plt.grid(True)

                    plt.subplot(2, 3, 2)
                    plt.plot(i_0, e_off, marker='o', ms=3, color='black')
                    plt.plot(i_0, e_off, label="Eoff vs I_0")
                    plt.xlabel('I_0(A)')
                    plt.ylabel('Eoff(uJ)')
                    plt.legend(loc='best', fancybox=True, shadow=True)
                    plt.grid(True)

                    plt.subplot(2, 3, 3)
                    plt.plot(i_0, dv_by_dt, marker='o', ms=3, color='black')
                    plt.plot(i_0, dv_by_dt, label="dv/dt vs I_0")
                    plt.xlabel('I_0(A)')
                    plt.ylabel('dv/dt(V/ns)')
                    plt.legend(loc='best', fancybox=True, shadow=True)
                    plt.grid(True)

                    plt.subplot(2, 3, 4)
                    plt.plot(i_0, di_by_dt, marker='o', ms=3, color='black')
                    plt.plot(i_0, di_by_dt, label="di/dt vs I_0")
                    plt.xlabel('I_0(A)')
                    plt.ylabel('di/dt(A/ns)')
                    plt.legend(loc='best', fancybox=True, shadow=True)
                    plt.grid(True)

                    plt.subplot(2, 3, 5)
                    plt.plot(i_0, vds_ov, marker='o', ms=3, color='black')
                    plt.plot(i_0, vds_ov, label="Vds_ov vs I_0")
                    plt.xlabel('I_0(A)')
                    plt.ylabel('Vds_ov(V)')
                    plt.legend(loc='best', fancybox=True, shadow=True)
                    plt.grid(True)

                    plt.show()

                    if event == "Cancel" or event == sg.WIN_CLOSED:
                        break

window_first.close()


