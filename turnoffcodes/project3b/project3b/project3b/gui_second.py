import PySimpleGUI as sg
import pandas as pd
import numpy as np
def second_window():
    layout_second = [[sg.Text('Vgs',font=(30)),
               sg.InputText(font=(30),do_not_clear=False,key="-INI-",justification='left', size=(5, 1)),
               sg.Text('(V)', font=(30)), sg.Text('',size=(0, 0)),
               sg.Text('Dataset for output Characteristics (Id(A)-Vds(V)) (.csv):', size=(50, 1), font=(30)),
               sg.Input(do_not_clear=False,key="-IN-"), sg.FileBrowse(font=(30)),sg.Button('Insert',font=(30))],
              [sg.Text('Please select a row inorder to insert a record  or delete a record .For a particular row,please enter values for both Vgs and the files.', font=(30))],
              [sg.Radio('1',key= 'RADIO1', enable_events = True,size=(1,1),group_id = 'R1'),sg.Text('Vgs:',font=(30)), sg.Text( enable_events=True,key='-OUTPUT1-',size=(20, 1)),
               sg.Text('The File is:',font=(30)), sg.Text( enable_events=True, size=(70, 1), key='-OUTPUT-'),
               sg.Button('Remove',font=(30))],[sg.Radio('2',key= 'RADIO2', enable_events = True,size=(1,1),group_id = 'R1'),sg.Text('Vgs:',font=(30)), sg.Text( enable_events=True,key='-2OUTPUT1-',size=(20, 1)),
               sg.Text('The File is:',font=(30)), sg.Text( enable_events=True, size=(70, 1), key='-2OUTPUT-'),
               sg.Button('Remove',font=(30))],
              [sg.Radio('3',key= 'RADIO3', enable_events = True,size=(1,1),group_id = 'R1'),sg.Text('Vgs:', font=(30)), sg.Text(enable_events=True, key='-3OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-3OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Radio('4', key='RADIO4', enable_events=True, size=(1, 1), group_id='R1'), sg.Text('Vgs:', font=(30)),
               sg.Text(enable_events=True, key='-4OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-4OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Radio('5', key='RADIO5', enable_events=True, size=(1, 1), group_id='R1'), sg.Text('Vgs:', font=(30)),
               sg.Text(enable_events=True, key='-5OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-5OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Radio('6', key='RADIO6', enable_events=True, size=(1, 1), group_id='R1'), sg.Text('Vgs:', font=(30)),
               sg.Text(enable_events=True, key='-6OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-6OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Radio('7', key='RADIO7', enable_events=True, size=(1, 1), group_id='R1'), sg.Text('Vgs:', font=(30)),
               sg.Text(enable_events=True, key='-7OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-7OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Radio('8', key='RADIO8', enable_events=True, size=(1, 1), group_id='R1'), sg.Text('Vgs:', font=(30)),
               sg.Text(enable_events=True, key='-8OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-8OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Radio('9', key='RADIO9', enable_events=True, size=(1, 1), group_id='R1'), sg.Text('Vgs:', font=(30)),
               sg.Text(enable_events=True, key='-9OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-9OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Radio('10', key='RADIO10', enable_events=True, size=(1, 1), group_id='R1'), sg.Text('Vgs:', font=(30)),
               sg.Text(enable_events=True, key='-10OUTPUT1-', size=(20, 1)),
               sg.Text('The File is:', font=(30)), sg.Text(enable_events=True, size=(70, 1), key='-10OUTPUT-'),
               sg.Button('Remove', font=(30))],
              [sg.Text('', font=(30), size=(20, 1)), sg.Submit(font=(30)), sg.Text('', font=(30), size=(10, 1))]]



    window = sg.Window('Second Form', layout_second)
    while True :
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == 'Cancel':
            break
        if (event == 'Submit'):
            break
        # if (event['SUBMIT'] == 'Submit'):
        #     break

        if (window['RADIO1'].get())== True and event == 'Insert':
            window['-OUTPUT-'].update(values["-IN-"])
            window['-OUTPUT1-'].update(values["-INI-"])
            window['RADIO1'].update(False)

        if (window['RADIO2'].get()) == True and event == 'Insert':
            window['-2OUTPUT-'].update(values["-IN-"])
            window['-2OUTPUT1-'].update(values["-INI-"])
            window['RADIO2'].update(False)

        if (window['RADIO3'].get()) == True and event == 'Insert':
            window['-3OUTPUT-'].update(values["-IN-"])
            window['-3OUTPUT1-'].update(values["-INI-"])
            window['RADIO3'].update(False)

        if (window['RADIO4'].get()) == True and event == 'Insert':
            window['-4OUTPUT-'].update(values["-IN-"])
            window['-4OUTPUT1-'].update(values["-INI-"])
            window['RADIO4'].update(False)

        if (window['RADIO5'].get()) == True and event == 'Insert':
            window['-5OUTPUT-'].update(values["-IN-"])
            window['-5OUTPUT1-'].update(values["-INI-"])
            window['RADIO5'].update(False)

        if (window['RADIO6'].get()) == True and event == 'Insert':
            window['-6OUTPUT-'].update(values["-IN-"])
            window['-6OUTPUT1-'].update(values["-INI-"])
            window['RADIO6'].update(False)

        if (window['RADIO7'].get()) == True and event == 'Insert':
            window['-7OUTPUT-'].update(values["-IN-"])
            window['-7OUTPUT1-'].update(values["-INI-"])
            window['RADIO7'].update(False)

        if (window['RADIO8'].get()) == True and event == 'Insert':
            window['-8OUTPUT-'].update(values["-IN-"])
            window['-8OUTPUT1-'].update(values["-INI-"])
            window['RADIO8'].update(False)

        if (window['RADIO9'].get()) == True and event == 'Insert':
            window['-9OUTPUT-'].update(values["-IN-"])
            window['-9OUTPUT1-'].update(values["-INI-"])
            window['RADIO9'].update(False)

        if (window['RADIO10'].get()) == True and event == 'Insert':
            window['-10OUTPUT-'].update(values["-IN-"])
            window['-10OUTPUT1-'].update(values["-INI-"])
            window['RADIO10'].update(False)



        if event == 'Remove' and window['RADIO1'].get() == True:
            window['-OUTPUT1-'].update('')

        if event == 'Remove' and window['RADIO1'].get() == True:
            window['-OUTPUT-'].update('')

        if event == 'Remove0' and window['RADIO2'].get() == True:
            window['-2OUTPUT1-'].update('')

        if event == 'Remove0' and window['RADIO2'].get() == True:
            window['-2OUTPUT-'].update('')


        if event == 'Remove1' and window['RADIO3'].get() == True:
            window['-3OUTPUT1-'].update('')

        if event == 'Remove1' and window['RADIO3'].get() == True:
            window['-3OUTPUT-'].update('')

        if event == 'Remove2' and window['RADIO4'].get() == True:
            window['-4OUTPUT1-'].update('')

        if event == 'Remove2' and window['RADIO4'].get() == True:
            window['-4OUTPUT-'].update('')

        if event == 'Remove3' and window['RADIO5'].get() == True:
            window['-5OUTPUT1-'].update('')

        if event == 'Remove3' and window['RADIO5'].get() == True:
            window['-5OUTPUT-'].update('')

        if event == 'Remove4' and window['RADIO6'].get() == True:
            window['-6OUTPUT1-'].update('')

        if event == 'Remove4' and window['RADIO6'].get() == True:
            window['-6OUTPUT-'].update('')

        if event == 'Remove5' and window['RADIO7'].get() == True:
            window['-7OUTPUT1-'].update('')

        if event == 'Remove5' and window['RADIO7'].get() == True:
            window['-7OUTPUT-'].update('')

        if event == 'Remove6' and window['RADIO8'].get() == True:
            window['-8OUTPUT1-'].update('')

        if event == 'Remove6' and window['RADIO8'].get() == True:
            window['-8OUTPUT-'].update('')

        if event == 'Remove7' and window['RADIO9'].get() == True:
            window['-9OUTPUT1-'].update('')

        if event == 'Remove7' and window['RADIO9'].get() == True:
            window['-9OUTPUT-'].update('')

        if event == 'Remove8' and window['RADIO10'].get() == True:
            window['-10OUTPUT1-'].update('')

        if event == 'Remove8' and window['RADIO10'].get() == True:
            window['-10OUTPUT-'].update('')

    v_gs11 = []
    id1_vds11 = []
    v_gs11.append(window['-OUTPUT1-'].get())
    v_gs11.append(window['-2OUTPUT1-'].get())
    v_gs11.append(window['-3OUTPUT1-'].get())
    v_gs11.append(window['-4OUTPUT1-'].get())
    v_gs11.append(window['-5OUTPUT1-'].get())
    v_gs11.append(window['-6OUTPUT1-'].get())
    v_gs11.append(window['-7OUTPUT1-'].get())
    v_gs11.append(window['-8OUTPUT1-'].get())
    v_gs11.append(window['-9OUTPUT1-'].get())
    v_gs11.append(window['-10OUTPUT1-'].get())

    id1_vds11.append(window['-OUTPUT-'].get())
    id1_vds11.append(window['-2OUTPUT-'].get())
    id1_vds11.append(window['-3OUTPUT-'].get())
    id1_vds11.append(window['-4OUTPUT-'].get())
    id1_vds11.append(window['-5OUTPUT-'].get())
    id1_vds11.append(window['-6OUTPUT-'].get())
    id1_vds11.append(window['-7OUTPUT-'].get())
    id1_vds11.append(window['-8OUTPUT-'].get())
    id1_vds11.append(window['-9OUTPUT-'].get())
    id1_vds11.append(window['-10OUTPUT-'].get())

    window.close()
    return v_gs11, id1_vds11



