import datetime as dt
import pyresito.io.ecl as ecl
import csv
import numpy as np
import sys


# setup_pipt.py: Script to generate input for the pipt toolbox
# Command line input to the program: data_file "data_dict", e.g.
# > python3 setup_pipt.py EG_MODEL_IOR "{'WOPRH':100**2, 'WWPRH':200**2, 'WGPRH':20000**2}"
# The case must have been run for the script to work. The output is:
# data_index.csv, true_data.csv, true_var.csv, data_types.csv, actnum.npz
# The actnum is based on the information in the .INIT file

def main(argv):

    # get the dates
    if len(argv) < 2:
        print('No input file specified')
        sys.exit()
    datafile = argv[1]
    case = ecl.EclipseCase(datafile)
    start = case.start_date()
    dates = case.report_dates()

    # get the assimilation time
    assim_time = []
    assim_time_str = []
    for i in range(len(dates) - 1):
        delta = dates[i + 1] - dates[0]
        assim_time.append(delta.days)
        assim_time_str.append(str(delta.days))

    # write the assimilation time
    with open('data_index.csv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(assim_time_str)

    # get the wells and producers
    wells = case.field_data('ZWEL', dates[-1])
    wells = [w for w in wells if w]
    ih = case.field_data('INTEHEAD', dates[-1])
    iw = case.field_data('IWEL', dates[-1])
    NWELLS = ih[16]
    NIWELZ = ih[24]
    producers = []
    for i in range(NWELLS):
        if iw[6 + i * NIWELZ] == 1:
            producers.append(wells[i])

    # get the data types that should be collected, and the variance for the data
    if len(argv) < 3:
        data_dict = {'WOPRH': 100**2, 'WWPRH': 200**2, 'WGPRH': 20000**2}
    else:
        data_dict = eval(argv[2])
    data_list = data_dict.keys()

    # write files
    f = open('true_data.csv', 'w', newline='')
    g = open('true_var.csv', 'w', newline='')
    h = open('data_types.csv', 'w', newline='')
    writer1 = csv.writer(f)
    writer2 = csv.writer(g)
    writer3 = csv.writer(h)
    list_datatyp = []
    for ind, time in enumerate(assim_time):
        tmp_data = []
        tmp_var = []
        for data in data_list:
            for well in producers:
                if _check_open(datafile, case.by_date[start + dt.timedelta(days=time)], well):
                    # noinspection PyBroadException
                    try:
                        single_data = case.summary_data(data + ' ' + well, start + dt.timedelta(days=time))
                    except:
                        print('Data type ' + data + ' for well ' + well, ' not available!')
                        sys.exit()
                    if len(single_data):
                        tmp_data.extend(single_data)
                        tmp_var.extend(['ABS', data_dict[data]])
                    else:
                        tmp_data.extend(['N/A'])
                        tmp_var.extend(['ABS', 0])
                else:
                    tmp_data.extend(['N/A'])
                    tmp_var.extend(['ABS', 0])
                if ind == 0:
                    list_datatyp.extend([data + ' ' + well])

        writer1.writerow(tmp_data)
        writer2.writerow(tmp_var)

    writer3.writerow(list_datatyp)

    f.close()
    g.close()
    h.close()

    # write the actnum file
    f = ecl.EclipseInit(datafile)
    v = f.actnum
    actnum = v.flatten()
    np.savez('actnum', actnum=actnum)


# sub function: check if well is open
def _check_open(datafile, time_int, well):
    f = ecl.EclipseFile(datafile, "X{0:04d}".format(time_int))
    wnames = list(filter(None, f.get('ZWEL')))  # list of defined wells
    niwelz = f.get('INTEHEAD')[24]  # number of data elements per well in IWEL array
    # nwells = f.get('INTEHEAD')[16] #  number of wells
    status = 0  # initialize as of
    for n, w in enumerate(wnames):
        if w == well:
            status = f.get('IWEL')[n * niwelz + 10]  # get well status
            break
    if status > 0:
        return True
    else:
        return False


main(sys.argv)
