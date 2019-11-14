__author__ = 'tubh'
import sys, os
sys.path.append('/home/tubh/Ekofisk/pipt/')
#from simulator.subsurf_flow import ecl_100
import datetime as dt
import input_output.ecl as ecl
import csv
import numpy as np

def main():
    start = dt.datetime(1970, 1, 5)
    case = ecl.EclipseCase('../RunEkofisk_Reference/PETW_DEC14_RESV_V4')

    assim_time = []

    with open('../RunEkofisk_Reference/data_index.csv') as f:
        reader = csv.reader(f,delimiter=',')
        for rec in reader:
            for el in rec:
                assim_time.append(int(el))

    with open('../RunEkofisk_Reference/well_names.txt','r') as f:
        prod_wells = [line.strip() for line in f]

    # RJL: This is synthetic!
    prod_data = ['WOPR', 'WWPR', 'WGPR']

    #inj_wells = ['C-4H', 'C-1H', 'C-2H', 'C-3H', 'F-1H', 'F-2H', 'F-3H', 'F-4H', 'C-4AH']
    #inj_data = ['WWIRH', 'WGIRH']


    all_var = {'WOPR': 100**2,
               'WWPR': 200**2,
               'WGPR': 20000**2,
    }

    f = open('true_data.csv', 'w', newline='')
    g = open('true_var.csv', 'w', newline='')
    writer1 = csv.writer(f)
    writer2 = csv.writer(g)
    for ind, time in enumerate(assim_time):
        tmp_data = []
        tmp_var = []
        list_datatyp = []
        for data in prod_data:
            for well in prod_wells:
                if _check_open(case.by_date[start + dt.timedelta(days=time)], well):
                    single_data = case.summary_data(data + ' ' + well, start + dt.timedelta(days=time))
                    if len(single_data):
                        tmp_data.extend(single_data)
                        tmp_var.extend(['ABS', all_var[data]])
                    else:
                        tmp_data.extend(['N/A'])
                        tmp_var.extend(['ABS', 0])
                else:
                    tmp_data.extend(['N/A'])
                    tmp_var.extend(['ABS', 0])
                list_datatyp.extend([data + ' ' + well])

        writer1.writerow(tmp_data)
        writer2.writerow(tmp_var)

 #       if time == 8:
 #           np.savez('data_typ_upd', list_datatyp=list_datatyp)
    f.close()
    g.close()

def _check_open(time_int,well):
    f = ecl.EclipseFile('../RunEkofisk_Reference/PETW_DEC14_RESV_V4', "X{0:04d}".format(time_int))
    wnames = list(filter(None, f.get('ZWEL')))  # list of defined wells
    niwelz = f.get('INTEHEAD')[24]  # number of data elements per well in IWEL array
    # nwells = f.get('INTEHEAD')[16] #  number of wells
    status = 0 # initiallize as of
    for n, w in enumerate(wnames):
        if w == well:
            status = f.get('IWEL')[n*niwelz + 10] # get well status
            break
    if status > 0:
        return True
    else:
        return False

main()
