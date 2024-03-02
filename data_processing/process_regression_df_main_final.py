#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 09:56:08 2021

@author: tannervarrelman
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import timedelta, date, datetime
import itertools
import matplotlib.dates as mdates
import statsmodels.api as sm
import os
import itertools
from functools import reduce
# import pylogit as pl


def vacc_eff_cli2(data, results, region, iso_3, cli_list, window):
    sympt_dict = {
        'b1_1': 'Fever',
        'b1_2': 'Cough',
        'b1_3': 'Difficulty breathing',
        'b1_4': 'Fatigue',
        'b1_5': 'Stuffy or runny nose',
        'b1_6': 'Aches or muscle pain',
        'b1_7': 'Sore throat',
        'b1_8': 'Chest pain',
        'b1_9': 'Nausea',
        'b1_10': 'Loss of smell or taste',
        'b1_12': 'Headache',
        'b1_13': 'Chills'
        }
    cli_1 = cli_list[0]
    cli_2 = cli_list[1]
    cli_label1 = sympt_dict.get(cli_1)
    cli_label2 = sympt_dict.get(cli_2)
    cli_label = '{0} & {1}'.format(cli_label1, cli_label2)
        # update window script
    for i in range(0, len(window)):
        wave = window.wave_p[i]
        start_date = window.start_date[i]
        end_date = window.end_date[i]
        df1 = data_agg(data, start_date, end_date)
        #have symptoms
        df1.loc[(df1[cli_1]==1) & (df1[cli_2]==1) & (df1['b2']<=14), 'synd_id'] = 1
        df1.loc[(df1[cli_1]!=1) | (df1[cli_2]!=1) | (df1['b2']>14), 'synd_id'] = 0
        df1.loc[(df1[cli_1]==1) & (df1[cli_2]==1) & (df1['b2']<=14), 'Symptom Present'] = 'Yes'
        df1.loc[(df1[cli_1]!=1) | (df1[cli_2]!=1) | (df1['b2']>14), 'Symptom Present'] = 'No'
        # 1 for 1dose 2 for 2dose
        df1.loc[df1['v2']==2, 'vacc_id'] = 1
        df1.loc[df1['v1']==2, 'vacc_id'] = 0
        df1.loc[df1['v2']==2, 'Vaccination Status'] = '2-Dose Vaccinated'
        df1.loc[df1['v1']==2, 'Vaccination Status'] = 'Unvaccinated'
        df1 = df1.dropna(subset = ['vacc_id'])
        # delta wave
        indi_list = ['Omicron Indicator']
        for indi in indi_list:
            if indi == 'Omicron Indicator':
                if wave == 'Delta':
                    df1['wave'] = [0]*len(df1)
                    # omicron wave
                if wave == 'Omicron':
                    df1['wave'] = [1]*len(df1)  
                df1['indicator'] = [indi]*len(df1)
            df1['label'] = [cli_label]*len(df1)
            df1['end_date'] = [end_date]*len(df1)
            df1_sub = df1[['synd_id', 'vacc_id', 'wave', 'iso_3', 'label', 'end_date', 'indicator', 'e3', 'e4', 'Vaccination Status', 'Symptom Present']]
            results = results.append(df1_sub)           
    return results

def process_data(file_path, file):
    iso_code = file[0:3]
    main_df = pd.read_csv(file_path + file)
    main_df['recordeddate'] = pd.to_datetime(main_df['recordeddate']).dt.date
    main_df['b3'] = main_df['b3'].replace([-99, -77, np.nan], 2) 
    main_df['b2'] = main_df['b2'].replace([-99, -77, np.nan], 9999)
    main_df.loc[main_df.b2 < 0, 'b2'] = 9999
    main_df['b1_1'] = main_df['b1_1'].replace([-99, -77, np.nan], 2)
    main_df['b1_2'] = main_df['b1_2'].replace([-99, -77, np.nan], 2)
    main_df['b1_3'] = main_df['b1_3'].replace([-99, -77, np.nan], 2)
    main_df['b1_4'] = main_df['b1_4'].replace([-99, -77, np.nan], 2)
    main_df['b1_5'] = main_df['b1_5'].replace([-99, -77, np.nan], 2)
    main_df['b1_6'] = main_df['b1_6'].replace([-99, -77, np.nan], 2)
    main_df['b1_7'] = main_df['b1_7'].replace([-99, -77, np.nan], 2)
    main_df['b1_8'] = main_df['b1_8'].replace([-99, -77, np.nan], 2)
    main_df['b1_9'] = main_df['b1_9'].replace([-99, -77, np.nan], 2)
    main_df['b1_10'] = main_df['b1_10'].replace([-99, -77, np.nan], 2)
    main_df['b1_12'] = main_df['b1_12'].replace([-99, -77, np.nan], 2)
    main_df['b1_13'] = main_df['b1_13'].replace([-99, -77, np.nan], 2)
    main_df['e3'] = main_df['e3'].replace([-99, -77], np.nan)
    main_df['e4'] = main_df['e4'].replace([-99, -77], np.nan)
    # main_df = main_df[main_df['recordeddate']>="2021-05-01"]
    # Add this step to ensure that the file name is for the correct country
    main_df = main_df[(main_df['iso_3']==iso_code)].reset_index(drop=True)
    if len(main_df) == 0:
        return print('error')
    else:
        return main_df, iso_code

# aggregate data based on the time windows for each wave
def data_agg(df, start_date, end_date):
    end_date = pd.Timestamp(end_date)
    start_date = pd.Timestamp(start_date)
    #print(start_date, end_date)
    data_sub = df[(df['recordeddate']>=start_date) & (df['recordeddate']<=end_date)]
    return data_sub

# get the list of UMD-CTIS line-list files for GTM, MEX, ZAF
file_path = "/Users/tannervarrelman/Documents/Comms_med_VE/data/Countries_2_22_22/"
dir_list = os.listdir(file_path)
file_list = [x for x in dir_list if '_22.csv' in x]

output_path = "/Users/tannervarrelman/Documents/Comms_med_VE/data/output/"
window_df = pd.read_csv(output_path + "GTM_MEX_ZAF_date_windows_median_1_15_23.csv")

sympt_list = ['b1_1', 'b1_2', 'b1_3', 'b1_4', 'b1_5', 'b1_6', 'b1_7',
              'b1_8', 'b1_9', 'b1_10', 'b1_12', 'b1_13']

# create a list of all combinations of symptoms
sympt_combo = list(itertools.combinations(sympt_list, 2))

results_df = pd.DataFrame()
for file in file_list:
    main_df, region = process_data(file_path, file)
    window_sub = window_df[window_df['Region']==region].reset_index(drop=True)
    for cli in sympt_combo:
        results_df = vacc_eff_cli2(main_df, results_df, region, region, cli, window_sub)

results_df.to_csv(output_path + "GTM_MEX_ZAF_regression_df_2dose_8_27_23.csv", header=True, index=False)





