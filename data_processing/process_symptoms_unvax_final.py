#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 15:31:05 2022

@author: tannervarrelman
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os
from itertools import groupby
from operator import itemgetter
import itertools


def cli_calc2(data_sub, df, cli_list, country, region, window_len):
    cli_dict = {
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
        'b1_13': 'Chills',
        'b3': 'CCLI'
        }
    cli_1 = cli_list[0]
    cli_2 = cli_list[1]
    cli_label1 = cli_dict.get(cli_1)
    cli_label2 = cli_dict.get(cli_2)
    cli_label = '{0} & {1}'.format(cli_label1, cli_label2)
    grouped1 = data_sub.groupby([cli_1, cli_2, 'v1', pd.Grouper(key='recordeddate', freq='1D')]).size().to_frame(name = 'count').reset_index().sort_values('recordeddate')
    grouped = grouped1[grouped1['v1']==2]
    tot_df = pd.DataFrame()
    for rec_date in grouped['recordeddate'].unique():
        init_sub = grouped[grouped['recordeddate']==rec_date]
        tot_yes = init_sub[(init_sub[cli_1]==1) & (init_sub[cli_2]==1)]
        tot_no = init_sub[(init_sub[cli_1]!=1) | (init_sub[cli_2]!=1)]
        if len(tot_yes) + len(tot_no) != len(init_sub):
            print('ERROR IN DATA')
        tot_df = tot_df.append({'recordeddate': rec_date, 'total_yes':tot_yes['count'].sum(),
                                'total_no':tot_no['count'].sum()}, ignore_index=True)
    tot_df = tot_df.reset_index().sort_values('recordeddate')
    fev_anos_yes = tot_df['total_yes'].rolling(window=window_len).sum().to_frame(name = 'rolling_yes')
    fev_anos_no = tot_df['total_no'].rolling(window=window_len).sum().to_frame(name = 'rolling_no')
    result = pd.merge(fev_anos_yes, tot_df, left_index=True, right_index=True, how='right').merge(fev_anos_no,left_index=True, right_index=True, how='left')
    for date in result['recordeddate'].unique():
        subset = result[result['recordeddate']==date]
        yes_sub = subset['rolling_yes'].item()
        no_sub = subset['rolling_no'].item()
        if (yes_sub+no_sub) == 0:
            new_val = 0
        else:
            new_val = yes_sub/(yes_sub + no_sub)
        df = df.append({'PropPos':new_val, 'date': date, 'iso_3':country,
                        'variable':(cli_1, cli_2), 'Region':region, 'label':cli_label,
                        'status': 'Unvaccinated'},
                       ignore_index=True)

    return df

def process_data(file_path, file):
    iso_code = file[0:3]
    main_df = pd.read_csv(file_path + file)
    main_df['recordeddate'] = pd.to_datetime(main_df['recordeddate'].tolist())
    main_df['b3'] = main_df['b3'].replace([-99, -77, np.nan], 2) 
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
    # Add this step to ensure that the file name is for the correct country
    main_df = main_df[main_df['iso_3']==iso_code].reset_index(drop=True)
    if len(main_df) == 0:
        return print('error')
    else:
        return main_df, iso_code
    


# path to the source files
file_path = "/Users/tannervarrelman/Documents/Comms_med_VE/data/Countries_2_22_22/"
# each file contains data for particular country
dir_list = os.listdir(file_path)
file_list = [x for x in dir_list if '22.csv' in x]
# list of UMD-CTIS symptoms
sympt_list = ['b1_1', 'b1_2', 'b1_3', 'b1_4', 'b1_5', 'b1_6', 'b1_7',
              'b1_8', 'b1_9', 'b1_10', 'b1_12', 'b1_13']
# pairwise combination of symptoms
sympt_combo = list(itertools.combinations(sympt_list, 2))
cli_df = pd.DataFrame()
for file in file_list:
    iso_df, country = process_data(file_path, file)
    for cli in sympt_combo:
        cli_df = cli_calc2(iso_df, cli_df, cli, country, country, 14)
output_path = "/Users/tannervarrelman/Documents/Comms_med_VE/data/output/"

cli_df.to_csv(output_path + 'GTM_MEX_ZAF_every_cli_unvax_2week_1_15_23.csv', header=True, index=False)


