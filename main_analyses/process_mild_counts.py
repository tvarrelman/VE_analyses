#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 09:56:08 2021, last edited: 8/27/23

@author: tannervarrelman
"""

import pandas as pd
import numpy as np
import itertools
import os
import csv

def vacc_eff_cli(data, results, region, iso_3, cli_mild, cli_severe, window):
    """
        This function takes in the pre-processed country data ('data'), a dataframe that results
        will be stored in ('results'), the region of interest ('region'), the 
        iso code for that region ('iso_3'), a tuple with the pairwise combo
        of symptoms ('cli_list'), and the time window data describing the waves ('window').
        The function itself aggregates data in preparation for the regression.
    """
    # dictionary of UMD-CTIS symptoms
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
    # unpack the pairwise symptoms passed to the function
    # get the name for the variables
    # format the pairwise label
    cli_label = 'mild infection'
    # iterate over the time windows (delta and omicron)
    with open('/Users/tannervarrelman/Documents/Comms_med_VE/data/output/{0}_regression_df_2dose_mild_9_24_23.csv'.format(iso_3), 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['synd_id', 'vacc_id', 'label', 'wave', 'iso_3', 'e3', 'e4'])
        for i in range(0, len(window)):
            wave = window.wave_p[i]
            start_date = window.start_date[i]
            end_date = window.end_date[i]
            agg_data = data_agg(data, start_date, end_date)
            df1 = agg_data.copy()
            for i in range(0, len(df1)):
                sub_df = df1.iloc[i]
                conditions = [(sub_df[cli_1] == 1) & (sub_df[cli_2] == 1) & (sub_df['b2'] <= 14) for cli_1, cli_2 in cli_mild]
                exclude_conditions = [(sub_df[cli_1] == 1) & (sub_df[cli_2] == 1) & (sub_df['b2'] <= 14) for cli_1, cli_2 in cli_severe]
                result = any(cond.any() for cond in conditions)
                exclude_result = any(cond.any() for cond in exclude_conditions)
                # Final result, including exclusion
                final_result = result and not exclude_result
                if final_result:
                    synd_id = 1
                else:
                    synd_id = 0
                if sub_df['v2'] == 2:
                    vacc_id = 1
                elif sub_df['v1'] == 2:
                    vacc_id = 0
                else:
                    continue
                if wave == 'Delta':
                    wave_id = 0
                elif wave == 'Omicron':
                    wave_id = 1
                csv_writer.writerow([synd_id, vacc_id, cli_label, wave_id, sub_df['iso_3'], sub_df['e3'], sub_df['e4']])

def process_data(file_path, file):
    """
        This function takes a file path, and name of the UMD-CTIS linelist data
        as input. This function reads in the data, and processes the nodata values
        in the linelist.
    """
    iso_code = file[0:3]
    # read the csv file for a country
    main_df = pd.read_csv(file_path + file)
    # turn date string into datetime
    main_df['recordeddate'] = pd.to_datetime(main_df['recordeddate']).dt.date
    # known CLI variable
    main_df['b3'] = main_df['b3'].replace([-99, -77, np.nan], 2) 
    # for how many days have you had at least one of these symptoms
    main_df['b2'] = main_df['b2'].replace([-99, -77, np.nan], 9999)
    # negative values to the duration question are assigned a large value
    # which result in No, for symptoms present within the 14 day timeperiod
    main_df.loc[main_df.b2 < 0, 'b2'] = 9999
    # non response for symptoms are treated as No
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
    # normalize non-response to age and gender questions
    main_df['e3'] = main_df['e3'].replace([-99, -77], np.nan)
    main_df['e4'] = main_df['e4'].replace([-99, -77], np.nan)
    # Add this step to ensure that the file name is for the correct country
    main_df = main_df[(main_df['iso_3']==iso_code)].reset_index(drop=True)
    if len(main_df) == 0:
        return print('error, country code does not match file name')
    else:
        return main_df, iso_code

def data_agg(df, start_date, end_date):
    """
        This function subsets the line-list UMD-CTIS data based on the 
        omicron and delta wave time periods
    """
    end_date = pd.Timestamp(end_date)
    start_date = pd.Timestamp(start_date)
    #print(start_date, end_date)
    data_sub = df[(df['recordeddate']>=start_date) & (df['recordeddate']<=end_date)]
    return data_sub

def main():
    """
        This function runs the data processing program.
    """
    # get the list of UMD-CTIS line-list files for GTM, MEX, ZAF
    file_path = "/Users/tannervarrelman/Documents/Comms_med_VE/data/Countries_2_22_22/"
    dir_list = os.listdir(file_path)
    file_list = [x for x in dir_list if '_22.csv' in x]
    print(file_list)

    output_path = "/Users/tannervarrelman/Documents/Comms_med_VE/data/output/"
    window_df = pd.read_csv(output_path + "GTM_MEX_ZAF_date_windows_median_1_15_23.csv")

    sympt_list = ['b1_1', 'b1_2', 'b1_3', 'b1_4', 'b1_5', 'b1_6', 'b1_7',
                  'b1_8', 'b1_9', 'b1_10', 'b1_12', 'b1_13']

    # create a list of all combinations of symptoms
    sympt_combo = list(itertools.combinations(sympt_list, 2))
    
    sev_list = ['b1_3', 'b1_8']
    mild_list = ['b1_1', 'b1_2', 'b1_6', 'b1_7', 'b1_9', 'b1_10','b1_12']

    # combinations must contain one COVID-19 mild symptom, and no severe symptoms
    mild_combo= [(x,y) for x,y in sympt_combo if (x in mild_list or y in mild_list) and (x not in sev_list and y not in sev_list)]
    # combinations must contain one COVID-19 severe symptom
    severe_combo = [x for x in sympt_combo if 'b1_3' in x or 'b1_8' in x]
    print('mild combinations: ', len(mild_combo))
    print('severe combinations: ', len(severe_combo))
    print('total combinations: ', len(sympt_combo))
    results_df = pd.DataFrame()
    # there is a UMD-CTIS line-list dataset for each country, so iterate over the country data
    for file in file_list:
        # perform pre-processing to handle missing data values
        main_df, region = process_data(file_path, file)
        print('Processing: {0}....'.format(region))
        # subset the window dataframe to find omicron and delta periods for the specified region
        window_sub = window_df[window_df['Region']==region].reset_index(drop=True)
        # iterate over the pairwise combination of symptoms
        vacc_eff_cli(main_df, results_df, region, region, mild_combo, severe_combo, window_sub)

    results_df.to_csv(output_path + "GTM_MEX_ZAF_regression_df_2dose_mild_8_27_23.csv", header=True, index=False)
    print('Successfully processed the UMD-CTIS line-list')
    
if __name__ == '__main__':
    main()

mex_df = pd.read_csv('/Users/tannervarrelman/Documents/Comms_med_VE/data/output/MEX_regression_df_2dose_mild_9_24_23.csv')
gtm_df = pd.read_csv('/Users/tannervarrelman/Documents/Comms_med_VE/data/output/GTM_regression_df_2dose_mild_9_24_23.csv')
zaf_df = pd.read_csv('/Users/tannervarrelman/Documents/Comms_med_VE/data/output/ZAF_regression_df_2dose_mild_9_24_23.csv')

final_df = pd.concat([mex_df, zaf_df, gtm_df])
#final_df.to_csv('/Users/tannervarrelman/Documents/Comms_med_VE/data/output/GTM_MEX_ZAF_regression_df_2dose_mild_9_24_23.csv')

final_df = final_df.dropna()
delta_sub = final_df[final_df['wave']==0]

delta_pos = len(delta_sub[delta_sub['synd_id']==1])
delta_neg = len(delta_sub[delta_sub['synd_id']==0])

omi_sub = final_df[final_df['wave']==1]

omi_pos = len(omi_sub[omi_sub['synd_id']==1])
omi_neg = len(omi_sub[omi_sub['synd_id']==0])





