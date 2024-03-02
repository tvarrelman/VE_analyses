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


def vacc_eff_cli(data, results, region, iso_3, cli_list, window, severity):
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
    cli_1 = cli_list[0]
    cli_2 = cli_list[1]
    # get the name for the variables
    cli_label1 = sympt_dict.get(cli_1)
    cli_label2 = sympt_dict.get(cli_2)
    # format the pairwise label
    cli_label = '{0} & {1}'.format(cli_label1, cli_label2)
    # iterate over the time windows (delta and omicron)
    for i in range(0, len(window)):
        wave = window.wave_p[i]
        start_date = window.start_date[i]
        end_date = window.end_date[i]
        agg_data = data_agg(data, start_date, end_date)
        df1 = agg_data.copy()
        # filter the symptom data (for mild illness, ensure that b1_3 and b1_8 aren't present)
        if severity == 'mild':
            df1.loc[((df1[cli_1]==1) & (df1[cli_2]==1) & (df1['b2']<=14) & (df1['b1_3']!=1) & (df1['b1_8']!=1)), 'synd_id'] = 1
            df1.loc[((df1[cli_1]!=1) | (df1[cli_2]!=1) | (df1['b2']>14) | (df1['b1_3']==1) | (df1['b1_8']==1)), 'synd_id'] = 0
        else:
            df1.loc[((df1[cli_1]==1) & (df1[cli_2]==1) & (df1['b2']<=14)), 'synd_id'] = 1
            df1.loc[((df1[cli_1]!=1) | (df1[cli_2]!=1) | (df1['b2']>14)), 'synd_id'] = 0
        df1.loc[((df1[cli_1]==1) & (df1[cli_2]==1) & (df1['b2']<=14)), 'Symptom Present'] = 'Yes'
        df1.loc[((df1[cli_1]!=1) | (df1[cli_2]!=1) | (df1['b2']>14)), 'Symptom Present'] = 'No'
        # 1 for 1 dose 2 for 2dose
        # v2 == 2 is 2-dose vaccinated in the survey 
        df1.loc[df1['v2']==2, 'vacc_id'] = 1
        # v1 == 2 is unvaccinated in the survey
        df1.loc[df1['v1']==2, 'vacc_id'] = 0
        # create a new column w/ the description of vaccination status
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
            # subset the dataframe to only columns that will be used for analysis
            df1_sub = df1[['synd_id', 'vacc_id', 'wave', 'iso_3', 'label', 'end_date', 'indicator', 'e3', 'e4', 'Vaccination Status', 'Symptom Present']]
            results = results.append(df1_sub)           
    return results

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
    data_sub = df[(df['recordeddate']>=start_date) & (df['recordeddate']<=end_date)]
    return data_sub

def main():
    """
        This function runs the data processing program.
    """
    # get the list of UMD-CTIS line-list files for GTM, MEX, ZAF
    file_path = "../data/"
    dir_list = os.listdir(file_path)
    file_list = [x for x in dir_list if '_22.csv' in x]

    output_path = "../data/output/"
    window_df = pd.read_csv(output_path + "GTM_MEX_ZAF_date_windows_median.csv")

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
    mild_results_df = pd.DataFrame()
    severe_results_df = pd.DataFrame()
    # there is a UMD-CTIS line-list dataset for each country, so iterate over the country data
    for file in file_list:
        # perform pre-processing to handle missing data values
        main_df, region = process_data(file_path, file)
        print('Processing: {0}....'.format(region))
        # subset the window dataframe to find omicron and delta periods for the specified region
        window_sub = window_df[window_df['Region']==region].reset_index(drop=True)
        # iterate over the pairwise combination of symptoms
        for cli in mild_combo:
            mild_results_df = vacc_eff_cli(main_df, mild_results_df, region, region, cli, window_sub, 'mild')
        for cli in severe_combo:
            severe_results_df = vacc_eff_cli(main_df, severe_results_df, region, region, cli, window_sub, 'severe')

    mild_results_df.to_csv(output_path + "GTM_MEX_ZAF_regression_df_2dose_mild.csv", header=True, index=False)
    severe_results_df.to_csv(output_path + "GTM_MEX_ZAF_regression_df_2dose_severe.csv", header=True, index=False)
    print('Successfully processed the UMD-CTIS line-list')
    
if __name__ == '__main__':
    main()





