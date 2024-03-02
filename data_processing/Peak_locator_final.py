#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 10:03:13 2022

@author: tannervarrelman
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import argrelextrema
from datetime import timedelta, date, datetime

output_path = "../data/output/"
sympt_df = pd.read_csv(output_path + 'GTM_MEX_ZAF_every_cli_unvax_2week.csv')

window_df = pd.DataFrame()
peak_list = []
no_peak = []
for sympt in sympt_df['label'].unique():
    for region in sympt_df['Region'].unique():
        # subset by symptom and region
        st_sub = sympt_df[(sympt_df['label']==sympt) & (sympt_df['Region']==region)]
        st_sub['date'] = pd.to_datetime(st_sub['date'])
        # set date as index. Important for indexing peak values
        st_sub.set_index('date', inplace=True) 
        # we use freq_df to simply plot the CTIS time-series (PropPos is smoothed prop. of respondents)
        freq_df = st_sub[['PropPos']].reset_index()
        # locate the local max (order=70 consistently provides 2 peaks)
        ilocs_max = argrelextrema(freq_df.PropPos.values, np.greater_equal, order=70)[0]
        # find the peak values diven the index
        max_df = freq_df.iloc[ilocs_max].reset_index()
        # check to make sure that we only locate two peaks

        if len(max_df) == 2:
            peak_list.append(1)
        else:
            no_peak.append(1)
        # loop through peak values and append to df 
        for i in range(0, len(max_df)):
            if i == 0:
                wave_p = 'Delta'
            if i == 1:
                wave_p = 'Omicron'
            window = max_df['date'][i]
            value = max_df['PropPos'][i]
            window_df = window_df.append({'windows': window, 'Region': region, 
                                          'symptom':sympt, 'value': value, 'wave_p': wave_p}, 
                                         ignore_index=True)
        # Plot the local max on the time-series to ensure that everything worked as expected
        fig, ax = plt.subplots(figsize=(20,8))
        freq_df.plot(x='date', y='PropPos', alpha=.3, ax=ax)
        max_df.plot(x='date', y='PropPos', style='.', lw=10, color='red', marker="v", ax=ax)
        title = '{0}: {1}'.format(region, sympt)
        fig.suptitle(title, fontsize=16)

window_filter = pd.DataFrame()        
for region in window_df.Region.unique():
    for wave in window_df.wave_p.unique():
        window_sub = window_df[(window_df.Region==region) & (window_df.wave_p==wave)]
        med_window = window_sub.windows.median()
        start_date = pd.to_datetime(med_window) - timedelta(weeks=2)
        window_filter = window_filter.append({'end_date':med_window,
                                              'start_date':start_date,
                                              'Region': region,
                                              'wave_p': wave}, ignore_index=True)
   
window_filter.to_csv(output_path + 'GTM_MEX_ZAF_date_windows_median.csv', index=False, header=True)  


