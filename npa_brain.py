#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

df = pd.read_csv("npadata.csv")

races = [
    'American Indian or Alaska Native',
    'Asian',
    'Black or African American',
    'Native Hawaiian or Other Pacific Islander',
    'White',
    'Other',
    'Prefer not to answer',
]

mult_arr = []
for i in range(1,8):
    arr = df["prace___{0}".format(i)].to_numpy()
    arr[pd.isnull(df["prace___{0}".format(i)])] = 0
    mult_arr.append(arr)
mult_arr = np.array(mult_arr)   


np.unique(np.sum(mult_arr, axis=0), return_counts = True)
# Counts of pts reporting 0, 1, 2, 3 races
# (array([0., 1., 2., 3.]), array([7589, 4077,   23,    1]))


pt_race = np.array([None] * mult_arr.shape[1])

for i in range(7):
    for j in np.where(mult_arr[i])[0]:
        if pt_race[j] is not None:
            pt_race[j] = pt_race[j] + ", " + races[i]
        else:
            pt_race[j] = races[i]


pt_race[pd.isnull(pt_race)] = "NA"

print(np.array(np.unique(pt_race, return_counts=True)).T)
# Counts of pts reporting each race or combination of races
"""
[['American Indian or Alaska Native' 16]
 ['American Indian or Alaska Native, Asian' 1]
 ['American Indian or Alaska Native, White' 4]
 ['Asian' 156]
 ['Asian, Native Hawaiian or Other Pacific Islander, White' 1]
 ['Asian, Other' 1]
 ['Asian, White' 6]
 ['Black or African American' 457]
 ['Black or African American, Native Hawaiian or Other Pacific Islander' 1]
 ['Black or African American, White' 6]
 ['NA' 7589]
 ['Native Hawaiian or Other Pacific Islander' 15]
 ['Native Hawaiian or Other Pacific Islander, White' 1]
 ['Other' 171]
 ['Prefer not to answer' 322]
 ['White' 2940]
 ['White, Other' 3]]
 """

df.insert(6, "prace", pt_race)

# save to CSV
df.to_csv("npadata_race.csv")




