# Author: Lilianne Nakazono
# This file retrieves the full path to the data and results folders

import os

#%%
mycwd = os.path.abspath(os.getcwd())
print(mycwd.split(os.sep)[-1])

if mycwd.split(os.sep)[-1] == "hostless": 
    parent_cwd = mycwd
    codes_path = os.path.join(parent_cwd,"codes")

if mycwd.split(os.sep)[-1] == "codes":
    # os.chdir('..')
    codes_path = os.path.abspath(os.getcwd())
    os.chdir('..')
    parent_cwd = os.path.abspath(os.getcwd())
    os.chdir(mycwd)

if mycwd.split(os.sep)[-1] == "notebooks":
    os.chdir('..')
    codes_path = os.path.abspath(os.getcwd())
    os.chdir('..')
    parent_cwd = os.path.abspath(os.getcwd())
    os.chdir(mycwd)
#%%

data_path = parent_cwd + '/_data'
if not os.path.exists(data_path):
    os.makedirs(data_path)

result_path = parent_cwd + '/_results'
if not os.path.exists(result_path):
    os.makedirs(result_path)