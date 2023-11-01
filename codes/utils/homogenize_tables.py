# Author: Robert Beck 
# Refactored by: Lilianne Nakazono

import os
import sys
import io
import gzip

from pathlib import Path
import pandas as pd
import numpy as np
import astropy.io.fits as fits

from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

if str(os.path.join(Path.cwd(),"codes")) not in sys.path:
    sys.path.append(str(os.path.join(Path.cwd(),"codes")))
from settings.column_names import TIME_SERIES_COLUMNS, STAMP_COLUMNS, REDUCE_COLUMNS

def process_simbad_tns_pickle(path_obj:str)->pd.DataFrame:
    """_summary_ 

    Args:
        path_obj (str): path to a .pkl file

    Returns:
        pd.DataFrame: _description_
    """

    df = pd.read_pickle(str(path_obj))
    
    id_col = 'i:objectId'
    for col in TIME_SERIES_COLUMNS:
        if col not in df:
            df[col] = None

    # Sort data by julian date (! not tested)
    df.sort_values(by='i:jd', inplace=True)

    # Aggregate data from same object id into a tuple
    ts_df = df[[id_col] + TIME_SERIES_COLUMNS].groupby(id_col, as_index=False).aggregate(tuple)
    
    # Convert tuple to numpy array
    for col in TIME_SERIES_COLUMNS:
        ts_df[col] = ts_df[col].apply(np.array)
    
    # Get first and last stamp data by julian date
    jd_groupby = df[[id_col, 'i:jd']].groupby(id_col)
    first_df = df.loc[jd_groupby.idxmin().values.flatten(), [id_col] + STAMP_COLUMNS]
    last_df = df.loc[jd_groupby.idxmax().values.flatten(), [id_col] + STAMP_COLUMNS]
    
    # Convert first and lat stamp data to numpy array
    for col in STAMP_COLUMNS:
        first_df[col] = first_df[col].apply(np.array)
        last_df[col] = last_df[col].apply(np.array)

    # Rename columns to identify first and last stamps
    first_df = first_df.rename(columns={col: f'{col}_first' for col in STAMP_COLUMNS})
    last_df = last_df.rename(columns={col: f'{col}_last' for col in STAMP_COLUMNS})
    
    # Merge data by object ID
    obj_df = ts_df.merge(first_df, on=id_col).merge(last_df, on=id_col)
    return obj_df

def read_spicy_fits_image(bytes_str:str):
    """_summary_

    Args:
        bytes_str (str): _description_

    Returns:
        _type_: _description_
    """
    hdu_list = fits.open(gzip.open(io.BytesIO(bytes_str)))
    primary_hdu = hdu_list[0]
    return primary_hdu.data

def homogenize_spicy_df(xmatch_spicy:pd.DataFrame)->pd.DataFrame:
    """_summary_

    Args:
        xmatch_spicy (pd.DataFrame): _description_

    Returns:
        pd.DataFrame: _description_
    """
    for col in TIME_SERIES_COLUMNS:
        if col not in xmatch_spicy:
            xmatch_spicy[col] = xmatch_spicy['i:aimage'].apply(lambda x: np.array([None]*len(x)))
    for col in STAMP_COLUMNS:
        orig_col = f'{col}_small'
        xmatch_spicy[f'{col}_first'] = xmatch_spicy[orig_col].apply(lambda x: read_spicy_fits_image(x[0]))
    for col in STAMP_COLUMNS:
        orig_col = f'{col}_small'
        xmatch_spicy[f'{col}_last'] = xmatch_spicy[orig_col].apply(lambda x: read_spicy_fits_image(x[-1]))
        xmatch_spicy = xmatch_spicy.drop(columns=[orig_col])
    xmatch_spicy = xmatch_spicy.rename(columns={'objectId': 'i:objectId'})
    return xmatch_spicy

def filter_na_get_unique(x, raise_on_duplicate=False):
    """_summary_

    Args:
        x (_type_): _description_
        raise_on_duplicate (bool, optional): _description_. Defaults to False.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    x_sr = pd.Series(x)
    uniques = np.unique(x[(x_sr != 'None').values & (x_sr != 'nan').values & (~pd.isna(x))])
    if len(uniques) > 1:
        message = f'Multiple matches {uniques} contained in array {x}'
        if raise_on_duplicate:
            raise ValueError(message)
        else:
            print(message)
            return None
    elif len(uniques) == 1:
        return uniques[0]
    else:
        return None

def get_classification_tuple(x):
    """_summary_

    Args:
        x (_type_): _description_

    Raises:
        KeyError: _description_
        KeyError: _description_

    Returns:
        _type_: _description_
    """
    count_available = 0
    count_incl_unknown = 0
    most_frequent_value = None
    most_frequent_count = 0
    
    for value, count in zip(*np.unique(x, return_counts=True)):
        if (value != 'None') and (value != 'nan') and not pd.isna(value):
            count_incl_unknown += count
            if (value != 'Unknown'):
                count_available += count
                if count > most_frequent_count:
                    most_frequent_value = value
                    most_frequent_count = count
    
    return (most_frequent_value, most_frequent_count/max(1, count_available), most_frequent_count/max(1,count_incl_unknown))

def add_reduced_classification_cols(df):
    """_summary_

    Args:
        df (_type_): _description_
    """
    class_col = 'v:classification'
    out_cols = [
        'v:classification_best',
        'v:classification_prob',
        'v:classification_prob_incl_unknown',
    ]
    df[out_cols[0]], df[out_cols[1]], df[out_cols[2]] = zip(*df[class_col].apply(get_classification_tuple))
    return

def homogenize_tables(simbad_path:str='/media3/CRP7/hosts/data/SIMBAD/Apr2023/obj_info', tns_path:str='/media3/CRP7/hosts/data/TNS/Apr2023/objects.pickle', spicy_path:str='/media3/CRP7/hosts/data/SPICY/SPICY_CROSSMATCHED_2_ASEC_SMALL/', output_path:str='/media3/CRP7/hosts/data/homogenized/', n_workers=8, verbose=True):
    """_summary_

    Args:
        simbad_path (str, optional): _description_. Defaults to '/media3/CRP7/hosts/data/SIMBAD/Apr2023/obj_info'.
        tns_path (str, optional): _description_. Defaults to '/media3/CRP7/hosts/data/TNS/Apr2023/objects.pickle'.
        spicy_path (str, optional): _description_. Defaults to '/media3/CRP7/hosts/data/SPICY/SPICY_CROSSMATCHED_2_ASEC_SMALL/'.
        output_path (str, optional): _description_. Defaults to '/media3/CRP7/hosts/data/homogenized/'.
        n_workers (int, optional): _description_. Defaults to 8.
        verbose (bool, optional): _description_. Defaults to True.

    Raises:
        KeyError: _description_
        KeyError: _description_
    """
    if verbose:
        print(f'Started processing TNS tables at {datetime.now()}')
    tns_path_objs = Path(tns_path).glob('TNS*.pickle')
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        obj_df_iterator = executor.map(process_simbad_tns_pickle, tns_path_objs)
        xmatch_tns = pd.concat(obj_df_iterator, ignore_index=True)
    if verbose:
        print(f'Finished processing TNS table at {datetime.now()}')

        print(f'Started processing SIMBAD tables at {datetime.now()}')
    simbad_path_objs = Path(simbad_path).glob('finkclass=*/*.pickle')
    with ProcessPoolExecutor(max_workers=8) as executor:
        obj_df_iterator = executor.map(process_simbad_tns_pickle, simbad_path_objs)
        xmatch_simbad = pd.concat(obj_df_iterator, ignore_index=True)
    if verbose:
        print(f'Finished processing SIMBAD table at {datetime.now()}')

        print(f'Started processing SPICY table at {datetime.now()}')
    xmatch_spicy = pd.read_parquet(spicy_path)
    xmatch_spicy = homogenize_spicy_df(xmatch_spicy)
    if verbose:
        print(f'Finished processing SPICY table at {datetime.now()}')

    try:
        list(xmatch_tns.columns) == list(xmatch_simbad.columns) == list(xmatch_spicy.columns)
    except:
        raise KeyError('Error in homogenizing tables. Columns of xmatch_tns, xmatch_simbad, and xmatch_spicy are not the same due to a possible error in the script.')

    for col in REDUCE_COLUMNS:
        xmatch_tns[col] = xmatch_tns[col].apply(filter_na_get_unique)
        xmatch_simbad[col] = xmatch_simbad[col].apply(filter_na_get_unique)
        xmatch_spicy[col] = xmatch_spicy[col].apply(filter_na_get_unique)

    add_reduced_classification_cols(xmatch_tns)
    add_reduced_classification_cols(xmatch_simbad)
    add_reduced_classification_cols(xmatch_spicy)

    try:
        list(xmatch_tns.columns) == list(xmatch_simbad.columns) == list(xmatch_spicy.columns)
    except:
        raise KeyError('Error in reducing columns. Columns of xmatch_tns, xmatch_simbad, and xmatch_spicy are not the same due to a possible error in the script.')

    os.makedirs(output_path, exist_ok=True)
    with open(os.path.join(output_path, 'tns.pkl'), 'wb') as out_file:
        xmatch_tns.to_pickle(out_file)
    with open(os.path.join(output_path, 'simbad.pkl'), 'wb') as out_file:
        xmatch_simbad.to_pickle(out_file)
    with open(os.path.join(output_path, 'spicy.pkl'), 'wb') as out_file:
        xmatch_spicy.to_pickle(out_file)