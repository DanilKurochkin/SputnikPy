from itertools import product
from copy import deepcopy
import pandas as pd
import numpy as np

labels = ["Зенит", "Надир", "Перпенд 1", "Перпенд 2", "Теневая", "Солнечная"]
args = ["min", "avg", "max"]

new_columns_for_T = []
for arg in args:
    for label in labels:
        new_columns_for_T.append(f"{label}_{arg}")


def recursive_comb(l, n):
    res_arr = []
    
    def rec(l, curr, n):
        if n == 0:
            res_arr.append(curr)
            return

        for e in l:
            curr_iter = deepcopy(curr)
            curr_iter.append(e)
            rec(l, curr_iter, n-1)
    
    rec(l, [], n)
    
    return res_arr

def all_combinations(x, y):
    comb = []
    
    for x_val in x:
        for y_val in y:
            comb.append([x_val, y_val])

    
    return comb

def extract(matrix):
    x_mat = []
    y_mat = []
    
    for array in matrix:
        x_array = []
        y_array = []
        for i in range(len(array)):
            x_array.append(array[i][0])
            y_array.append(array[i][1])
        x_mat.append(x_array)
        y_mat.append(y_array)
    
    return x_mat, y_mat

def construct_As_e_df(As, e, n):
    combinations = recursive_comb(all_combinations(As, e), 6)

    As_array, e_array = extract(combinations)

    df = pd.concat([pd.DataFrame(As_array), pd.DataFrame(e_array)], axis=1)
    new_columns_name_As = [f"As_{i}" for i in range(6)]
    new_columns_name_e = [f"e_{i}" for i in range(6)]
    new_columns_name = new_columns_name_As + new_columns_name_e

    df.columns = new_columns_name

    df[new_columns_for_T] = float(0)
    
    return df, new_columns_name_As, new_columns_name_e

