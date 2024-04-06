import numpy as np

from classes import light_solver as ls
from classes import utils as utils

solver = ls.Solver(
    width=0.005,
    Lx=1,
    Ly=1,
    Lz=1,
    c=800,
    p=2700,
    R=2,
    orbit_height=2,
    generation_heat=50*np.ones(6))

df, As_cols, e_cols = utils.construct_As_e_df([0.1, 0.9], [0.1, 0.9], 6)

print(df)

for i in range(len(df)):
    if i % 100 == 0:
        print(f"Iteration : {i} of {len(df)}")
    #средняя
    df.loc[i, utils.new_columns_for_T] = solver.temperature_diapozone(df[As_cols].iloc[i].values, df[e_cols].iloc[i].values)

df.to_csv('result.csv')