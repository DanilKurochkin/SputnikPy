import numpy as np
from classes import light_solver as ls

Tmin_avg = 265
Tmax_avg = 281
Tmin = 250
Tmax = 300

labels = ["Зенит", "Надир", "Перпенд 1", "Перпенд 2", "Теневая", "Солнечная"]
args = ["min", "avg", "max"]

new_columns_for_T = []
for arg in args:
    for label in labels:
        new_columns_for_T.append(f"{label}_{arg}")

solver = ls.Solver(
    width=0.005, Lx=1, Ly=1, Lz=1,
    c=800, p=2700, R=2,
    orbit_height=500,
    generation_heat=50*np.ones(6))

df, As_cols, e_cols = solver.construct_df(
    [0.1, 0.366, 0.633, 0.9],
    [0.1, 0.366, 0.633, 0.9],
    Tmin_avg, Tmax_avg)
print(f"Finded: {len(df)} combination")

df[new_columns_for_T] = np.vectorize(
    solver.temperature_diapozone,
    signature='(m), (m)->(n)')(df[As_cols], df[e_cols])

df['Tmin'] = df[new_columns_for_T].min(axis=1)
df['Tmax'] = df[new_columns_for_T].max(axis=1)
df.to_csv('result/pre_result.csv')

df = df[(df['Tmin'] > Tmin) & (df['Tmax'] > Tmax)]
print(f"Matching with requiremnts count: {len(df)}")
df.to_csv('result/result.csv')