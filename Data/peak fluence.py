import numpy as np

P = 5 * 1e-3
f = 500
dx = 440 * 1e-6
dy = 360 * 1e-6
A = np.pi * dx * dy/4
Fpeak = np.log(2) * P / (A * f)
Fpeak_tg = Fpeak * 4

print(f'F single: {Fpeak}\nF tg: {Fpeak_tg}')