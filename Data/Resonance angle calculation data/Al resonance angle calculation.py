import numpy as np
import matplotlib.pyplot as plt

wl = 750e-9
n = 1.958355454
k = 7.571403838

pitch = 1620e-9 # grating pitch
N = 1
ed = 1
em = (n**2 - k**2)
theta = np.nan
sameDirection = True

print("Surface plasmon resonance angle calculation:")
print(f"The resonance angle is: {theta:.2f} rad, {theta/(2*np.pi)*360:.2f} deg.")

if em < -1:
    sin1 = np.sqrt(em * ed / (em + ed)) + wl * N / pitch
    sin2 = np.sqrt(em * ed / (em + ed)) - wl * N / pitch
    if np.abs(sin1) < 1:
        theta = np.arcsin(sin1)
        print('Direction is the same.')
        #print(sin1)
    elif np.abs(sin2) < 1:
        theta = np.arcsin(sin2)
        sameDirection = False
        print('Direction is the opposite.')
        #print(sin2)
        
