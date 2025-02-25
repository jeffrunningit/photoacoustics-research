import numpy as np

PowerW = 5 * 1e-3
RepetitionFreq = 500
spotSizedx = 440 * 1e-6
spotSizedy = 360 * 1e-6
spotArea = np.pi * spotSizedx * spotSizedy/4
FluencePeak = np.log(2) * PowerW / (spotArea * RepetitionFreq)
FluencePeak_transientGrating = FluencePeak * 4


print("Peak fluence calculation:")
print(f'Power: {PowerW} W\nLaser repetition rate: {RepetitionFreq} Hz\nSpot size\ndx: {spotSizedx:.2}\ndy: {spotSizedy:.2}\nSpot area: {spotArea:.2e}')
print("-----------")
print(f'Peak Fluence single: {FluencePeak:.2f} J/m^2\nPeak Fluence transient grating: {FluencePeak_transientGrating:.2f} J/m^2')