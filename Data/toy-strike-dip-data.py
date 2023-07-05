
import numpy as np
import csv

# strike azumuths at 10 degree intervals
regular_strikes = list(np.repeat(np.linspace(0,360,num=36,endpoint=False),10))

# dips at 10 degree intervals for each strike azumuth
regular_dips = list(np.linspace(0,90,num=10,endpoint=False)) * 36

merged_list = list(zip(regular_strikes,regular_dips))
print(merged_list)

with open('toy-strike-dip-data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['strike','dip'])
    writer.writerows(merged_list)