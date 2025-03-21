# %%
import numpy as np

import matplotlib.pyplot as plt

# Function to read the rr variable from the files
def read_file(filename):
    rr_values = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 5:
                rr = float(parts[4])
                rr_values.append(rr)
    return np.array(rr_values)

# Reading rr from each file
outer_long = read_file("outer_anode_long.txt")
outer_short = read_file("outer_anode_short.txt")
inner_long = read_file("inner_anode_long.txt")
inner_short = read_file("inner_anode_short.txt")


plt.figure(figsize=(10, 4))
plt.hist(outer_long, bins=150, alpha=0.7, label='Outer Long')
plt.hist(inner_long, bins=150, alpha=0.7, label='Inner long')
plt.xlabel('Energy')
plt.ylabel('count')
plt.title('Outer Long vs Inner Long')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 4))
plt.hist(outer_short, bins=150, alpha=0.7, label='Outer Short')
plt.hist(inner_short, bins=150, alpha=0.7, label='Inner Short')
plt.xlabel('Energy')
plt.ylabel('count')
plt.title('Outer Short vs Inner Short')
plt.legend()
plt.grid(True)
plt.show()
