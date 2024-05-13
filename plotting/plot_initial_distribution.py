"""Module to plot the particle positions from a txt file."""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import seaborn as sns
sns.set_style("dark")


# Get path from argument
args = sys.argv
PATH = '.'

if len(args) > 1:
    PATH = args[1]

# Load data from CSV file
df = pd.read_csv(f'{PATH}/particles.csv')

initial_cong = df[df['time']==1]

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlim((0,10))
ax.set_ylim((0,10))

scat = ax.scatter([], [], c=[], s=1, cmap='PRGn')
scat.set_offsets(np.c_[initial_cong['pos_x'], initial_cong['pos_y']])
scat.set_array(initial_cong['vorticity'])  # Set colors based on vorticity

# Show animation
# plt.show()
plt.savefig(f'{PATH}/init_config.png')