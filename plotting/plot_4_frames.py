"""Module to plot the particle positions from a txt file."""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import seaborn as sns
sns.set_style("dark")

# Set the font size
plt.rcParams.update({'font.size': 12})

# Get path from argument
args = sys.argv
PATH = '.'

if len(args) > 1:
    PATH = args[1]

# Load data from CSV file
df = pd.read_csv(f'{PATH}/particles.csv')

# Create a figure with colorbar and 4 subplots
fig, ax = plt.subplots(1, 4, figsize=(15, 4))

fig.subplots_adjust()
for i, frame in enumerate([1, 300, 600, 900]):
    initial_cong = df[df['time']==frame]
    scat = ax[i].scatter(
        initial_cong['pos_x'],
        initial_cong['pos_y'],
        c=initial_cong['vorticity'],
        s=1, cmap='viridis',
    )
    # Add grid
    ax[i].grid(True)

    # Delete x and y axis
    ax[i].set_xlim((0,10))
    ax[i].set_ylim((0,10))
    ax[i].set_title(f'Time: {frame}', loc='right')

    # Divide grid in 4
    ax[i].set_xticks(np.arange(0, 10, 10/4))
    ax[i].set_yticks(np.arange(0, 10, 10/4))
    ax[i].set_xticklabels([])
    ax[i].set_yticklabels([])


# Add color bar with title
# cbar_ax = fig.add_axes([0.93, 0.1, 0.01, 0.8])
# cbar = fig.colorbar(scat, cax=cbar_ax)
# cbar.set_label('Vorticity', labelpad=15)

plt.tight_layout(h_pad=4)#, rect=[0, 0, 0.92, 1])
plt.savefig(f'{PATH}/4_frames.png')