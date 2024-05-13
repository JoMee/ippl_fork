"""Module to plot the particle positions from a txt file."""
import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns
sns.set_style("whitegrid")

# Get path from argument
args = sys.argv
PATH = '.'

if len(args) > 1:
    PATH = args[1]

# Load data from CSV file
df = pd.read_csv(f'{PATH}/energy.csv')

plt.plot(df['energy'])
plt.xlabel('Time step')
plt.ylabel('Energy')
plt.xlim(0, len(df['energy']))

# Show animation
# plt.show()
plt.tight_layout()
plt.savefig(f'{PATH}/energy.png')
