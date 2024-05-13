"""Module to plot the particle positions from a txt file."""
import matplotlib.pyplot as plt
import pandas as pd
import sys

# Get path from argument
args = sys.argv
PATH = '.'

if len(args) > 1:
    PATH = args[1]
    print(PATH)

# Load data from CSV file
df = pd.read_csv(f'{PATH}/energy.csv')

plt.plot(df['energy'])
# Show animation
# plt.show()
plt.savefig(f'{PATH}/energy.png')
