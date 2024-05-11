"""Module to plot the particle positions from a txt file."""
import matplotlib.pyplot as plt
import pandas as pd

# Change this path to the path of the txt file
PATH="../build_serial/alvine"

# Load data from CSV file
df = pd.read_csv(f'{PATH}/energy.csv')

plt.plot(df['energy'])

plt.savefig('energy.png')