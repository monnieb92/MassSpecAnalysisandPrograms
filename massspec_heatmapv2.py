#!/bin/python3

import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate and save a heatmap plot')
parser.add_argument('--info', action='help', help='--size nargs=2 default=[6, 16] Plot size --color default=vlag Color map; --csv CSV file path ; --save PNG file name to save; --rows default=93 number of rows minus 1')  # Add --help option
parser.add_argument('--size', type=float, nargs=2, default=[6, 16], help='Plot size')
parser.add_argument('--color', type=str, default='vlag', help='Color map')
parser.add_argument('--csv', type=str, help='CSV file path, Create a CSV file with headers containing Visible, Starred, Protein Name, Accession, etc')
parser.add_argument('--save', type=str, help='PNG file name to save')
parser.add_argument('--rows', type=int, default=93, help='Number of rows minus 1')
parser.add_argument('--fontsize', type=int, default=8, help='Font size for the annotations, default 8')
parser.add_argument('--fontsize_tick', type=int, default=8, help='Font size for the yaxis, default 8') ## version 2 addition 
parser.add_argument('--locat', type=str, default='9:15', help='Location of spectral count columns for the heatmap, default 9:15 (This assumes ParentalA, ParentalB, ParentalC, SampleA, SampleB, SampleC)')
args = parser.parse_args()

# Read CSV file
df = pd.read_csv(args.csv)
print("Top 10 rows of the csv file")
print(df.head(10))

# Extract GeneName
df['GeneName'] = df.iloc[:,3].str.extract(r'GN=(\w+)')

# Prepare heatmap data
dfheatmapv1 = df['GeneName'].to_frame(name='GeneName')
loc_start, loc_end = map(int, args.locat.split(':'))
dfheatmap = df.iloc[:, loc_start:loc_end]  # Select columns 9 to 14
dfheatmap = pd.concat([dfheatmapv1,dfheatmap], axis = 1)
print(dfheatmap.head(5))
dfheatmap = dfheatmap.iloc[0:args.rows]
dfheatmap.set_index('GeneName', inplace=True)
print("Top 5 rows of the new df for heatmap")
print(dfheatmap.head(5))
dfheatmap_filled = dfheatmap.fillna(0)

# Plot heatmap
offset = 1e-1
dfheatmap_offset = dfheatmap_filled + offset
fig, ax = plt.subplots(figsize=args.size)
p1 = sns.heatmap(dfheatmap_offset, fmt='.0f', annot=dfheatmap_filled, annot_kws={"size": args.fontsize, "weight": "bold"}, cmap=args.color, linecolor='white', linewidth='0.5', norm=LogNorm())
ax.set_yticklabels(ax.get_yticklabels(), size=args.fontsize_tick, weight='bold') ## version 2 addition 
plt.savefig(args.save)

plt.show()
