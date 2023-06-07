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
## argument of size of the figure; You may have to troubleshoot the more proteins/rows you add
parser.add_argument('--size', type=float, nargs=2, default=[6, 16], help='Plot size, default is 6 16 (should be entered as two integers with a space between them)')
## argument of color map of the heatmap
parser.add_argument('--color', type=str, default='vlag', help='Color map, default is vlag')
## argument for the path to the .csv file to use
parser.add_argument('--csv', type=str, help='CSV file path, Create a CSV file with headers containing Visible, Starred, Protein Name, Accession, etc')
## argument for the path to the .png file to save
parser.add_argument('--save', type=str, help='PNG file name to save')
## argument of the numbers of proteins to include, aka the number of rows (The python index starts at 0 instead of 1, so this is the # of rows minus 1)
parser.add_argument('--rows', type=int, default=93, help='Number of rows aka proteins to iclude minus 1 (python index starts at 0 instead of 1), default 93')
## argument of the fontsize for the annotation inside each heatmap box
parser.add_argument('--fontsize', type=int, default=8, help='Font size for the annotations aka spectral counts, default 8')
## argument of the fontsize for the y-axis proteins 
parser.add_argument('--fontsize_tick', type=int, default=8, help='Font size for the y-axis, default 8') ## version 2 addition 
## argument of the location of your sample spectral count columns (again this is the column minus 1 as python starts at 0 instead of 1)
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
print("Top 5 rows of the new concatenated data frame for the heatmap")
print(dfheatmap.head(5))
dfheatmap = dfheatmap.iloc[0:args.rows]
## Set the GeneName column as the row names for each spectral count 
dfheatmap.set_index('GeneName', inplace=True)
print("Top 5 rows of the new df for heatmap with the rows labeled as GeneNames")
print(dfheatmap.head(5))
## Making any NAs into 0s 
dfheatmap_filled = dfheatmap.fillna(0)

# Plot heatmap
offset = 1e-1
## Offsetting the heatmap by 0.01 to help adjust the coloring in the heatmap because the logNorm of 0 is -inf
dfheatmap_offset = dfheatmap_filled + offset
## figure heatmap size 
fig, ax = plt.subplots(figsize=args.size)
## plotting heatmap with the offset df as the coloring because the logNorm of 0 is -inf, the annotated spectral counts as dfheatmap_filled (the actual # of spectral counts) and no decimal point fmt='.0f'; cmap is the coloring map/palette you chose or the default; linecolor is always white with a width of 0.5 between each heatmap square; performing the log normalization of the spectral counts for proper coloring  
p1 = sns.heatmap(dfheatmap_offset, fmt='.0f', annot=dfheatmap_filled, annot_kws={"size": args.fontsize, "weight": "bold"}, cmap=args.color, linecolor='white', linewidth='0.5', norm=LogNorm())
## adjust the size and boldness of the y-axis labeling aka the protein/Gene names 
ax.set_yticklabels(ax.get_yticklabels(), size=args.fontsize_tick, weight='bold') ## version 2 addition 
## Saving the final heatmap 
plt.savefig(args.save)

plt.show()
