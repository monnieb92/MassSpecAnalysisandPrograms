# MassSpecAnalysisandPrograms
Different scripts to use to analyze or visualize mass spec data sets
## massspec_heatmap.py
### Start with the CSV file
Use the following screenshot below from excel as to how your table should look before saving as a .csv file:
![Screen Shot 2023-06-07 at 11 26 06 AM](https://github.com/monnieb92/MassSpecAnalysisandPrograms/assets/60197768/9fabe202-bcee-4149-82c3-b54f0fa9cc61)

### Determine on the location of the Spectral count columns 
In this case they are 9:15, which is the default 

### Determine on the how many proteins you want to show in the heatmap 
In this case they are 0:93 (total of 94), which is the default 

### Select what color to make the heatmap or leave with default 

### Write a name for the output file .png 

### Defining the code: 

import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from matplotlib.colors import LogNorm

#### Parse command-line arguments to allow you to customize the size of the figure, color of the heatmap, name, number of proteins, and number of samples in your heatmap 
parser = argparse.ArgumentParser(description='Generate and save a heatmap plot')
parser.add_argument('--info', action='help', help='--size nargs=2 default=[6, 16] Plot size --color default=vlag Color map; --csv CSV file path ; --save PNG file name to save; --rows default=93 number of rows minus 1')  # Add --help option
parser.add_argument('--size', type=float, nargs=2, default=[6, 16], help='Plot size')
parser.add_argument('--color', type=str, default='vlag', help='Color map')
parser.add_argument('--csv', type=str, help='CSV file path, Create a CSV file with headers containing Visible, Starred, Protein Name, Accession, etc')
parser.add_argument('--save', type=str, help='PNG file name to save')
parser.add_argument('--rows', type=int, default=93, help='Number of rows minus 1')
parser.add_argument('--fontsize', type=int, default=8, help='Font size for the annotations, default 8')
parser.add_argument('--locat', type=str, default='9:15', help='Location of spectral count columns for the heatmap, default 9:15 (This assumes ParentalA, ParentalB, ParentalC, SampleA, SampleB, SampleC)')
args = parser.parse_args()

#### Read CSV file (printing out the header will allow you to see if you csv file looks correct)
df = pd.read_csv(args.csv)
print("Top 10 rows of the csv file")
print(df.head(10))

#### Extract GeneName (this extracts the gene name from the Identified protein column) ex: Protein CBFA2T3 OS=Homo sapiens OX=9606 GN=CBFA2T3 PE=1 SV=2![image](https://github.com/monnieb92/MassSpecAnalysisandPrograms/assets/60197768/42c81fea-7e5e-466b-abbf-5709491b3d37)

df['GeneName'] = df.iloc[:,3].str.extract(r'GN=(\w+)')

#### Prepare heatmap data, this selects only the the columns of the Samples with the spectral counts, the GeneNames, and the number of proteins to show in the rows. 
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

#### Plot heatmap, this plots the heatmap using the log normalization for the spectral counts to better visualize the color, and showing the raw spectral counts in each box 
offset = 1e-1
dfheatmap_offset = dfheatmap_filled + offset
fig, ax = plt.subplots(figsize=args.size)
p1 = sns.heatmap(dfheatmap_offset, fmt='.0f', annot=dfheatmap_filled, annot_kws={"size": args.fontsize, "weight": "bold"}, cmap=args.color, linecolor='white', linewidth='0.5', norm=LogNorm())

plt.savefig(args.save)

plt.show()
