import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Data input (averaged for simplicity)
data = {
    'Integration': [
        '5\'_3wk', '5\'_3wk', '5\'_3wk',
        '5\'_8wk', '5\'_8wk', '5\'_8wk',
        '3\'_3wk', '3\'_3wk', '3\'_3wk',
        '3\'_8wk', '3\'_8wk', '3\'_8wk'
    ],
    'Type': [
        'Forward_F9', 'Reverse_F9', 'ITR',
        'Forward_F9', 'Reverse_F9', 'ITR',
        'Forward_F9', 'Reverse_F9', 'ITR',
        'Forward_F9', 'Reverse_F9', 'ITR'
    ],
    'Average_Percentage': [
        # 5' Integration, 3wk averages
        (1.1994003 + 1.03425986 + 0.86289549 + 0.59016393) / 4,  # Forward_F9
        (1.34932534 + 0.06464124 + 1.05465005 + 0.72131148) / 4,  # Reverse_F9
        (0.67466267 + 0 + 0.86289549 + 0.91803279) / 4,  # ITR
        # 5' Integration, 8wk averages
        (1.37875101 + 1.6064257 + 1.17647059 + 1.1011011) / 4,  # Forward_F9
        (1.45985401 + 0.46854083 + 0.48442907 + 0.1001001) / 4,  # Reverse_F9
        (0.81103001 + 0.46854083 + 0 + 0) / 4,  # ITR
        # 3' Integration, 3wk averages
        (1.1994003 + 1.098901099 + 3.451581975 + 1.967213115) / 4,  # Forward_F9
        (0.299850075 + 0 + 1.534036433 + 1.901639344) / 4,  # Reverse_F9
        (0.224887556 + 0 + 0.287631831 + 0.786885246) / 4,  # ITR
        # 3' Integration, 8wk averages
        (3.487429035 + 2.409638554 + 2.422145329 + 9.90990991) / 4,  # Forward_F9
        (1.540957015 + 1.004016064 + 0.207612457 + 1.301301301) / 4,  # Reverse_F9
        (0.648824006 + 0.401606426 + 0.207612457 + 0.600600601) / 4  # ITR
    ]
}

# Create DataFrame
df = pd.DataFrame(data)

# Create pivot table for the heatmap using named parameters
pivot_df = df.pivot(index='Type', columns='Integration', values='Average_Percentage')

# Plotting the heatmap
plt.figure(figsize=(10, 6), dpi=600)
ax = sns.heatmap(pivot_df, annot=True, cmap="OrRd", fmt=".3f", cbar_kws={'label': 'Average Percentage of Integration'})
plt.title('Heatmap of Gene Integration by Type and Time Point')
plt.savefig('/Users/gwisna/Desktop/heatmap_DNA_int.png', format='png', dpi=600)  # Save as PNG with high DPI
plt.savefig('/Users/gwisna/Desktop/heatmap_DNA_int.svg', format='svg')
plt.show()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype

# Detailed data points matched to categories
data = {
    'Integration': (
        ['5\'_3wk']*4 + ['5\'_3wk']*4 + ['5\'_3wk']*4 +
        ['5\'_8wk']*4 + ['5\'_8wk']*4 + ['5\'_8wk']*4 +
        ['3\'_3wk']*4 + ['3\'_3wk']*4 + ['3\'_3wk']*4 +
        ['3\'_8wk']*4 + ['3\'_8wk']*4 + ['3\'_8wk']*4
    ),
    'Type': (
        ['Forward_F9']*4 + ['Reverse_F9']*4 + ['ITR']*4 +
        ['Forward_F9']*4 + ['Reverse_F9']*4 + ['ITR']*4 +
        ['Forward_F9']*4 + ['Reverse_F9']*4 + ['ITR']*4 +
        ['Forward_F9']*4 + ['Reverse_F9']*4 + ['ITR']*4
    ),
    'Data': [
        1.1994003, 1.03425986, 0.86289549, 0.59016393,  # 5'_3wk Forward_F9
        1.34932534, 0.06464124, 1.05465005, 0.72131148, # 5'_3wk Reverse_F9
        0.67466267, 0, 0.86289549, 0.91803279,          # 5'_3wk ITR
        1.37875101, 1.6064257, 1.17647059, 1.1011011,   # 5'_8wk Forward_F9
        1.45985401, 0.46854083, 0.48442907, 0.1001001,  # 5'_8wk Reverse_F9
        0.81103001, 0.46854083, 0, 0,                   # 5'_8wk ITR
        1.1994003, 1.098901099, 3.451581975, 1.967213115, # 3'_3wk Forward_F9
        0.299850075, 0, 1.534036433, 1.901639344,       # 3'_3wk Reverse_F9
        0.224887556, 0, 0.287631831, 0.786885246,       # 3'_3wk ITR
        3.487429035, 2.409638554, 2.422145329, 9.90990991, # 3'_8wk Forward_F9
        1.540957015, 1.004016064, 0.207612457, 1.301301301, # 3'_8wk Reverse_F9
        0.648824006, 0.401606426, 0.207612457, 0.600600601  # 3'_8wk ITR
    ]
}

# Create DataFrame
df = pd.DataFrame(data)

# Define the desired order of categories
type_order = CategoricalDtype(['Forward_F9', 'Reverse_F9', 'ITR'], ordered=True)
df['Type'] = df['Type'].astype(type_order)

# Calculate the median values for each group
median_df = df.groupby(['Integration', 'Type']).median().reset_index()

# Create pivot table for the heatmap using the median values
pivot_df = median_df.pivot(index='Type', columns='Integration', values='Data')

# Set font to Arial
plt.rcParams["font.family"] = "Arial"

# Plotting the heatmap
plt.figure(figsize=(10, 6), dpi=600)
ax = sns.heatmap(pivot_df, annot=True, cmap="OrRd", fmt=".3f", cbar_kws={'label': 'Median Percentage of Integration'})
plt.title('Heatmap of Gene Integration by Type and Time Point')
plt.savefig('/Users/gwisna/Desktop/heatmap_DNA_int_median.pdf', format='pdf', dpi=600)  # Save as PDF with high DPI
plt.show()