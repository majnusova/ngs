import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


data_path = '/home/majnusova/Desktop/distr.xlsx' 
df = pd.read_excel(data_path)


df_cleaned = df[['Species*', 'Cilium', 'RAQ', 'RAZ', 'GOR3P']].rename(columns={'Species*': 'Organism'})
#df_cleaned = df_cleaned.fillna('no')  # Nan = means no
df_cleaned.replace({'yes': 1, 'no': 0}, inplace=True)

# yes/no -> 1/0
df_numeric = df_cleaned.set_index('Organism').replace({'yes': 1, 'no': 0})

# plotting the map
plt.figure(figsize=(3, 60))  # map size

# cmap = color, annot = labels (1/0)
ax = sns.heatmap(df_numeric, annot=False, cmap='Pastel2_r', cbar=False, linewidths=.5, linecolor='grey') #cpam binary
ax.set_xticklabels(['Cilium', 'RAQ', 'RAZ', 'GOR3P'], fontsize=15, rotation=45)  # rotating labels
ax.xaxis.tick_top()  # moving column names to the top
ax.set_yticklabels(ax.get_yticklabels(), fontstyle='italic') # row labels -> italics

current_labels = [label.get_text() for label in ax.get_yticklabels()]
ax.set_yticklabels(current_labels, fontstyle='italic', fontsize=15)


plt.savefig('/home/majnusova/Desktop/my_heatmap.pdf', bbox_inches='tight')
plt.show()

