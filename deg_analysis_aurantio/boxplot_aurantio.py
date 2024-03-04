#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt #pro prejmenovani os

sns.set_style("whitegrid")

norm_counts = pd.read_csv('/home/majnusova/all/projects/jotnarlogs/scripts/aurantio_norm_counts_.csv', index_col=0)

specimen = pd.Series({"motile1": "motile", "motile2": "motile", "motile3": "motile", "nonmotile1": "nonmotile","nonmotile2": "nonmotile","nonmotile3": "nonmotile"}, name="specimen")

genes = norm_counts.transpose()
genes = pd.concat([genes], axis=1)
df = pd.concat([genes, specimen], axis=1)


df

sns.boxplot(y="e_gw1.9.395.1", x="specimen", data=df)
plt.ylabel("RAQ")

sns.boxplot(y="fgenesh1_pg.22_#_86", x="specimen", data=df)
plt.ylabel("RAW")

sns.boxplot(y="estExt_Genemark1.C_2_t10283", x="specimen", data=df)
plt.ylabel("GOR3P")

