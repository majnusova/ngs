#!/usr/bin/env python
# coding: utf-8

# In[124]:


# https://moshi4.github.io/pyGenomeViz/getting_started/

from pygenomeviz import GenomeViz

gv = GenomeViz()
track = gv.add_feature_track(name="bla", size=30000)  arrow_shaft_ratio=1
fig = gv.plotfig()  # or gv.savefig("test.png") or gv.savefig_html("test.html")


# In[145]:


# vizualizace jednoho genomu
from pygenomeviz import GenomeViz

gv = GenomeViz(tick_style="axis")
track = gv.add_feature_track(name="Chara_n866", size=22855)

track.add_feature(start=1, end=59, strand=1, label="TIR", facecolor="black")
track.add_feature(start=4576, end=5043, strand=1, plotstyle = "bigarrow", label="cys prot", facecolor="orange", linewidth=1)
track.add_feature(start=5614, end=6411, strand=-1, plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1)
track.add_feature(start=7029, end=8642, strand=1,plotstyle = "bigarrow", label="RVE integrase", facecolor="yellow", linewidth=1)
track.add_feature(start=9730, end=10806, strand=1, plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1)
track.add_feature(start=11909, end=12424, strand=1, plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1)
track.add_feature(start=14153, end=16789, strand=1, plotstyle = "bigarrow", label="primase-helicase", facecolor="lightblue", linewidth=1, labelvpos="top")
track.add_feature(start=17024, end=18265, strand=1,plotstyle = "bigarrow", label="MCP1", facecolor="darkgreen", linewidth=1)
track.add_feature(start=18294, end=19541, strand=1, plotstyle = "bigarrow", label="MCP2", facecolor="darkgreen", linewidth=1)
track.add_feature(start=19568, end=20275, strand=1,plotstyle = "bigarrow", label="mCP", facecolor="lightgreen", linewidth=1)
track.add_feature(start=3893, end=4564, strand=-1,plotstyle = "bigarrow", label="pATPase", facecolor="red", linewidth=1, labelvpos="top")
track.add_feature(start=3418, end=3873, strand=-1,plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1)
track.add_feature(start=1962, end=3398, strand=-1,plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1)
track.add_feature(start=22796, end=22855, strand=1, label="TIR", facecolor="black")

fig = gv.plotfig()


# In[147]:


# vizualizace jednoho genomu
from pygenomeviz import GenomeViz

gv = GenomeViz(tick_style="axis")
track = gv.add_feature_track(name="Chara_n866", size=22855)

track.add_feature(start=1, end=59, strand=1, label="TIR", facecolor="black")
track.add_feature(start=4576, end=5043, strand=1, plotstyle = "bigarrow", label="cys prot", facecolor="orange", linewidth=1,  arrow_shaft_ratio=1)
track.add_feature(start=5614, end=6411, strand=-1, plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=7029, end=8642, strand=1,plotstyle = "bigarrow", label="RVE integrase", facecolor="yellow", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=9730, end=10806, strand=1, plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=11909, end=12424, strand=1, plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=14153, end=16789, strand=1, plotstyle = "bigarrow", label="primase-helicase", facecolor="lightblue", linewidth=1, arrow_shaft_ratio=1, labelvpos="top")
track.add_feature(start=17024, end=18265, strand=1,plotstyle = "bigarrow", label="MCP1", facecolor="darkgreen", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=18294, end=19541, strand=1, plotstyle = "bigarrow", label="MCP2", facecolor="darkgreen", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=19568, end=20275, strand=1,plotstyle = "bigarrow", label="mCP", facecolor="lightgreen", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=3893, end=4564, strand=-1,plotstyle = "bigarrow", label="pATPase", facecolor="red", linewidth=1, arrow_shaft_ratio=1, labelvpos="top")
track.add_feature(start=3418, end=3873, strand=-1,plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=1962, end=3398, strand=-1,plotstyle = "bigarrow", label="", facecolor="lightgrey", linewidth=1, arrow_shaft_ratio=1)
track.add_feature(start=22796, end=22855, strand=1, label="TIR", facecolor="black")

fig = gv.plotfig()


# In[125]:


# vizualizace vice genomu najednou s jednim meritkem 
from pygenomeviz import GenomeViz

genome_list = (
    {
        "name": "genome 01", 
        "size": 1000, # total length of genome to be visualized in nt
        "cds_list": [  # each ORF parameters: start, end, strand + parameters (color, shape, label)
            {"range": (150, 300, 1), "facecolor": "red", "plotstyle" : "arrow", "label" : "XX"},
            {"range": (500, 700, -1), "facecolor": "blue", "plotstyle": "bigarrow", "label": "CDS 2"},
            {"range": (750, 950, 1), "facecolor": "green", "plotstyle": "bigarrow", "label": "CDS 3"}
        ]
    },
    {
        "name": "genome 02", 
        "size": 1300, 
        "cds_list": [
            {"range": (50, 200, 1), "facecolor": "purple", "plotstyle": "arrow", "label": "CDS 1"},
            {"range": (350, 450, 1), "facecolor": "orange", "plotstyle": "box", "label": "CDS 2"},
            {"range": (700, 900, -1), "facecolor": "cyan", "plotstyle": "bigbox", "label": "CDS 3"},
            {"range": (950, 1150, -1), "facecolor": "magenta", "plotstyle": "bigrbox", "label": "CDS 4"}
        ]
    },
    {
        "name": "genome 03", 
        "size": 1200, 
        "cds_list": [
            {"range": (150, 300, 1), "facecolor": "yellow", "plotstyle": "arrow", "label": "CDS 1"},
            {"range": (350, 450, -1), "facecolor": "lightblue", "plotstyle": "box", "label": "CDS 2"},
            {"range": (500, 700, -1), "facecolor": "grey", "plotstyle": "bigbox", "label": "CDS 3"},
            {"range": (701, 900, -1), "facecolor": "black", "plotstyle": "bigrbox", "label": "CDS 4"}
        ]
    }
)

gv = GenomeViz(tick_style="axis")
for genome in genome_list:
    name, size, cds_list = genome["name"], genome["size"], genome["cds_list"]
    track = gv.add_feature_track(name, size)
    for cds in cds_list:
        start, end, strand = cds["range"]
        facecolor = cds["facecolor"]
        plotstyle = cds["plotstyle"]
        label = cds["label"]
        track.add_feature(start, end, strand, plotstyle=plotstyle, label=label, facecolor=facecolor, labelcolor="black", labelsize=10, labelvpos="top", labelrotation=0) #arrow_shaft_ratio=1.0 

fig = gv.plotfig()

