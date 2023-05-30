
import sys
import numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

pwd = sys.argv[1]
structure_list_file = sys.argv[2]

vh = np.array([0.000,0.325,0.839])  # plddt > 90, hex: 0053d6
h  = np.array([0.396,0.749,0.953])  # 90 > plddt > 70, hex: 65bff3
l  = np.array([1.000,0.859,0.075])  # 70 > plddt > 50, hex: ffdb13
vl = np.array([1.000,0.490,0.271])  # plddt < 50, hex: ff7d45


# ----------------------------------------
# filling the color dictionaries
# ----------------------------------------
all_colorids = list(range(33,1057))
plddt_range = np.linspace(0,100,len(all_colorids))

color_defs = {}
for colorid, plddt in zip(all_colorids,plddt_range):
    index = colorid - all_colorids[0]
    if plddt < 50.:
        color_defs[colorid] = f'color change rgb {colorid} {vl[0]} {vl[1]} {vl[2]}\n'
    elif plddt < 70.:
        color_defs[colorid] = f'color change rgb {colorid} {l[0]} {l[1]} {l[2]}\n'
    elif plddt < 90.:
        color_defs[colorid] = f'color change rgb {colorid} {h[0]} {h[1]} {h[2]}\n'
    else:
        color_defs[colorid] = f'color change rgb {colorid} {vh[0]} {vh[1]} {vh[2]}\n'


# ----------------------------------------
# MAIN
# ----------------------------------------
# gathering models' paths for models to be analyzed/visualized
with open(structure_list_file,'r') as lst_file:
    model_paths = [Path(line.split()[0]) for line in lst_file.readlines()]

for model_path in model_paths:
    # get the output path
    output_path = model_path.parent
    model_name = model_path.stem
    # open and write the vis state file
    with open(f'{str(output_path)}/{model_name}_quality_vis.vmd','w') as vmd_out:
        # writing preamble
        vmd_out.write('#!/usr/local/bin/vmd\nset viewplist {}\nset fixedlist {}\n\n')
        ### prepping the AF model molecule and reps
        af_model_string = ''
        af_model_string += '### AF MODEL ###\n'
        af_model_string += 'mol new ' + str(model_path) + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'
        af_model_string += 'mol rename top ' + str(model_name) + '_AF_model\n'
        af_model_string += 'mol delrep 0 top\n'
        # starting representation
        af_model_string += 'mol representation NewCartoon 0.160000 50.000000 4.100000 0\n'
        af_model_string += 'mol color Beta\n'
        af_model_string += 'mol selection {all}\n'
        af_model_string += 'mol material AOEdgy\n'
        af_model_string += 'mol addrep top\n'
        vmd_out.write(af_model_string)
        # writing color range out to file
        vmd_out.write('# setting colorid rgb values\n')
        for i in all_colorids:
            vmd_out.write(color_defs[i])


# ----------------------------------------
# creating the mpl cmap
# ----------------------------------------
plddt_cdict = {'red': [],'green': [], 'blue': []}
plddt_range = np.linspace(50,100,len(all_colorids))

# vl is the first element in the cdict
plddt_cdict['red'].append((0,vl[0],vl[0]))
plddt_cdict['green'].append((0,vl[1],vl[1]))
plddt_cdict['blue'].append((0,vl[2],vl[2]))

for plddt in plddt_range:
    if plddt < 70.:
        plddt_cdict['red'].append(  ((plddt-50)/50.,l[0],l[0]))
        plddt_cdict['green'].append(((plddt-50)/50.,l[1],l[1]))
        plddt_cdict['blue'].append( ((plddt-50)/50.,l[2],l[2]))
    elif plddt < 90.:
        plddt_cdict['red'].append(  ((plddt-50)/50.,h[0],h[0]))
        plddt_cdict['green'].append(((plddt-50)/50.,h[1],h[1]))
        plddt_cdict['blue'].append( ((plddt-50)/50.,h[2],h[2]))
    else:
        plddt_cdict['red'].append(  ((plddt-50)/50.,vh[0],vh[0]))
        plddt_cdict['green'].append(((plddt-50)/50.,vh[1],vh[1]))
        plddt_cdict['blue'].append( ((plddt-50)/50.,vh[2],vh[2]))

cmap = mpl.colors.LinearSegmentedColormap('my_cmap',plddt_cdict,len(all_colorids))
fig, ax = plt.subplots(figsize=(2,8))
fig.subplots_adjust(right=0.5)
norm = mpl.colors.Normalize(vmin=50,vmax=100)
cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,spacing='uniform',norm=norm,orientation='vertical',ticks=[50,70,90,100],extend='min')
cb.set_label('pLDDT Score',size=16)
cb.set_ticklabels(['50','70','90','100'],size=14)
plt.savefig(f'alphafold_plddt_colorbar.png',dpi=600,transparent=True)
plt.close()


