
import time
import pickle
from pathlib import Path
import sys
import MDAnalysis
import numpy as np
import pandas
import logging
import traceback

import matplotlib as mpl
import matplotlib.pyplot as plt
import PIL
from PIL import Image, ImageDraw, ImageFont, ImageColor
import ast
import math

import usalign_parser

pwd = sys.argv[1]
pdb70_metadata_pkl_file = sys.argv[2]
ranking_file_name = sys.argv[3]
tmscore_cutoff = float(sys.argv[4])


###################
# FUNCTIONS
###################

def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s    %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


# ----------------------------------------
# setting logging file
# ----------------------------------------

# set up the main logger file and list all relevant parameters.
main_logger = setup_logger('vis_state_creator','vis_state_creator.log')
start_time = time.time()
main_logger.info(f'Starting analysis.')
main_logger.info(f'Working directory: {pwd}')
main_logger.info(f'TMscore cutoff: {tmscore_cutoff}')
main_logger.info(f'Making color scales.')


# ----------------------------------------
# setting up color and radius ranges
# ----------------------------------------
color_defs = {}
#high_counts_color = np.array([0.0,0.0,1.0])         # blue
#low_counts_color  = np.array([0.827,0.827,0.827])   # light gray
#high_cons_color = np.array([1.0,0.0,0.0])           # red
#low_cons_color  = np.array([0.827,0.827,0.827])     # light gray
#high_counts_color = np.array([0.302,0.686,0.290])   # green #4daf4a
#low_counts_color  = np.array([0.827,0.827,0.827])   # light gray
#high_cons_color = np.array([0.596,0.306,0.639])     # purple #984ea3
#low_cons_color  = np.array([0.827,0.827,0.827])     # light gray
high_counts_color = np.array([1.0,0.0,0.0])         # red
low_counts_color  = np.array([0.827,0.827,0.827])   # light gray
high_cons_color = np.array([0.0,0.0,1.0])           # blue
low_cons_color  = np.array([0.827,0.827,0.827])     # light gray
min_radius = 0.1
max_radius = 1.5
radius_diff = max_radius - min_radius
x_space = 0.01
bg = (255, 255, 255, 0)     # transparent white


# ----------------------------------------
# filling the color dictionaries
# ----------------------------------------
all_colorids = list(range(33,1057))
# feature counts
feature_colorids = list(range(33,545))
nFeatureColorids = len(feature_colorids)
cmap_positions = np.linspace(0,1,nFeatureColorids)

feature_cdict = {'red':[], 'green':[], 'blue':[]}
for colorid in feature_colorids:
    index = colorid - feature_colorids[0]
    thiscolor = np.abs(low_counts_color + index * (high_counts_color - low_counts_color)/(nFeatureColorids-1))
    ### VMD Stuff
    color_defs[colorid] = f"color change rgb {colorid} {thiscolor[0]} {thiscolor[1]} {thiscolor[2]}\n"
    ### MATPLOTLIB Stuff
    feature_cdict['red'].append((cmap_positions[index],thiscolor[0],thiscolor[0]))
    feature_cdict['green'].append((cmap_positions[index],thiscolor[1],thiscolor[1]))
    feature_cdict['blue'].append((cmap_positions[index],thiscolor[2],thiscolor[2]))

# conservation counts
cons_colorids = list(range(545,1057))
nConsColorids = len(cons_colorids)
cmap_positions = np.linspace(0,1,nConsColorids)

cons_cdict = {'red':[], 'green':[], 'blue':[]}
for colorid in cons_colorids:
    index = colorid - cons_colorids[0]
    thiscolor = np.abs(low_cons_color + index * (high_cons_color - low_cons_color)/(nConsColorids-1))
    ### VMD Stuff
    color_defs[colorid] = f"color change rgb {colorid} {thiscolor[0]} {thiscolor[1]} {thiscolor[2]}\n"
    ### MATPLOTLIB Stuff
    cons_cdict['red'].append((cmap_positions[index],thiscolor[0],thiscolor[0]))
    cons_cdict['green'].append((cmap_positions[index],thiscolor[1],thiscolor[1]))
    cons_cdict['blue'].append((cmap_positions[index],thiscolor[2],thiscolor[2]))


# ----------------------------------------
# beginning the analysis
# ----------------------------------------

main_logger.info(f'Loading the metadata pickle file.')
with open(pdb70_metadata_pkl_file,'rb') as pkl_file:
    pdbid_chainid_metadata = pickle.load(pkl_file)

main_logger.info(f'Finding genes to visualize.')
# ----------------------------------------
# NOTE: this only works because of the set directory naming convention
#genes = [direc.name for direc in Path(pwd).glob('Sphm*/')]
genes = [direc.name for direc in Path(pwd).glob('Sphm02G160700.*/')]
# ----------------------------------------

main_logger.info(f'Looping over {len(genes)} genes.')
for gene in genes:
    # setting up gene-specific logging file
    gene_logger = setup_logger('{gene}',f'{pwd}/{gene}/vis_state_creator.log')
    gene_st = time.time()
    gene_logger.info(f'Parsing alignment results and creating VMD vis state files for {gene}.')
    main_logger.info(f'Parsing alignment results and creating VMD vis state files for {gene}.')
    
    # ----------------------------------------
    # NOTE: this only works because of the set directory naming convention
    # get the model pdb file, the ranking file, and the path to the 
    # subdirectory where alignment hits are stored
    gene_model = list(Path(f'{pwd}/{gene}/').glob('model_*pdb'))[0]
    gene_ranking = Path(f'{pwd}/{gene}/{ranking_file_name}')
    gene_aln_path = Path(f'{pwd}/{gene}/structural_aln_results/')
    # ----------------------------------------
    
    gene_logger.info(f'AF Model file: {gene_model}')
    gene_logger.info(f'AF Ranking file: {gene_ranking}')
    gene_logger.info(f'AF Alignment subdir: {gene_aln_path}')
    
    # ----------------------------------------
    ## NOTE: this only works because of the set directory naming convention
    ## I specifically numbered the pdb files by their rank... this should change
    ## get the alignment hit pdb files, sorted by quality, high to low
    #hits = sorted(list(Path(f'{pwd}/{gene}/structural_aln_results/').glob('*/*pdb')),key=lambda x: x.stem)
    
    # read ranking csv and apply the tmscore_cutoff to the avgTMscore value
    df = pandas.read_csv(str(gene_ranking), usecols = ['Target Path','avgTMscore'])
    hits = df[df['avgTMscore'] > tmscore_cutoff]

    # test to see if df is empty; if it is, clean logger and move on
    if hits.empty:
        gene_logger.info(f"No hits to be visualized. Closing out. Took {time.time() - gene_st} seconds.")
        clean_logger(gene_logger)
        continue

    # ----------------------------------------

    gene_logger.info(f'Grabbing the Center of Geometry (CoG) coordinates of residues in {gene_model}.')
    # get the model's residues' center of mass to use for visualizing feature counts
    af_model_resid_coord = {}
    af_model_features_dict = {}
    u = MDAnalysis.Universe(gene_model)
    u_all = u.select_atoms('all')
    for i in range(u_all.n_residues):
        resid = str(u_all.residues[i].resid)
        af_model_resid_coord[resid] = u_all.residues[i].atoms.center_of_geometry()
        af_model_features_dict[resid] = {'cons_counts':0,
                                         'feat_counts':0,
                                         'feat_string':"",
                                         'feat_types':[]}

    gene_logger.info(f'Opening the vis state file at {pwd}/{gene}/aln_hits.vmd .')
    # open vmd vis state file to be written to
    with open(f'{pwd}/{gene}/aln_hits.vmd','w') as vmd_out:
        # writing preamble
        vmd_out.write('#!/usr/local/bin/vmd\nset viewplist {}\nset fixedlist {}\n\n')
        
        ### prepping the AF model molecule and reps
        af_model_string = ''
        af_model_string += '### AF MODEL ###\n'
        af_model_string += 'mol new ' + str(gene_model) + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'
        af_model_string += 'mol rename top ' + str(gene) + '_AF_model\n'
        af_model_string += 'mol delrep 0 top\n'
        # starting representation
        af_model_string += 'mol representation NewCartoon 0.160000 50.000000 4.100000 0\n'
        af_model_string += 'mol color ColorID 26\n'  # violet2
        af_model_string += 'mol selection {all}\n'
        af_model_string += 'mol material AOEdgy\n'
        af_model_string += 'mol addrep top\n'

        gene_logger.info(f"Beginning to loop over alignment hits to parse their UniProt features. Will denote 'important' features with ***. If the feature's residues are within the alignment landmark mapping, then including those residues in the counts metrics.")
        # loop over the alignment hit structures gather the feature details and write out to vis state file
        for index, hit in hits.iterrows():
        #for hit in hits:
        #    # ----------------------------------------
        #    # NOTE: may have to change this if the filenames change
        #    pdbid_chainid = hit.stem[4:]    # only because the files are named this way
        #    # need to analyze the usalign log file to get the mapping
        #    aln_log = str(hit.parent / 'usalign.log')
        #    # ----------------------------------------
            
            hit_pdbid_chainid = Path(hit['Target Path']).stem
            aln_log = Path(f'{pwd}/{gene}/structural_aln_results/{hit_pdbid_chainid}/usalign.log')
            aln_pdb = Path(f'{pwd}/{gene}/structural_aln_results/{hit_pdbid_chainid}/{hit_pdbid_chainid}.pdb')
            # gather the parsed alignment results
            results = usalign_parser.parse_usalign_file(aln_log,'sNS')
            aln_results = results[0]
            # aln_results is a dictionary w/ format:
            # keys: "struct1", "struct2", "struct1_chainID", "struct2_chainID", "tmscore1",  "tmscore2", "Len1",  "Len2", "d0_1",  "d0_2", "map_1_to_2" dictionary
            # we really care about the "map_1_to_2" dictionary
            # this dict's keys are query structure resids, values are a tuple w/ format:
            # elem[0] --> target structure residue *** important to us
            # elem[1] --> query resname     *** important to us
            # elem[2] --> target resname    *** important to us
            # elem[3] --> *optional* alignment distance
            gene_logger.info('\n------------------------------------------------------------------------------------------------')
            gene_logger.info(f"Parsing {hit_pdbid_chainid}. Alignment scores are: {[aln_results['tmscore1'],aln_results['tmscore2']]}")
            #if np.mean([aln_results['tmscore1'],aln_results['tmscore2']]) < tmscore_cutoff:
            #    gene_logger.info(f"Average alignment scores is below the TMscore cutoff {tmscore_cutoff}, moving on without parsing this hit.\n")
            #    continue

            # will be used to aggregate conserved residues
            hit_mapping_conservation = []

            # prep the basic molecule and starting representation for this hit
            hit_string = ''
            hit_string += 'mol new ' + str(aln_pdb) + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'
            hit_string += 'mol delrep 0 top\n'
            hit_string +=f'mol rename top {index}_{hit_pdbid_chainid}\n'
            hit_string += 'mol representation NewCartoon 0.160000 50.000000 4.100000 0\n'
            hit_string += 'mol color ColorID 19\n'   # green2
            hit_string += 'mol selection {all}\n'
            hit_string += 'mol material AOEdgy\n'
            hit_string += 'mol addrep top\n'
           
            # parse hit's uniprot metadata to get feature information
            metadata = pdbid_chainid_metadata[hit_pdbid_chainid]
            # check to see if a UniProt ID was associated with the pdbid_chainid
            if 'uniprotid_aln' in metadata:
                uniprotid = list(metadata['uniprotid_aln'].keys())[0]
                
                # avoid an issue if the uniprotid == None
                if not uniprotid:
                    continue

                gene_logger.info(f"    {hit_pdbid_chainid} is associated with {uniprotid}")
                # appending the starting representation to a comment line that contains meta-information
                hit_string = f'### {hit_pdbid_chainid} - {uniprotid} ###\n' + hit_string
                    
                # check to see if the features list is filled with elements 
                try:
                    if metadata['uniprotid_aln'][uniprotid]['features']:
                        for feature in metadata['uniprotid_aln'][uniprotid]['features']:
                            # feature list w/ format:
                            # elem[0] --> Type of feature
                            # elem[1] --> Residue indices
                            # elem[2] --> evidence lines
                            
                            ######################
                            # NOTE: CHANGE THIS IF STATEMENT TO IGNORE OR FOCUS ON SPECIFIC FEATURES
                            if feature[0].upper() not in ['ACTIVE_SITE','ACT_SITE','BINDING']:
                                gene_logger.info(f"        {feature}")
                                continue
                            ######################
                            
                            gene_logger.info(f"        ***{feature}***")
                            feature_string_q = ''
                            feature_string_t = ''
                            resids = list(range(feature[1][0],feature[1][-1]+1))
                            # loop over the resids found in the feature metadata
                            for resid in resids:
                                # check to see if the resid is within the mapping dictionary's keys; if the resid isn't mapped, then ignore it
                                resid_ = str(resid)
                                if resid_ in aln_results['map_1_to_2']:
                                    # grabbing target resids
                                    target_resid = aln_results['map_1_to_2'][resid_][0]

                                    feature_string_q += f'{resid_} '
                                    feature_string_t += f'{target_resid} '

                                    # check to see if the two residues have the same resname
                                    if aln_results['map_1_to_2'][resid_][1] == aln_results['map_1_to_2'][resid_][2]:
                                        hit_mapping_conservation.append(target_resid)
                                    
                                    # add to counter for each target_resid; to be used for drawing spheres representing the frequency of this resid in "feature space"
                                    # a single resid can be associated with multiple features
                                    af_model_features_dict[target_resid]['feat_counts'] += 1
                                    # add feature's information to the target_resid string list to be printed to the log file.
                                    af_model_features_dict[target_resid]['feat_string'] = af_model_features_dict[target_resid].get('feat_string') + f'###{hit_pdbid_chainid} {resid_} resid is listed as {feature}\n'
                                    # add feature's type info to the target_resid list of types
                                    af_model_features_dict[target_resid]['feat_types'].append(feature[0])

                            
                            # check to see if the string is not empty; if it is, then no representation will be written for this feature
                            if feature_string_q:
                                # prep the representation for the feature, query structure
                                hit_string += f'### {feature}\n'
                                hit_string += 'mol representation Licorice 0.100000 90.000000 90.000000\n' 
                                hit_string += 'mol color Name\n' 
                                hit_string += 'mol selection {resid ' + feature_string_q +'}\n'
                                hit_string += 'mol material FigRBD\n' 
                                hit_string += 'mol addrep top\n' 
                except Exception as e:
                    print(traceback.format_exc())

            # hit structure does not have a Uniprot ID so no need to dig for features
            else:
                hit_string = f'### {hit_pdbid_chainid} ###\n' + hit_string
           
            hit_mapping_conservation = list(set(hit_mapping_conservation))
            for resid in hit_mapping_conservation:
                af_model_features_dict[resid]['cons_counts'] += 1

            # hide the hit's molecular representations to expedite loading the structures in VMD
            hit_string += 'mol off top\n'
            # write the hit's lines to the vmd vis state file
            vmd_out.write(hit_string)

        gene_logger.info('\n------------------------------------------------------------------------------------------------')
        gene_logger.info(f"Done parsing the alignment hit metadata.\n\n")

        # get a list of resids and their respective feat_counts values, if they are observed at least once
        features_list = [[resid,af_model_features_dict[resid]['feat_counts']] for resid in af_model_features_dict if af_model_features_dict[resid]['feat_counts'] != 0]
        # rank the features_list by the 'feat_counts' values
        features_list = sorted(features_list, key=lambda x: x[1])
        
        rep_count = 1
        # loop over residues in the features_list
        for resid_ in features_list[::-1]:
            resid = resid_[0]
            # if the residue is never observed in a feature, skip it
            if af_model_features_dict[resid]['feat_counts'] == 0:
                continue
            
            # for each residue in a feature space element:
            # a set of comment lines will be written out holding information about the hit's feature types
            # one representation will be written that has the feature_type(s) in the atom_selection string; will be hidden
            # one representation will be written that has the resid in the atom_selection string
           
            gene_logger.info('\n------------------------------------------------------------------------------------------------')
            # writing feature comments out to human-readable log file
            gene_logger.info(f"Model resid {resid} - feature counts: {af_model_features_dict[resid]['feat_counts']}, conserved counts: {af_model_features_dict[resid]['cons_counts']}\n{af_model_features_dict[resid]['feat_string']}")

            # add comment lines to the vis state file listing the RESID meta-data 
            af_model_string += f"### RESID {resid} - feature counts: {af_model_features_dict[resid]['feat_counts']}, conserved counts: {af_model_features_dict[resid]['cons_counts']}\n"
            af_model_string += af_model_features_dict[resid]['feat_string']

            # get list of unique feature types
            feat_types = list(set(af_model_features_dict[resid]['feat_types']))
            # add molecular representation that should always remain hidden
            af_model_string += 'mol representation Licorice 0.000000 1.000000 1.000000\n'
            af_model_string += 'mol selection {type FeatCounts ' + str(af_model_features_dict[resid]['feat_counts']) + ' ConsCounts ' + str(af_model_features_dict[resid]['cons_counts']) + ' ' + ' '.join(feat_types) + '}\n'
            af_model_string += 'mol addrep top\n' 
            af_model_string += f'mol showrep top {rep_count} off\n' 
            rep_count += 1
            # add molecular representation highlighting the residue
            af_model_string += 'mol representation Licorice 0.100000 90.000000 90.000000\n' 
            af_model_string += 'mol color Name\n' 
            af_model_string += 'mol selection {resid ' + resid +'}\n'
            af_model_string += 'mol material FigRBD\n' 
            af_model_string += 'mol addrep top\n' 
            af_model_string += f'mol showrep top {rep_count} off\n' 
            rep_count += 1


        # ----------------------------------------
        # setting up the graphics molecules, spheres showing metrics of residues
        # ----------------------------------------
        # writing the colorbars out to vis state file
        vmd_out.write('# setting colorid rgb values\n')
        for i in all_colorids:
            vmd_out.write(color_defs[i])
      
        vmd_out.write('mol new\n # setting up feature_counts molecule\n')
        
        # gather and rank the feature counts results
        if len(features_list) > 0:
            feature_max_counts = np.max([elem[1] for elem in features_list])
            color_delta = feature_max_counts/(nFeatureColorids-1)

            # drawing spheres of the feature counts, normalized to the desired range.
            for elem in features_list[::-1]:
                resid = elem[0]
                counts = elem[1]
                coord = af_model_resid_coord[resid]
                radius = min_radius + ((counts/feature_max_counts)*radius_diff)
                color = feature_colorids[0] + int(counts/color_delta)
                vmd_out.write(f'# resid {resid} with {counts} feature counts\ngraphics top color {color}\n' + 'draw sphere {' + str(coord[0]) + ' ' + str(coord[1]) + ' ' + str(coord[2]) + '} resolution 85 radius ' + str(radius) + '\n')
            vmd_out.write('mol delrep 0 top\nmol rename top feature_counts\n')

            vmd_out.write('mol new\n # setting up conservation_counts molecule\n')
        
            cons_counts = [[resid,af_model_features_dict[resid]['cons_counts']] for resid in af_model_features_dict if af_model_features_dict[resid]['cons_counts'] != 0]
            if len(cons_counts) > 0:
                # gather and rank the conservation counts results
                cons_counts = sorted(cons_counts,key=lambda x: x[1])
                cons_max_counts = np.max([elem[1] for elem in cons_counts])
                color_delta = cons_max_counts/(nConsColorids-1)

                # drawing spheres of the conservation counts, normalized to the desired range.
                for elem in cons_counts[::-1]:
                    resid = elem[0]
                    counts = elem[1]
                    coord = af_model_resid_coord[resid]
                    radius = min_radius + ((counts/cons_max_counts)*radius_diff)
                    color = cons_colorids[0] + int(counts/color_delta)
                    vmd_out.write(f'# resid {resid} with {counts} conservation counts\ngraphics top color {color}\n' + 'draw sphere {' + str(coord[0]) + ' ' + str(coord[1]) + ' ' + str(coord[2]) + '} resolution 85 radius ' + str(radius) + '\n')
                vmd_out.write('mol delrep 0 top\nmol rename top cons_counts\n')
        
        # ----------------------------------------
        # writing the AF model's molecule to the vis state file
        # ----------------------------------------
        vmd_out.write(af_model_string)
        gene_logger.info('\n------------------------------------------------------------------------------------------------')
        gene_logger.info(f"Done writing the VMD vis-state file.")
    
    # FEATURES - plot colorbar 
    cmap = mpl.colors.LinearSegmentedColormap('my_cmap',feature_cdict,nFeatureColorids)
    fig, ax = plt.subplots(figsize=(2,8))
    fig.subplots_adjust(right=0.5)
    norm = mpl.colors.Normalize(vmin=1,vmax=feature_max_counts)
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,spacing='uniform',norm=norm,orientation='vertical',ticks=[1,feature_max_counts])
    cb.set_label('Feature Counts',size=16)
    cb.set_ticklabels(['1',str(feature_max_counts)])
    plt.savefig(f'{pwd}/{gene}/feature_color_bar.png',dpi=600,transparent=True)
    plt.close()
    gene_logger.info(f"Done making the colorbar for the feature_counts metric.")

    # FEATURES - plot radius and colorbar legend...
    fig = plt.figure(dpi=600)
    ax = plt.gca()
    # gotta actually plot something on the plot to force the equal aspect apparently...
    plt.plot([-max_radius,max_radius],[0,0],'w-')
    plt.plot([0,0],[-max_radius,max_radius],'w-')
    plt.xlim((-max_radius-x_space,max_radius+x_space))
    plt.ylim((-max_radius-x_space,max_radius+x_space))
    ax.set_aspect('equal')
    
    # get the bounding box object associated with the axis
    axbb = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # axbb.height and axbb.width should very nearly be equal, but stil, grab the maximum value of the pair
    axbb_dim = np.max([axbb.height,axbb.width]) # dimension in inches 
    width = 2*max_radius    # dimension in data
    px_dim = axbb_dim * fig.dpi * width / (ax.get_xlim()[1] - ax.get_xlim()[0])
    px_int = math.ceil(px_dim)
    # prep the canvas to be drawn on
    canvas = PIL.Image.new('RGBA',(px_int,px_int),bg)
    draw = PIL.ImageDraw.Draw(canvas)
    # loop over possible colorids
    # do this in reverse because cdict is organized small to large counts/radiuses
    # but I need to draw the circles large to small
    radiuses = np.linspace(min_radius,max_radius,nFeatureColorids)
    for colorid in range(nFeatureColorids)[::-1]:
        color = (math.ceil(feature_cdict['red'][colorid][1]*255),
                 math.ceil(feature_cdict['green'][colorid][1]*255),
                 math.ceil(feature_cdict['blue'][colorid][1]*255))
        # relative radius to use for defining pixel positions
        rel_radius_dif = (1 - radiuses[colorid]/max_radius)/2
        start_pixel = math.ceil(rel_radius_dif*px_int)
        end_pixel   = math.ceil((1-rel_radius_dif)*px_int)
        draw.ellipse((start_pixel,start_pixel,end_pixel,end_pixel),fill=color,outline=color)
    
    # draw circle assocaited with radiuses below the minimum radius; color with white
    final_radius = min_radius - (radiuses[1]-radiuses[0])
    rel_radius_dif = (1 - final_radius/max_radius)/2
    color = (255,255,255)
    start_pixel = math.ceil(rel_radius_dif*px_int)
    end_pixel   = math.ceil((1-rel_radius_dif)*px_int)
    draw.ellipse((start_pixel,start_pixel,end_pixel,end_pixel),fill=color,outline=color)

    # imshow to put the pillow numpy array on the plot
    ax.imshow(np.asarray(canvas), 
              extent = (-max_radius, max_radius, -max_radius, max_radius), 
              aspect = 'equal', 
              interpolation = 'lanczos',
              zorder= 3)

    plt.xlim((0,max_radius+x_space))
    plt.ylim((0,max_radius+x_space))

    plt.xlabel(r'Radius ($\AA$)',size=16)
    plt.xticks(ticks=[0,min_radius,max_radius],labels=[0,min_radius,max_radius])
    
    plt.ylabel('Feature Counts',size=16)
    plt.yticks(ticks=[0,min_radius,max_radius],labels=[0,1,feature_max_counts])
    
    plt.tight_layout()
    plt.savefig(f'{pwd}/{gene}/feature_radius_colorbar.png',dpi=600,transparent=True)
    plt.close()
    gene_logger.info(f"Done making the plot for the feature_counts metric.")

    # CONSERVATION - plot colorbar 
    cmap = mpl.colors.LinearSegmentedColormap('my_cmap',cons_cdict,nConsColorids)
    fig, ax = plt.subplots(figsize=(2,8))
    fig.subplots_adjust(right=0.5)
    norm = mpl.colors.Normalize(vmin=1,vmax=cons_max_counts)
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,spacing='uniform',norm=norm,orientation='vertical',ticks=[1,cons_max_counts])
    cb.set_label('Conservation Counts',size=16)
    cb.set_ticklabels(['1',str(cons_max_counts)])
    plt.savefig(f'{pwd}/{gene}/conservation_color_bar.png',dpi=600,transparent=True)
    plt.close()
    gene_logger.info(f"Done making the colorbar for the cons_counts metric.")

    # CONSERVATION -plot radius and colorbar legend...
    fig = plt.figure(dpi=600)
    ax = plt.gca()
    # gotta actually plot something on the plot to force the equal aspect apparently...
    plt.plot([-max_radius,max_radius],[0,0],'w-')
    plt.plot([0,0],[-max_radius,max_radius],'w-')
    plt.xlim((-max_radius-x_space,max_radius+x_space))
    plt.ylim((-max_radius-x_space,max_radius+x_space))
    ax.set_aspect('equal')
    
    # get the bounding box object associated with the axis
    axbb = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # axbb.height and axbb.width should very nearly be equal, but stil, grab the maximum value of the pair
    axbb_dim = np.max([axbb.height,axbb.width]) # dimension in inches 
    width = 2*max_radius    # dimension in data
    px_dim = axbb_dim * fig.dpi * width / (ax.get_xlim()[1] - ax.get_xlim()[0])
    px_int = math.ceil(px_dim)
    # prep the canvas to be drawn on
    canvas = PIL.Image.new('RGBA',(px_int,px_int),bg)
    draw = PIL.ImageDraw.Draw(canvas)
    # loop over possible colorids
    # do this in reverse because cdict is organized small to large counts/radiuses
    # but I need to draw the circles large to small
    radiuses = np.linspace(min_radius,max_radius,nConsColorids)
    for colorid in range(nConsColorids)[::-1]:
        color = (math.ceil(cons_cdict['red'][colorid][1]*255),
                 math.ceil(cons_cdict['green'][colorid][1]*255),
                 math.ceil(cons_cdict['blue'][colorid][1]*255))
        # relative radius to use for defining pixel positions
        rel_radius_dif = (1 - radiuses[colorid]/max_radius)/2
        start_pixel = math.ceil(rel_radius_dif*px_int)
        end_pixel   = math.ceil((1-rel_radius_dif)*px_int)
        #print((start_pixel,start_pixel,end_pixel,end_pixel),color)
        draw.ellipse((start_pixel,start_pixel,end_pixel,end_pixel),fill=color,outline=color)
    
    # draw circle assocaited with radiuses below the minimum radius; color with white
    final_radius = min_radius - (radiuses[1]-radiuses[0])
    rel_radius_dif = (1 - final_radius/max_radius)/2
    color = (255,255,255)
    start_pixel = math.ceil(rel_radius_dif*px_int)
    end_pixel   = math.ceil((1-rel_radius_dif)*px_int)
    draw.ellipse((start_pixel,start_pixel,end_pixel,end_pixel),fill=color,outline=color)

    # imshow to put the pillow numpy array on the plot
    ax.imshow(np.asarray(canvas), 
              extent = (-max_radius, max_radius, -max_radius, max_radius), 
              aspect = 'equal', 
              interpolation = 'lanczos',
              zorder= 3)

    plt.xlim((0,max_radius+x_space))
    plt.ylim((0,max_radius+x_space))

    plt.xlabel(r'Radius ($\AA$)',size=16)
    plt.xticks(ticks=[0,min_radius,max_radius],labels=[0,min_radius,max_radius])
    
    plt.ylabel('Conservation Counts',size=16)
    plt.yticks(ticks=[0,min_radius,max_radius],labels=[0,1,cons_max_counts])
    
    plt.tight_layout()
    plt.savefig(f'{pwd}/{gene}/conservation_radius_colorbar.png',dpi=600,transparent=True)
    plt.close()
    gene_logger.info(f"Done making the plot for the cons_counts metric.")
    gene_logger.info(f"Finished. Took {time.time() - gene_st} seconds.")
    clean_logger(gene_logger)

main_logger.info(f"Finished creating all of the vis state file and plots for each gene. Closing out after {time.time() - start_time} seconds.")
clean_logger(main_logger)

