
import sys
import pickle
import numpy as np
import pandas
import ast
import math
import matplotlib.pyplot as plt
import PIL
from PIL import Image, ImageDraw, ImageFont, ImageColor
import wordcloud
from pathlib import Path

# sys.argv[1] - alignment results csv file
# sys.argv[2] - pdb-to-uniprot mapping dict
# sys.argv[3] - uniprot metadata dict
# sys.argv[4] - pickle file containing parsed pdb70 ecid counts
# sys.argv[5] - pickle file containing parsed brenda ecid lists
# sys.argv[6] - avgTMscore cutoff
# sys.argv[7] - pid string
# sys.argv[8] - bootstrap iterations


###################
# PLOTTING VARIABLES
###################

colors = {'1': '#f00',      # pure red
          '2': '#0f0',      # pure green
          '3': '#00f',      # pure blue
          '4': '#ff8000',   # pure orange   # cause yellow is bad
          '5': '#0ff',      # pure cyan
          '6': '#8000ff',   # pure purple 
          '7': '#f0f'}      # pure pink


###################
# MATH FUNCTIONS
###################

def calc_rel_entropy_component(posterior_prop,prior_prop):
    """takes two probabilities, one from the posterior and the prior, and 
    calculate the Kullback-Leibler divergence component. 
    
    MATH:
    D_{KL}(P \parrallel Q) = \sum_{x \in X} P(x) \log_{2} \left( \frac{P(x)}{Q(x)} \right)
    where Q is the prior distribution, P is the posterior distribution, and x 
    is a single element in the full sample space, X, within which both P and Q 
    sample. 
    This function does the intra-summation calculation for a value of x. 

    INPUT:
        posterior_prop: float, probability value
        prior_prop: float, probability value

    RETURNS:
        KLD component for these two probability values
    """
    return posterior_prop * np.log2( posterior_prop / prior_prop )


def calc_prop(counts_dict,pseudocount=1):
    """take a counts dictionary and convert into a probability dictionary,
    accounting for pseudocounts if defined.

    INPUT:
        counts_dict: dict, keys are any qualitative classifier, values are the 
                     count values for observing said classifier
        pseudocount: int/float, added to all counts during the probability 
                     calcualtion; necessary to implement a pseudocount to avoid
                     zeros
    """
    prop_dict = {}
    nRealCounts = np.sum([counts_dict[key] for key in counts_dict])
    nCounts = nRealCounts + len(counts_dict)*pseudocount
    for key,value in counts_dict.items():
        prop_dict[key] =  (value + pseudocount)/(nCounts)
    return prop_dict


###################
# PLOTTING FUNCTIONS
###################

def text_draw_np(text, width, height, *, font = '/home/russ/Apps/anaconda3/pkgs/font-ttf-ubuntu-0.83-hab24e00_0/fonts/UbuntuMono-B.ttf', bg = (255, 255, 255, 0), color = (0, 0, 0, 255), remove_gaps = False, cache = {}):
    """
    modified from https://stackoverflow.com/questions/65609379/is-there-a-way-to-plot-text-with-specified-width-and-height-using-matplotlib
   
    all fields of rgba color definitions range from 0 to 255, including the opacity value. duh

    """

    def get_font(fname,size):
        """
        """
        key = ('font',fname,size)
        if key not in cache:
            cache[key] = PIL.ImageFont.truetype(fname, size = size, encoding = 'unic')
        return cache[key]

    # by default, bg is defined as the rgba tuple for transparent white (len of 4).
    # but if bg is defined as a string, then this will convert to a rgb tuple (len of 3) assuming full opacity.
    if type(bg) != tuple:
        bg = PIL.ImageColor.getrgb(bg)
    
    # by default, color is defined as the rgba tuple for opaque black (len of 4).
    # but if color is defined as a string, then this will convert to a rgb tuple (len of 3) assuming full opacity.
    if type(color) != tuple:
        color = PIL.ImageColor.getrgb(color)
    # if color is a tuple and has a len of 4, need to make sure that the transparency value != 0. What would be the point if it was?
    elif type(color) == tuple and len(color) == 4:
        assert color[3] != 0
   
    # check to make sure the user didn't define the same rgb(a) values for the bg and color variables. What would be the point if they were?
    if len(color) == len(bg):
        assert color != bg
    else:
        assert color[:3] != bg[:3]

    width, height = math.ceil(width), math.ceil(height)
    pil_font = get_font(font, 24)   # grab desired font with fs of 24
    text_width, text_height = pil_font.getsize(text)    # grab height and width of text for known fs    # getsize is throwing a depreciation warning
    #bbox = pil_font.getbbox(text)
    #text_width = bbox[2] - bbox[0]
    #text_height= bbox[3] - bbox[1]
    #print(bbox,text_width,text_height)
    pil_font = get_font(font, math.ceil(1.2 * 24 * max(width / text_width, height / text_height)))  # calculate the desired fontsize and get the text in desired font/scaling
    #bbox = pil_font.getbbox(text)
    #text_width = bbox[2] - bbox[0]
    #text_height= bbox[3] - bbox[1]
    #print(bbox,text_width,text_height)
    text_width, text_height = pil_font.getsize(text)    # grab height and width of text in new font # getsize is throwing a depreciation warning
    # prep the pillow image canvas, allow for transparency if desired
    if len(bg) == 4 or len(color) == 4:
        canvas = PIL.Image.new('RGBA', (text_width, text_height), bg)
    else:
        canvas = PIL.Image.new('RGB',  (text_width, text_height), bg)

    draw = PIL.ImageDraw.Draw(canvas)
    draw.text((0, 0), text, font = pil_font, fill = color)

    # removing excess white space on top and bottom of the text; 
    # removing excess white space on left and right currently commented out
    if remove_gaps: # and len(bg) == 3:
        a = np.asarray(canvas)  # canvas' rgb values in Height x Width x 3 matrix; include foreground and background colors
        # copy shape of the canvas (MxNx 3 or 4)
        b = np.zeros_like(a)
        # fill b's elements with the background color's rgb values
        b[:, :, 0] = bg[0]
        b[:, :, 1] = bg[1]
        b[:, :, 2] = bg[2]
        # check to see if canvas is "rgba" and `bg` has a transparency value
        if a.shape[-1] == 4 and len(bg) == 4:
            b[:, :, 3] = bg[3]
        # if canvas is "rgba" but `bg` is an rgb tuple, then `color` is rgba tuple
        # and we need to fill in the `b` array with the `a` array's transparency values
        elif a.shape[-1] == 4:
            b[:, :, 3] = a[:, :, 3]

        # (a != b) will create a 3d boolean array, True for any color element that does not match between the bg array and the canvas
        # .reshape(a.shape[0],-1) flattens the 3rd dimension, creating a MxN*3 (or 4) array
        # np.any(...,axis=-1) flattens the 2d bool array down to a 1d bool array of shape M (height)
        rows = np.any((a != b).reshape(a.shape[0], -1), axis = -1)  # for height
        ## .transpose(1, 0 , 2) switches axes 0 and 1 of the boolean array, now NxMx3 (or 4) array
        ## .reshape(a.shape[1],-1) flattens the 3rd dimension, creating a NxM*3 (or 4) array
        ## np.any(...,axis=-1) flattens the 2d bool array down to a 1d bool array of shape N (width)
        #columns = np.any((a != b).transpose(1, 0, 2).reshape(a.shape[1], -1), axis=-1)  # for width
        
        top = np.flatnonzero(rows)[0]
        bot = np.flatnonzero(rows)[-1]
        #lef = np.flatnonzero(columns)[0]
        #rig = np.flatnonzero(columns)[-1]
        # remove both top/bottom as well as right/left "gaps"
        #a = a[top : bot, lef : rig]
        # only remove top/bottom "gaps"
        a = a[top : bot]
        # convert the `a` array back to a canvas object
        canvas = PIL.Image.fromarray(a)
    # truly ensure the exact width and height by resizing canvas explicitly. 
    # PIL.Image.LANCZOS is a high quality filterfor doing this
    canvas = canvas.resize((width, height), PIL.Image.LANCZOS)
    return np.asarray(canvas)


def text_draw_mpl(fig, ax, text, offset_x, offset_y, width, height, zorder = 3, **nargs):
    axbb = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    pxw, pxh = axbb.width * fig.dpi * width / (ax.get_xlim()[1] - ax.get_xlim()[0]), axbb.height * fig.dpi * height / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.imshow(text_draw_np(text, pxw * 1.2, pxh * 1.2, **nargs), extent = (offset_x, offset_x + width, offset_y, offset_y + height), aspect = 'auto', interpolation = 'lanczos',zorder= zorder)


def get_ecid_color(word, font_size, position, orientation, font_path, random_state=None, **kwargs):
    '''
    '''
    temp = word.split('.')
    starting_color = PIL.ImageColor.getrgb(colors[temp[0]])
    
    M = np.max(starting_color)
    m = np.min(starting_color)
    # calculate the lightness, should have range of 0 to 1
    L = (M+m)/510
    
    # if all three rgb values are equal, problems arise...
    # d will equal zero, so the color will always turn out gray (whether desired or not)
    # H will have a nan denominator value causing problems
    # if L == 1, color will always be white, no further calcs needed (also avoids S calc error)
    # if L == 0, color will always be black, no further calcs needed
    if M == m or L == 1 or L == 0:
        H = 0
        S = 0
        return f"hsl({H}, {S}%, {L}%)"
    else:
        d = (M-m)/255
        
        # calculate saturation
        S = d/(1-np.abs(2*L-1))
        # S > 1 if L is very close to 1
        if S > 1:
            S = 1.
        
        # calculate Hue
        numerator  = (starting_color[0] - 0.5*starting_color[1] - 0.5*starting_color[2])
        denomenator = np.sqrt(starting_color[0]**2 + starting_color[1]**2 + starting_color[2]**2 - starting_color[0]*starting_color[1] - starting_color[0]*starting_color[2] - starting_color[1]*starting_color[2])
        element = numerator/denomenator
        if starting_color[1] >= starting_color[2]:
            H = int(180*np.arccos(element)/np.pi)
        else:
            H = int(360 - 180*np.arccos(element)/np.pi)

    L *= 100
    S *= 100
    L = int(L)
    S = int(S)

    if len(temp) ==4 and temp[-1] != '*':
        # add or subtract a random value from the Hue value
        H += np.random.randint(-5,6)
        # H cannot be negative, take the modulus with 360
        H = H % 360
        # set Saturation randomly between 33% and 100%; don't want things to get too gray
        S = np.random.randint(50,101)
        # set Lightness randomly between 25% and 75%; don't want things to get too black or white
        L = np.random.randint(25,76)

    return f"hsl({H}, {S}%, {L}%)"


def draw_logo(counts_dict,total_rel_entropy,figure_name):
    """
    """
    nCounts = np.sum([counts_dict[key] for key in counts_dict])
    logo_heights = [[ecid, total_rel_entropy*counts_dict[ecid]/nCounts] for ecid in counts_dict if counts_dict[ecid] != 0]
    logo_heights = sorted(logo_heights, key = lambda x: x[1])

    fig = plt.figure(dpi=600)
    ax = plt.gca()
    x_space = 0.01
    y_space = 0.01*total_rel_entropy
    plt.plot([0-x_space,1+x_space],[total_rel_entropy,total_rel_entropy],'b--',zorder=1)
    plt.xlim((0-x_space,1+x_space))
    plt.ylim((0-y_space,total_rel_entropy+y_space))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(plt.NullLocator())
    plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--',zorder=1)
    max_length = np.max([len(elem[0]) for elem in logo_heights])
    offset_y = 0
    offset_x = 0
    for i, ecid in enumerate(logo_heights):
        text_draw_mpl(fig, ax, f' {ecid[0]:<{max_length}} ', offset_x, offset_y, 1, ecid[1], color = colors[ecid[0].split('.')[0]], remove_gaps = True, zorder= i+1) #bg = 'white', 
        offset_y += ecid[1]
    
    plt.ylabel('Relative Entropy\n(bits)',size=16)
    plt.tight_layout()
    plt.savefig(figure_name,dpi=600,transparent=True)
    plt.close()


###################
# MAIN
###################

# reading structure alignment results csv
struct_aln_results = pandas.read_csv(sys.argv[1],usecols = ['Target Path','avgTMscore'])

# pdb2uniprot ID mapping
with open(sys.argv[2],'rb') as in_pkl:
    pdb2uni = pickle.load(in_pkl)

# Uniprot Metadata dictionary
with open(sys.argv[3],'rb') as in_pkl:
    uniprot_dict = pickle.load(in_pkl)

# reading the pdb70 base distribution of EC IDs
with open(sys.argv[4],'rb') as pkl_file:
    parsed_pdb70_ecid_results = pickle.load(pkl_file)
pdb70_active_ecids, pdb70_all_ecids, pdb70_ecid_reps = parsed_pdb70_ecid_results[5:8]

# calculate the prior distributions of active and all ECIDs
all_prior    = {}
for ecid in pdb70_all_ecids:
    all_prior[ecid] = pdb70_all_ecids[ecid]/pdb70_all_ecids['Total']

active_prior = {}
for ecid in pdb70_active_ecids:
    active_prior[ecid] = pdb70_active_ecids[ecid]/pdb70_active_ecids['Total']

# reading the BRENDA list of active/all EC IDs for checking on alignment hit 
# labeled EC IDs
with open(sys.argv[5],'rb') as pkl_file:
    parsed_brenda = pickle.load(pkl_file)
active_ecids, defunct_ecids = parsed_brenda[2:4]

# other input variables
avgTMscore_cutoff = float(sys.argv[6])
bs_iterations = int(sys.argv[8])

# parse the aln hits dataframe to pull PDBID_CHAINID strings, the respective 
# uniprot accession ID, and finally the associated EC ID(s), if present.
# EC IDs are then checked against the BRENDA list of active and defunct IDs.
# If the EC ID is defunct, check if there's a more recent ID to be used. If
# there is, use the more up-to-date ID, else gather the old ID but don't label
# it as "active", just categorize it in the "all". 
# If the original EC ID is incomplete, categorize it in the "all" list as well.
# 
nUniprots = 0
active_nECIDs = 0
all_nECIDs = {1:0,2:0,3:0,4:0}
active_counts = {}
all_counts    = {}
for hit in struct_aln_results.iloc:
    if hit['avgTMscore'] < avgTMscore_cutoff:
        print(f"{Path(hit['Target Path']).stem} has an avgTMscore less than the {avgTMscore_cutoff}, breaking out of hits loop")
        break
    
    # grab UniProt ID associated with alignment hit
    try:
        uniprotid = pdb2uni[Path(hit['Target Path']).stem]
    except Exception as e:
        print(f"{Path(hit['Target Path']).stem} has no uniprot ID or experienced some other error: {e}")
        continue
    # if uniprotid == None, skip to the next hit
    if not uniprotid:
        print(f"{Path(hit['Target Path']).stem} has no uniprot ID")
        continue 
    
    nUniprots += 1

    # parse the ecid string(s)
    try:
        ecids = uniprot_dict[uniprotid]['ecIDs']
        #ecids = [line.split()[1].split('EC=')[1].strip(',;:') for line in uniprot_dict[uniprotid]['ecIDs'] if line[:2] == 'DE']
    except Exception as e:
        print(f"{Path(hit['Target Path']).stem}, {uniprotid}, Error: {e}, error in parsing the ECID string(s)")
        continue
    
    # loop over the parsed ecid strings to gather relevant ecid prefix counts, etc
    for ecid in ecids:
        # check if the labeled EC ID is "active" in the BRENDA list; initially
        # assuming false
        active_bool = False
        if ecid not in active_ecids and ecid in defunct_ecids:
            # defunct_ecids[ecid] should be a list of strings; but it can also
            # be an empty list. Check to see if its an empty list. active_bool
            # will remain false
            if not defunct_ecids[ecid]:
                print(f"{Path(hit['Target Path']).stem}: {ecid} is listed in the BRENDA EC ID table as having been 'deleted' and has no replacement.")
            # it isn't empty, good. Now, loop over the elements in the list to
            # find the first "active" EC ID, if there is one. Else, active_bool 
            # will remain false and the original defunct EC ID will be passed 
            # on for the "all" category.
            else:
                for replacement in defunct_ecids[ecid]:
                    if replacement in active_ecids:
                        print(f"{Path(hit['Target Path']).stem}: {ecid} is listed in the BRENDA EC ID table as having been 'changed' to {defunct_ecids[ecid]}. Will use EC {replacement} in the logo figure.")
                        active_bool = True
                        ecid = replacement
                        break
        # if ecid not in either active_ecids nor defunct_ecids (assumed by the 
        # previous boolean test); active_bool will remain false
        elif ecid not in active_ecids:
            print(f"{Path(hit['Target Path']).stem}: {ecid} is not found within the active EC ID table. Likely an incomplete EC ID.")
        # ecid is in active ecids list (assumed by the two previous boolean tests)
        else:
            active_bool = True

        # parse the ecid
        digits = ecid.split('.')
        # ignore digits if they are '-'; keep preliminary labels
        digit_bools = [digit != '-' for digit in digits]
        # if there are no "good" digits, skip this ecid
        if np.sum(digit_bools) == 0:
            print(f"{Path(hit['Target Path']).stem}: {ecid} has no parsable digits. Skipping.")
            continue

        ############################
        # DEALING WITH ALL ECIDS
        # loop over the first three digits, creating prefices if they do not
        # include '-'; add one count to the respective dict key
        for i in range(3):
            # check the bool value
            # NOTE: USED TO HAVE A TRY EXCEPT BLOCK HERE... NOT SURE WHY THOUGH
            if digit_bools[i]:
                # create the prefix string
                digit_string = '.'.join(digits[:i+1])+'.*'
                # add 1 to the all_counts array for each prefix
                all_counts[digit_string] = all_counts.get(digit_string,0) + 1
                all_nECIDs[i+1] += 1
        
        # add 1 to the all_counts array for the ecid
        if digit_bools[-1]:
            all_counts[ecid] = all_counts.get(ecid,0) + 1
            all_nECIDs[4] += 1

        ############################

        ############################
        # DEALING WITH ONLY ACTIVE ECIDS
        # loop over the first three digits, creating prefices add one count
        # to the respective dict key
        if active_bool:
            active_nECIDs += 1
            # add 1 to the active_counts array for the full ecid
            active_counts[ecid] = active_counts.get(ecid,0) + 1
            for i in range(3):
                # create the prefix string
                digit_string = '.'.join(digits[:i+1])+'.*'
                # add 1 to the active_counts array for each prefix
                active_counts[digit_string] = active_counts.get(digit_string,0) + 1
        ############################


print(f'For this protein, {active_nECIDs} active, full ECIDs were found from the alignment results:')
print(f'{active_counts}')
print(f'Considering incomplete/defunct ECIDs, the number ECIDs were found from the alignment results:\n{all_nECIDs[1]} with at least one digit,\n{all_nECIDs[2]} with at least two digits,\n{all_nECIDs[3]} with at least three digits,\n{all_nECIDs[4]} full ECIDs.')
print(f'{all_counts}\n\n')

# no need to go into plotting if the all_counts dictionary is completely empty
if not all_counts:
    exit()

###################
# RELATIVE ENTROPY MATH
###################

# Relative Entropy
# total height: D_{KL}= $\Sigma_{n=1}^{N} P_{n} log2(\frac(P_{n})(Q_{n})); 
#   where n is one EC ID in the set of N EC IDs within the PDB70 structure library's EC IDlist,
#   Q_{n} is n's proportion in the PDB70 structure library's EC ID distribution (prior),
#   and P_{n} is the observed proportion of n in the structural alignment analysis (posterior).
# To avoid numerical errors, a pseudocount must be added to either distribution if the other samples a new EC ID (not possible for Q_{n}); where P_{n} = 0, add a pseudocount
#
# Relative height of each stack element: total height * $p_{n}$
# width: equal to 1 for the time being

# convert counts to smoothed proportion;
# theta_i = (x_i + alpha)/(N + alpha*d)
#       theta_i is Laplace smoothed posterior probability for EC ID i
#       x_i is counts of observed EC ID i
#       alpha is the pseudocount value
#       N is the total number of observed EC IDs
#       d is the total number of possible EC IDs, i
#       the denominator is the total number of counts when pseudocounts are accounted for

###################
# ACTIVE EC ID LOGO FIGURES
###################

# check to see if any EC IDs are in the active_counts dict, if not skip
if active_counts:
    
    ###################
    # prep work for pseudocount applications and bootstrapping
    ###################
    # get the list of full ecids possible in the pdb70 structure library
    full_ecids = [ecid for ecid in pdb70_active_ecids if ecid[-1] != '*' and ecid != 'Total']
    # create a flat list of all possible active ecids, used to randomly pull 
    # ids from during bootstrapping
    active_ecid_list = []
    for ecid in full_ecids:
        active_ecid_list += [ecid]*pdb70_active_ecids[ecid]
    
    # set total pseudocount amount
    total_pseudocount  = 1  # the sum of all applied pseudocounts

    ###################
    # OBSERVED REL ENTROPY CALCULATIONS
    ###################
    obs_rel_entropies = {1:0,2:0,3:0,4:0}
    obs_rel_entropy_components = {1:[],2:[],3:[],4:[]}
    
    ###
    # first digit
    ###
    pseudocount = total_pseudocount / pdb70_ecid_reps['*']
    obs_counts = {}
    prior = {}
    # loop over ecids in the pdb70_ecid_reps dictionary, just using the keys not the values
    for ecid in pdb70_ecid_reps:
        # grab only the first digit prefices
        if len(ecid.split('.')) == 2 and ecid[0] != '*':
            # if prefix in active_counts, grab the active_counts value, else grab 0
            obs_counts[ecid] = active_counts.get(ecid,0)
    true_posterior = calc_prop(obs_counts, pseudocount = pseudocount)
    obs_rel_entropy_components[1] =  [[ecid, calc_rel_entropy_component(true_posterior[ecid],active_prior[ecid])] for ecid in obs_counts]
    obs_rel_entropies[1] = np.sum([elem[1] for elem in obs_rel_entropy_components[1]])
    draw_logo(obs_counts,obs_rel_entropies[1],f'{sys.argv[7]}_active_1.png')

    ###
    # second digit
    ###
    pseudocount = total_pseudocount / pdb70_ecid_reps['*.*']
    obs_counts = {}
    prior = {}
    # loop over ecids in the pdb70_ecid_reps dictionary, just using the keys not the values
    for ecid in pdb70_ecid_reps:
        # grab only the second digit prefices
        if len(ecid.split('.')) == 3 and ecid[0] != '*':
            # if prefix in active_counts, grab the active_counts value, else grab 0
            obs_counts[ecid] = active_counts.get(ecid,0)
    true_posterior = calc_prop(obs_counts, pseudocount = pseudocount)
    obs_rel_entropy_components[2] =  [[ecid, calc_rel_entropy_component(true_posterior[ecid],active_prior[ecid])] for ecid in obs_counts]
    obs_rel_entropies[2] = np.sum([elem[1] for elem in obs_rel_entropy_components[2]])
    draw_logo(obs_counts,obs_rel_entropies[2],f'{sys.argv[7]}_active_2.png')
    
    ###
    # three digit
    ###
    pseudocount = total_pseudocount / pdb70_ecid_reps['*.*.*']
    obs_counts = {}
    prior = {}
    # loop over ecids in the pdb70_ecid_reps dictionary, just using the keys not the values
    for ecid in pdb70_ecid_reps:
        # grab only the third digit prefices
        if len(ecid.split('.')) == 4 and ecid[-1] == '*' and ecid[0] != '*':
            # if prefix in active_counts, grab the active_counts value, else grab 0
            obs_counts[ecid] = active_counts.get(ecid,0)
    true_posterior = calc_prop(obs_counts, pseudocount = pseudocount)
    obs_rel_entropy_components[3] =  [[ecid, calc_rel_entropy_component(true_posterior[ecid],active_prior[ecid])] for ecid in obs_counts]
    obs_rel_entropies[3] = np.sum([elem[1] for elem in obs_rel_entropy_components[3]])
    draw_logo(obs_counts,obs_rel_entropies[3],f'{sys.argv[7]}_active_3.png')

    ###
    # full ECIDs
    ###
    pseudocount = total_pseudocount / pdb70_ecid_reps['*.*.*.*']
    obs_counts = {}
    prior = {}
    # loop over ecids in the pdb70_active_ecids dictionary
    for ecid in pdb70_active_ecids:
        # grab only the full ECIDs
        if len(ecid.split('.')) == 4 and ecid[-1] != '*':
            # if prefix in active_counts, grab the active_counts value, else grab 0
            obs_counts[ecid] = active_counts.get(ecid,0)
    true_posterior = calc_prop(obs_counts, pseudocount = pseudocount)
    obs_rel_entropy_components[4] =  [[ecid, calc_rel_entropy_component(true_posterior[ecid],active_prior[ecid])] for ecid in obs_counts]
    obs_rel_entropies[4] = np.sum([elem[1] for elem in obs_rel_entropy_components[4]])
    draw_logo(obs_counts,obs_rel_entropies[4],f'{sys.argv[7]}_active_4.png')

    ###################
    # BOOTSTRAPPING
    ###################
    
    print(f'Performing {bs_iterations} iterations of bootstrapping, where a random set of full, active EC IDs are gathered (w/o replacement), parsed, and the assocaited relative entropies are calculated.')
    full_bs_klds = {1:[],2:[],3:[],4:[]}
    full_bs_kld_components = {1:[],2:[],3:[],4:[]}
    for i in range(bs_iterations):
        bs_counts = {1:{},2:{},3:{},4:{}}
        # use numpy's default random number generator to get indices of the active_ecid_list
        rng = np.random.default_rng()
        indexes = rng.choice(len(active_ecid_list),size=active_nECIDs,replace=False,shuffle=False)
        bs_ecid_list = [active_ecid_list[idx] for idx in indexes]
        # count the random sampling of EC IDs
        for ecid in bs_ecid_list:
            bs_counts[4][ecid] = bs_counts[4].get(ecid,0) + 1
            digits = ecid.split('.')
            for i in range(3):
                digit_string = '.'.join(digits[:i+1])+'.*'
                bs_counts[i+1][digit_string] = bs_counts[i+1].get(digit_string,0) + 1
        # fill in missing ecids with zeros
        for ecid in pdb70_active_ecids:
            digits = ecid.split('.')
            if len(digits) == 2:
                bs_counts[1][ecid] = bs_counts[1].get(ecid,0)
            elif len(digits) == 3:
                bs_counts[2][ecid] = bs_counts[2].get(ecid,0)
            elif len(digits) == 4 and digits[-1] == '*':
                bs_counts[3][ecid] = bs_counts[3].get(ecid,0)
            elif len(digits) == 4 and digits[-1] != '*':
                bs_counts[4][ecid] = bs_counts[4].get(ecid,0)

        # calculate bootstrapped posterior distributions, the rel entropy components, and finally the total rel entropy
        # first digits
        bs_prop = calc_prop(bs_counts[1],pseudocount=total_pseudocount / pdb70_ecid_reps['*'])
        bs_rel_entropies = [calc_rel_entropy_component(bs_prop[ecid],active_prior[ecid]) for ecid in bs_counts[1]]
        bs_rel_entropy = np.sum(bs_rel_entropies)
        full_bs_klds[1].append(bs_rel_entropy)
        full_bs_kld_components[1].append(bs_rel_entropies)
        
        # second digits
        bs_prop = calc_prop(bs_counts[2],pseudocount=total_pseudocount / pdb70_ecid_reps['*.*'])
        bs_rel_entropies = [calc_rel_entropy_component(bs_prop[ecid],active_prior[ecid]) for ecid in bs_counts[2]]
        bs_rel_entropy = np.sum(bs_rel_entropies)
        full_bs_klds[2].append(bs_rel_entropy)
        full_bs_kld_components[2].append(bs_rel_entropies)

        # third digits
        bs_prop = calc_prop(bs_counts[3],pseudocount=total_pseudocount / pdb70_ecid_reps['*.*.*'])
        bs_rel_entropies = [calc_rel_entropy_component(bs_prop[ecid],active_prior[ecid]) for ecid in bs_counts[3]]
        bs_rel_entropy = np.sum(bs_rel_entropies)
        full_bs_klds[3].append(bs_rel_entropy)
        full_bs_kld_components[3].append(bs_rel_entropies)

        # four digits
        bs_prop = calc_prop(bs_counts[4],pseudocount=total_pseudocount / pdb70_ecid_reps['*.*.*.*'])
        bs_rel_entropies = [calc_rel_entropy_component(bs_prop[ecid],active_prior[ecid]) for ecid in bs_counts[4]]
        bs_rel_entropy = np.sum(bs_rel_entropies)
        full_bs_klds[4].append(bs_rel_entropy)
        full_bs_kld_components[4].append(bs_rel_entropies)

    ###################
    # PRINTING SUMMARY
    ###################
    print(f"ACTIVE ECID RELATIVE ENTROPY SUMMARY:\n{active_nECIDs} EC IDs were observed in the alignment hits.")
    print(f"First Digits, Total relative entropy = {obs_rel_entropies[1]} bits.\nBootstrapped Rel Entropy: Mean = {np.mean(full_bs_klds[1])} bits, 100th percentile is {np.max(full_bs_klds[1])} bits.\nBootstrapped Rel Entropy components: Mean = {np.mean(full_bs_kld_components[1])} bits, 100th percentile is {np.max(full_bs_kld_components[1])} bits.\nObserved EC IDs that surpass the 100th percentile:")
    cutoff = np.max(full_bs_kld_components[1])
    for elem in obs_rel_entropy_components[1]:
        if elem[1] > cutoff:
            print(elem[0],elem[1])

    print(f"\nSecond Digits, Total relative entropy = {obs_rel_entropies[2]} bits.\nBootstrapped Rel Entropy: Mean = {np.mean(full_bs_klds[2])} bits, 100th percentile is {np.max(full_bs_klds[2])} bits.\nBootstrapped Rel Entropy components: Mean = {np.mean(full_bs_kld_components[2])} bits, 100th percentile is {np.max(full_bs_kld_components[2])} bits.\nObserved EC IDs that surpass the 100th percentile:")
    cutoff = np.max(full_bs_kld_components[2])
    for elem in obs_rel_entropy_components[2]:
        if elem[1] > cutoff:
            print(elem[0],elem[1])

    print(f"\nThird Digits, Total relative entropy = {obs_rel_entropies[3]} bits.\nBootstrapped Rel Entropy: Mean = {np.mean(full_bs_klds[3])} bits, 100th percentile is {np.max(full_bs_klds[3])} bits.\nBootstrapped Rel Entropy components: Mean = {np.mean(full_bs_kld_components[3])} bits, 100th percentile is {np.max(full_bs_kld_components[3])} bits.\nObserved EC IDs that surpass the 100th percentile:")
    cutoff = np.max(full_bs_kld_components[3])
    for elem in obs_rel_entropy_components[3]:
        if elem[1] > cutoff:
            print(elem[0],elem[1])

    print(f"\nFull ECIDs, Total relative entropy = {obs_rel_entropies[4]} bits.\nBootstrapped Rel Entropy: Mean = {np.mean(full_bs_klds[4])} bits, 100th percentile is {np.max(full_bs_klds[4])} bits.\nBootstrapped Rel Entropy components: Mean = {np.mean(full_bs_kld_components[4])} bits, 100th percentile is {np.max(full_bs_kld_components[4])} bits.\nObserved EC IDs that surpass the 100th percentile:")
    cutoff = np.max(full_bs_kld_components[4])
    for elem in obs_rel_entropy_components[4]:
        if elem[1] > cutoff:
            print(elem[0],elem[1])

    with open(sys.argv[7] + '_active_ecid_data.pkl','wb') as out_pkl:
        pickle_string = 'Organization: this message, observed counts for active EC IDs, Total Relative Entropy values for all completeness of ECIDs, Relative Entropy components for all completenesses of ECIDs, bootstrapped relative entropy values for all completenesses of ECIDS, bootstrapped relative entropy components for all completenesses of ECIDs.'
        pickle.dump([pickle_string,active_counts,obs_rel_entropies,obs_rel_entropy_components,full_bs_klds,full_bs_kld_components],out_pkl)

###################
# ANY/ALL EC ID LOGO FIGURES
###################

### FULL ECIDs

# prep list of relevant ecids possible in the pdb70 structure library
ecids_list  = [ecid for ecid in pdb70_all_ecids if ecid[-1] != '*' and len(ecid.split('.')) == 4]
nECIDs = len(ecids_list)
total_pseudocount = 1
pseudocount = total_pseudocount / nECIDs

# fill counts with zeros or the observed
counts = {}
for ecid in ecids_list:
    counts[ecid] = all_counts.get(ecid,0)
# get number of relevant EC IDs counted, whether defunct or incomplete or active
nCounts = np.sum([counts[ecid] for ecid in ecids_list])

# calculate the posterior probabilities with pseudocounts
posterior = calc_prop(counts,pseudocount=pseudocount)
rel_entropy_components = [[ecid, calc_rel_entropy_component(posterior[ecid],all_prior[ecid])] for ecid in ecids_list]
rel_entropy = np.sum([elem[1] for elem in rel_entropy_components])

# gather relative entropy componenets of actually observed EC IDs
rel_entropies = [elem for elem in rel_entropy_components if counts[elem[0]] != 0]
rel_entropies = sorted(rel_entropies, key=lambda x: x[1])
    
# calculate heights for logo figure
rel_heights = [[ecid, rel_entropy*counts[ecid]/nCounts] for ecid in counts if counts[ecid] != 0]
rel_heights = sorted(rel_heights,key=lambda x: x[1])

print(f'\n\nALL ECIDs:\n{nCounts} full EC IDs were observed, including defunct and active IDs.')
print(f"Total relative entropy = {rel_entropy} bits. Broken down into components:\n")
for elem in rel_entropies:
    print(elem[0],elem[1])

fig = plt.figure(dpi=600)
ax = plt.gca()
x_space = 0.01
y_space = 0.01*rel_entropy
plt.plot([0-x_space,1+x_space],[rel_entropy,rel_entropy],'b--',zorder=1)
plt.xlim((0-x_space,1+x_space))
plt.ylim((0-y_space,rel_entropy+y_space))
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_major_locator(plt.NullLocator())
plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--',zorder=1)
max_length = np.max([len(ecid[0]) for ecid in rel_heights])
offset_y = 0
offset_x = 0
for i, ecid in enumerate(rel_heights):
    text_draw_mpl(fig, ax, f' {ecid[0]:<{max_length}} ', offset_x, offset_y, 1, ecid[1], color = colors[ecid[0].split('.')[0]], remove_gaps = True, zorder= i+1) #bg = 'white', 
    offset_y += ecid[1]

plt.ylabel('Relative Entropy\n(bits)',size=16)
plt.tight_layout()
plt.savefig(f'{sys.argv[7]}_all_4.png',dpi=600,transparent=True)
plt.close()

### Out to 3rd digit

# prep list of relevant ecids possible in the pdb70 structure library
ecids_list  = [ecid for ecid in pdb70_all_ecids if ecid[-1] == '*' and len(ecid.split('.')) == 4]
nECIDs = len(ecids_list)
total_pseudocount = 1
pseudocount = total_pseudocount / nECIDs

# fill counts with zeros or the observed
counts = {}
for ecid in ecids_list:
    counts[ecid] = all_counts.get(ecid,0)
# get number of relevant EC IDs counted, whether defunct or incomplete or active
nCounts = np.sum([counts[ecid] for ecid in ecids_list])

# calculate the posterior probabilities with pseudocounts
posterior = calc_prop(counts,pseudocount=pseudocount)
rel_entropy_components = [[ecid, calc_rel_entropy_component(posterior[ecid],all_prior[ecid])] for ecid in ecids_list]
rel_entropy = np.sum([elem[1] for elem in rel_entropy_components])

# gather relative entropy componenets of actually observed EC IDs
rel_entropies = [elem for elem in rel_entropy_components if counts[elem[0]] != 0]
rel_entropies = sorted(rel_entropies, key=lambda x: x[1])
    
# calculate heights for logo figure
rel_heights = [[ecid, rel_entropy*counts[ecid]/nCounts] for ecid in counts if counts[ecid] != 0]
rel_heights = sorted(rel_heights,key=lambda x: x[1])

print(f'\n{nCounts} third-digit EC IDs were observed, including defunct and active IDs.')
print(f"Total relative entropy = {rel_entropy} bits. Broken down into components:\n")
for elem in rel_entropies:
    print(elem[0],elem[1])

fig = plt.figure(dpi=600)
ax = plt.gca()
x_space = 0.01
y_space = 0.01*rel_entropy
plt.plot([0-x_space,1+x_space],[rel_entropy,rel_entropy],'b--',zorder=1)
plt.xlim((0-x_space,1+x_space))
plt.ylim((0-y_space,rel_entropy+y_space))
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_major_locator(plt.NullLocator())
plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--',zorder=1)
max_length = np.max([len(ecid[0]) for ecid in rel_heights])
offset_y = 0
offset_x = 0
for i, ecid in enumerate(rel_heights):
    text_draw_mpl(fig, ax, f' {ecid[0]:<{max_length}} ', offset_x, offset_y, 1, ecid[1], color = colors[ecid[0].split('.')[0]], remove_gaps = True, zorder= i+1) #bg = 'white', 
    offset_y += ecid[1]

plt.ylabel('Relative Entropy\n(bits)',size=16)
plt.tight_layout()
plt.savefig(f'{sys.argv[7]}_all_3.png',dpi=600,transparent=True)
plt.close()

### Out to 2nd digit

# prep list of relevant ecids possible in the pdb70 structure library
ecids_list  = [ecid for ecid in pdb70_all_ecids if ecid[-1] == '*' and len(ecid.split('.')) == 3]
nECIDs = len(ecids_list)
total_pseudocount = 1
pseudocount = total_pseudocount / nECIDs

# fill counts with zeros or the observed
counts = {}
for ecid in ecids_list:
    counts[ecid] = all_counts.get(ecid,0)
# get number of relevant EC IDs counted, whether defunct or incomplete or active
nCounts = np.sum([counts[ecid] for ecid in ecids_list])

# calculate the posterior probabilities with pseudocounts
posterior = calc_prop(counts,pseudocount=pseudocount)
rel_entropy_components = [[ecid, calc_rel_entropy_component(posterior[ecid],all_prior[ecid])] for ecid in ecids_list]
rel_entropy = np.sum([elem[1] for elem in rel_entropy_components])

# gather relative entropy componenets of actually observed EC IDs
rel_entropies = [elem for elem in rel_entropy_components if counts[elem[0]] != 0]
rel_entropies = sorted(rel_entropies, key=lambda x: x[1])
    
# calculate heights for logo figure
rel_heights = [[ecid, rel_entropy*counts[ecid]/nCounts] for ecid in counts if counts[ecid] != 0]
rel_heights = sorted(rel_heights,key=lambda x: x[1])

print(f'\n{nCounts} second-digit EC IDs were observed, including defunct and active IDs.')
print(f"Total relative entropy = {rel_entropy} bits. Broken down into components:\n")
for elem in rel_entropies:
    print(elem[0],elem[1])

fig = plt.figure(dpi=600)
ax = plt.gca()
x_space = 0.01
y_space = 0.01*rel_entropy
plt.plot([0-x_space,1+x_space],[rel_entropy,rel_entropy],'b--',zorder=1)
plt.xlim((0-x_space,1+x_space))
plt.ylim((0-y_space,rel_entropy+y_space))
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_major_locator(plt.NullLocator())
plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--',zorder=1)
max_length = np.max([len(ecid[0]) for ecid in rel_heights])
offset_y = 0
offset_x = 0
for i, ecid in enumerate(rel_heights):
    text_draw_mpl(fig, ax, f' {ecid[0]:<{max_length}} ', offset_x, offset_y, 1, ecid[1], color = colors[ecid[0].split('.')[0]], remove_gaps = True, zorder= i+1) #bg = 'white', 
    offset_y += ecid[1]

plt.ylabel('Relative Entropy\n(bits)',size=16)
plt.tight_layout()
plt.savefig(f'{sys.argv[7]}_all_2.png',dpi=600,transparent=True)
plt.close()

### Out to 1st digit

# prep list of relevant ecids possible in the pdb70 structure library
ecids_list  = [ecid for ecid in pdb70_all_ecids if ecid[-1] == '*' and len(ecid.split('.')) == 2]
nECIDs = len(ecids_list)
total_pseudocount = 1
pseudocount = total_pseudocount / nECIDs

# fill counts with zeros or the observed
counts = {}
for ecid in ecids_list:
    counts[ecid] = all_counts.get(ecid,0)
# get number of relevant EC IDs counted, whether defunct or incomplete or active
nCounts = np.sum([counts[ecid] for ecid in ecids_list])

# calculate the posterior probabilities with pseudocounts
posterior = calc_prop(counts,pseudocount=pseudocount)
rel_entropy_components = [[ecid, calc_rel_entropy_component(posterior[ecid],all_prior[ecid])] for ecid in ecids_list]
rel_entropy = np.sum([elem[1] for elem in rel_entropy_components])

# gather relative entropy componenets of actually observed EC IDs
rel_entropies = [elem for elem in rel_entropy_components if counts[elem[0]] != 0]
rel_entropies = sorted(rel_entropies, key=lambda x: x[1])
    
# calculate heights for logo figure
rel_heights = [[ecid, rel_entropy*counts[ecid]/nCounts] for ecid in counts if counts[ecid] != 0]
rel_heights = sorted(rel_heights,key=lambda x: x[1])

print(f'\n{nCounts} second-digit EC IDs were observed, including defunct and active IDs.')
print(f"Total relative entropy = {rel_entropy} bits. Broken down into components:\n")
for elem in rel_entropies:
    print(elem[0],elem[1])

fig = plt.figure(dpi=600)
ax = plt.gca()
x_space = 0.01
y_space = 0.01*rel_entropy
plt.plot([0-x_space,1+x_space],[rel_entropy,rel_entropy],'b--',zorder=1)
plt.xlim((0-x_space,1+x_space))
plt.ylim((0-y_space,rel_entropy+y_space))
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_major_locator(plt.NullLocator())
plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--',zorder=1)
max_length = np.max([len(ecid[0]) for ecid in rel_heights])
offset_y = 0
offset_x = 0
for i, ecid in enumerate(rel_heights):
    text_draw_mpl(fig, ax, f' {ecid[0]:<{max_length}} ', offset_x, offset_y, 1, ecid[1], color = colors[ecid[0].split('.')[0]], remove_gaps = True, zorder= i+1) #bg = 'white', 
    offset_y += ecid[1]

plt.ylabel('Relative Entropy\n(bits)',size=16)
plt.tight_layout()
plt.savefig(f'{sys.argv[7]}_all_1.png',dpi=600,transparent=True)
plt.close()


###################
# WORDCLOUD - FULL ECIDS
###################

#wordcloud_dict = {}
#for ecid in prop:
#    if len(ecid.split('.')) == 4 and ecid.split('.')[-1] != '*':
#        wordcloud_dict[ecid] = counts[ecid]
#print(wordcloud_dict)

print(f"Making the ECID word cloud.")

all_counts_keys = list(all_counts.keys())
all_counts_keys = sorted(all_counts_keys, key = lambda x: len(x))
print(all_counts_keys)
temp = {}
for key in all_counts_keys[::-1]:
    temp[key] = all_counts[key]

wc = wordcloud.WordCloud(background_color=None,
                         mode='RGBA',
                         width=400,
                         height=400,
                         include_numbers=True,
                         prefer_horizontal=1.0,
                         relative_scaling=1.0,
                         max_font_size=400)    # ,max_words=100
#wc.generate_from_frequencies(wordcloud_dict)
wc.generate_from_frequencies(temp)
#wc.generate_from_frequencies(all_counts)
wc.recolor(color_func=get_ecid_color)

plt.imshow(wc,interpolation='bilinear')
plt.axis('off')
plt.savefig(f'{sys.argv[7]}_wordcloud.png',dpi=600,transparent=True)
plt.close()


