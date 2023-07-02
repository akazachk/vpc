# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib

matplotlib.use('MACOSX')

scale=2
DPI = 300
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.spines', **{'bottom':True, 'left':True, 'right':False, 'top':False})
plt.rc('axes', titlesize=12*scale)
plt.rc('axes', labelsize=8*scale)
plt.rc('xtick', labelsize=8*scale)
plt.rc('ytick', labelsize=8*scale)
plt.rc("legend", fontsize=8*scale)
plt.rc("figure", figsize=[6*scale,4*scale])
# plt.rc("figure", figsize=[6,4])
#plt.rc("figure", figsize=[3,2])
plt.rc("savefig", dpi=DPI)

# %%
# Read in the csv files
#path = r'.' # use your path

# Print current directory
print("Current directory: ", os.getcwd())

# Read in value of VPC_DIR from environment variable; if empty, use parent directory
if 'VPC_DIR' not in os.environ:
    VPC_DIR = '..'
    print("VPC_DIR not in environment variable; using parent directory: ", os.path.abspath(VPC_DIR))
else:
    VPC_DIR = os.environ['VPC_DIR']
    print("VPC_DIR set to: ", os.path.abspath(VPC_DIR))
# path = VPC_DIR/results/2023-07-01
path = VPC_DIR + '/results/2023-07-01'

# Convert path to absolute path
path = os.path.abspath(path)

print("Path set to: ", path)

# %%
all_files = []
for file in os.listdir(path):
    if file.endswith(".csv"):
        # print(os.path.join(path, file))
        all_files.append(os.path.join(path, file))

# Sort all_files lexicographically
all_files.sort()
print(all_files)

# li = []
PREFIX = "cutrounds_"
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

def adjust_lightness(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> adjust_lightness('g', 0.3)
    >> adjust_lightness('#F034A3', 0.6)
    >> adjust_lightness((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

MARKERS = ['x', 'o', None, '*', None, None, None, None, None]
LINESTYLE_LIST = ['solid', 'dashed', 'dotted', 'dashdot']
ALPHA_LIST = [1.0, 0.75, 0.5, 0.25]
prev_instance = ""
fig = None
ax1 = None
ax2 = None
ax3 = None
ax4 = None
axes = [ax1]
color_index = 0
depth_index = 0
for filename in all_files:
    # The instance name is everything after "path/PREFIX" and before "_dXXX.csv" in the filename
    instance = filename

    # Check that PREFIX is in the filename
    if PREFIX not in instance:
        continue



    # Delete everything before and including "PREFIX"
    instance = instance[instance.find(PREFIX) + len(PREFIX):]
    
    # Delete everything after ".csv"
    instance = instance[:instance.find(".csv")]

    # Pull everything after the last underscore
    depth_stub = instance[instance.rfind("_d"):] 
    depth = int(depth_stub[2:])
    instance = instance[:instance.rfind("_d")]
    
    print("Instance: ", path + '/' + instance)
    print("Depth: {:d}".format(depth))

    # Read in the csv file; skip the last column (it is created due to comma at end of every row)
    df = pd.read_csv(filename, engine='python')
    
    # Delete the last column
    del df['Unnamed: 10']

    # Rename first column to 'Round'
    df.rename(columns={df.columns[0]: 'Round'}, inplace=True)
    
    # If we are in a Jupyter notebook, display the dataframe
    import sys
    if 'IPython' in sys.modules:
        from IPython.display import display
        display(df.head())

    # Plot for this instance the LPTime across the Round columns and print the LP bound on a second axis
    # On a third axis, which goes below the plot into the negative numbers, put a bar plot of the number of cuts added using column 'Cuts'
    # fig, ax1 = plt.subplots()
    if (prev_instance != instance):
        # Save old plot if prev_instance is not empty
        if prev_instance != "":
            # Add the legend
            plt.legend(legend_handles, legend_labels, loc='lower right', bbox_to_anchor=(1,-0.13), ncol=4, fontsize=5*scale)
            
            # Add a title
            plt.title(r"\ttfamily{"+prev_instance+r"}", fontsize=8*scale)

            fig.tight_layout()  # otherwise the right y-label is slightly clipped

            # Save the plot to a file with the instance name
            # To prevent it from being blank, save it as a pdf
            plt.savefig(path + '/' + prev_instance + ".pdf", bbox_inches='tight', pad_inches=0, dpi=300)
            
        prev_instance = instance
        fig = plt.figure(figsize=(10, 6))
        axes[0] = fig.add_subplot(111)
        depth_index = 0
    else:
        depth_index += 1

    # Plot the gmic_bound on the first axis
    color_index = 0
    curr_ax = axes[0]
    color = CB_color_cycle[color_index]
    marker = MARKERS[color_index]
    linestyle = LINESTYLE_LIST[depth_index]
    alphastyle = ALPHA_LIST[depth_index]
    
    curr_ax.set_xlabel('Rounds of cuts', horizontalalignment='left', x=0.0)
    curr_ax.set_ylabel('Bound', color='black')
    x_ticks = df['Round'].to_numpy()+1
    
    curr_ax.plot(x_ticks, df['bound_gmic'], color=color, label='gmic'+depth_stub, marker=marker, linestyle=linestyle, alpha=alphastyle)
    
    # Plot vpc bound
    if depth > 0:
        color_index += 1
        color = CB_color_cycle[color_index]
        marker = MARKERS[color_index]
        curr_ax.plot(x_ticks, df['bound_vpc'], color=color, label='vpc'+depth_stub, marker=marker, linestyle=linestyle, alpha=alphastyle, markersize=4)

    # curr_ax.tick_params(axis='y', labelcolor=color)
    # curr_ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))


    continue


    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    curr_ax = ax2
    color_index += 1

    color = CB_color_cycle[color_index]
    marker = MARKERS[color_index]
    curr_ax.set_ylabel('LP time (s)', color=color)
    curr_ax.plot(x_ticks, df['LPTime'], color=color, label='LP time (s)', marker=marker, linestyle='dashed')
    curr_ax.tick_params(axis='y', labelcolor=color)

    # Create the third y-axis and plot the number of cuts as a bar plot
    # Do not show the y-axis ticks or label
    # Make this be cumulative number of cuts
    ax3 = ax1.twinx()
    curr_ax = ax3
    color_index += 1
    color = CB_color_cycle[color_index]

    df['CumulativeCuts'] = df['Cuts'].cumsum()
    curr_ax.bar(x_ticks, df['CumulativeCuts'], color=color, label='Cumulative cuts (\#)', alpha=0.25)
    # curr_ax.set_ylabel("Number of Cuts", color=color)
    curr_ax.tick_params(axis='y', labelcolor=color)

    # Overlay a label with the number of cuts on the last bar
    # Get the last bar
    rects = curr_ax.patches
    # Get the height of the last bar
    height = rects[-1].get_height()
    # Add the text label
    curr_ax.text(rects[-1].get_x() + rects[-1].get_width() / 2, height + 5, str(int(height)), ha='center', va='bottom', color=color)

    # Hide the ticks on the third axis
    curr_ax.set_yticks([])
    curr_ax.set_yticklabels([])

    # Add fourth axis same as the first one in which the approximation for number of cuts (harmonic series) is used relative to the first round's bound improvement
    ax4 = ax1
    curr_ax = ax4
    color_index += 1
    color = adjust_lightness(CB_color_cycle[0], 1.5)
    marker = MARKERS[color_index]
    
    first_bound = df['LPBound'][0]
    first_improvement = df['LPBound'][1] - first_bound
    improvement = [first_bound]
    for i in range(1,len(df['Round'])):
        prev_bound = improvement[-1]
        improvement += [prev_bound + first_improvement / (i)]
    # print("First bound = ", first_bound)
    # print("First improvement = ", first_improvement)
    # print("First bound + first improvement = ", first_bound + first_improvement)
    # print("Improvement[0] ", improvement[0])
    # print("Improvement[1] ", improvement[1])
    # print("Second bound = ", df['LPBound'][1])
    curr_ax.plot(x_ticks, improvement, color=color, label='Predicted bound', marker=marker, linestyle='dotted')
    # curr_ax.tick_params(axis='y', labelcolor=color)

    # Prepare a single legend box
    # Concatenate the legend handles and labels
    legend_handles = ax1.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0] + ax3.get_legend_handles_labels()[0]
    legend_labels = ax1.get_legend_handles_labels()[1] + ax2.get_legend_handles_labels()[1] + ax3.get_legend_handles_labels()[1]

    # Add the legend
    plt.legend(legend_handles, legend_labels, loc='lower right', bbox_to_anchor=(1,-0.13), ncol=4, fontsize=5*scale)

    # Add a title
    plt.title(r"\ttfamily{"+instance+r"}", fontsize=8*scale)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    # Save the plot to a file with the instance name
    # To prevent it from being blank, save it as a pdf
    plt.savefig(path + '/' + instance + ".pdf", bbox_inches='tight', pad_inches=0, dpi=300)

# %%
## Plot the last instance

# Prepare a single legend box
# Concatenate the legend handles and labels
legend_handles = []
legend_labels = []
for ax in axes:
    legend_handles += ax.get_legend_handles_labels()[0]
    legend_labels += ax.get_legend_handles_labels()[1]

# Add the legend
# plt.legend(legend_handles, legend_labels, loc='lower right', bbox_to_anchor=(1,-0.13), ncol=4, fontsize=5*scale)
plt.legend(legend_handles, legend_labels, loc='lower right', bbox_to_anchor=(1.,-0.13), ncol=5, fontsize=5*scale)

# Add a title
plt.title(r"\ttfamily{"+prev_instance+r"}", fontsize=8*scale)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

# Save the plot to a file with the instance name
# To prevent it from being blank, save it as a pdf
plt.savefig(path + '/' + prev_instance + ".pdf", bbox_inches='tight', pad_inches=0, dpi=300)
# %%
