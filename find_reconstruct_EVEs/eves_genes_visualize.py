import pandas as pd
import matplotlib.pyplot as plt

def visualize_gene_presence_transposed(input_file_path, column_colors, x_spacing, y_spacing, font_size, italic):
    # Load the data from the provided Excel file
    df = pd.read_excel(input_file_path)
    
    # Replace 'yes' with True and 'no' with False
    df = df.replace({'yes': True, 'no': False})

    # Reverse DataFrame rows to maintain the order in the plot as in the Excel file
    df = df.iloc[::-1].reset_index(drop=True)
    
    # Get organism names and gene names for labeling
    organisms = df.iloc[:, 0].tolist()
    genes = df.columns[1:].tolist()

    

    # Calculate figure size dynamically based on the number of items and spacing, plus static margins
    fig_width = len(genes) * x_spacing + 2  # Adding static margin space to width
    fig_height = len(organisms) * y_spacing + 2  # Adding static margin space to height
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Plotting each cell as a circle with adjustable marker size
    marker_size = 20  # Adjusted marker size for better visibility
    for i, gene in enumerate(genes):
        for j, organism in enumerate(organisms):
            color = column_colors[i] if df.iloc[j, df.columns.get_loc(gene)] else 'white'
            ax.plot(i * x_spacing + 1, j * y_spacing + 1, 'o', color=color, markeredgecolor='black', markersize=marker_size)

    # Set axis labels with genes and organisms, move the gene labels to the top
    ax.set_xticks([i * x_spacing + 1 for i in range(len(genes))])
    ax.set_xticklabels(genes, rotation=45, fontsize=font_size, style='italic' if italic else 'normal')
    ax.xaxis.set_label_position('top') 
    ax.xaxis.tick_top()
    ax.set_yticks([j * y_spacing + 1 for j in range(len(organisms))])
    ax.set_yticklabels(organisms, fontsize=font_size, style='italic' if italic else 'normal')

    # Hide the spines and ticks for a cleaner look
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=0)

    # Set the aspect of the plot to equal, so circles are shown as circles
    ax.set_aspect('equal')

    plt.tight_layout()


    plt.show()

    # Save the plot as a file
    output_filepath = input_file_path.replace('.xlsx', '_transposed_circle_plot.svg')
    fig.savefig(output_filepath, bbox_inches='tight', dpi=300)

    return output_filepath

# Define the colors for each gene column
column_colors = ['red', 'darkgreen', 'lightgreen', 'orange', 'magenta', 'cyan']

# Set the spacing for the circles along the x and y axes
x_spacing = 0.5
y_spacing = 0.7

# Parameters for font size and italic
font_size = 15  
italic = False  

output_file_transposed = visualize_gene_presence_transposed('/home/majnusova/Desktop/eves.xlsx', column_colors, x_spacing, y_spacing, font_size, italic)
output_file_transposed

