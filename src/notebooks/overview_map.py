import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap

# Domain parameters
L = 100.0  # Domain size
N = 101  # Number of cells
h1 = 25.0  # Left head boundary (example value)
h2 = 24.5  # Right head boundary (example value)
q = -500  # Pumping rate (example value)
con_max = 100.0  # Maximum concentration (example value)

# Create grid
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)

# Initialize background grid
domain = np.zeros((N, N))

# Mark constant head boundaries (left and right)
domain[:, 0] = 25  # Left boundary
domain[:, -1] = 24.5  # Right boundary

# Well locations
wells = [    (41, 26),
    (34, 26),
    (26, 26)]


# Contamination sources
contamination_sources =[
    (41, 26),
      (37, 26),
       (43, 26),
        (34, 26),
         (43, 24),
          (34, 23)
           ]

# Observation well
obs_well = (26, 56)  # (row, col)

# Create the plot
plt.figure(figsize=(12, 10))

# Create custom colormap
cmap = ListedColormap(['white', 'lightblue', 'red', 'green', 'purple', 'black'])
boundaries =[0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
norm = plt.cm.colors.BoundaryNorm(boundaries, cmap.N)

# Plot domain background
plt.imshow(domain, cmap=cmap, norm=norm,
extent=[0, L, 0, L], origin='lower', alpha=0.3)

# Add features
#plt.text(0, L / 2, f'CHD: h={h1}', fontsize = 10, ha = 'right', va = 'center', rotation = 90)
#plt.text(L, L / 2, f'CHD: h={h2}', fontsize=10, ha='left', va='center', rotation=-90)

# Mark wells
for well in wells:
    plt.scatter(x[well[1]], y[well[0]], s=100, c='red', marker='o', edgecolor='k', label='Well')

# Mark contamination sources
for source in contamination_sources:
    plt.scatter(x[source[1]], y[source[0]], s=100, c='green', marker='^', edgecolor='k', label='Contamination')

# Mark observation well
plt.scatter(x[obs_well[1]], y[obs_well[0]], s=150, c='purple', marker='*', edgecolor='k', label='Observation')

# Add legend with unique entries
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='upper right')

# Add grid and labels
plt.grid(True, alpha=0.3)
plt.title('Groundwater Model Domain', fontsize=16)
plt.xlabel('X Distance (m)', fontsize=12)
plt.ylabel('Y Distance (m)', fontsize=12)
plt.xticks(np.arange(0, L + 1, 10))
plt.yticks(np.arange(0, L + 1, 10))

# Add scale bar
plt.plot([L - 20, L - 10], [5, 5], 'k-', lw=2)
plt.text(L - 15, 6, '10 m', ha='center', fontsize=10)

plt.tight_layout()
plt.savefig('groundwater_model_domain.png', dpi=300)
plt.show()