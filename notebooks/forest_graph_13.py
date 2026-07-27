# This will be a forest plot

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np



# This is just data pulled from earlier runs:
# https://pypi.org/project/forestplot/

findings = [
    # (Label, HR, CI_low, CI_high, p-value, cohort), this is all data for the forest plot to summarize my findings
    ("Grade 3 (TCGA, PFI)", 2.14, 1.06, 4.36, 0.035, "TCGA"),
    ("Grade 3 (CGGA, OS)", 2.48, 1.38, 4.46, 0.002, "CGGA validation"),
    ("Age at diagnosis (TCGA)", 1.03, 1.00, 1.07, 0.050, "TCGA"),
    ("Age at diagnosis (CGGA)", 1.03, 1.00, 1.06, 0.052, "CGGA validation"),
    ("Secondary oncogene mutation\n(NOTCH1/PIK3CA, TCGA)", 2.55, 1.27, 5.13, 0.008, "TCGA"),
    ("Transcriptomic risk score\n(9-gene, TCGA)", 1.81, 1.44, 2.27, 0.0001, "TCGA"),
    ("Transcriptomic risk score\n(9-gene, CGGA validation, log-score p)", 1.22, 0.73, 2.05, 0.005, "CGGA validation"),
    ("PAX5 promoter methylation\n(TCGA, multivariable)", 1.04, 1.00, 1.08, 0.046, "TCGA"),
] # CGGA doesnt have methylation validation

fig, ax = plt.subplots(figsize=(11, 8))

y_positions = np.arange(len(findings))[::-1]  # reverse so first item is on top, when we plot we want it to be practically chronogically correct


for i, (label, hr, ci_low, ci_high, p, source) in enumerate(findings):
    y = y_positions[i] # just pairing the index to its y-slot (this is why we reversed)
    color = '#d7191c' if source == "TCGA" else '#2c7bb6'
    marker = 'o' # just circle all
    
    ax.plot([ci_low, ci_high], [y, y], color=color, linewidth=2, alpha=0.8)
    ax.scatter(hr, y, color=color, s=120, marker=marker, zorder=3, # I had trouble with formatting the dots to put them above. Gemini helped with adding zorder = 3, which brings it to the front (z-leve, basically)
               edgecolors='black', linewidth=0.8)
    
    sig_marker = "*" if p < 0.06 else "" # just marking the 'more' significant ones
    ax.text(ci_high + 0.15, y, f'HR={hr:.2f} [{ci_low:.2f}-{ci_high:.2f}], p={p:.3f}{sig_marker}',
            va='center', fontsize=9) # This is to write the text, +0.15 is the standard safety offset I used 

ax.axvline(x=1, color='black', linestyle='--', linewidth=1) # null hypothesis (no change expected from anything)
ax.set_yticks(y_positions)
ax.set_yticklabels([f[0] for f in findings], fontsize=10)
ax.set_xlabel('Hazard Ratio (log scale)', fontsize=12) # Log scaled as below
ax.set_xscale('log') # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xscale.html necessary, just logging scale
ax.set_xlim(0.5, 12)
ax.set_title('Summary of Validated Findings\nOligodendroglioma Multi-Cohort Survival Analysis', 
             fontsize=13, fontweight='bold', pad=15)



# https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html, 
# I had trouble with labels and legends so I had to use Patch and Line2D for custom labels.

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Patch(facecolor='#d7191c', label='TCGA findings cohort'),
    Patch(facecolor='#2c7bb6', label='CGGA validation cohort'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markeredgecolor='black', markersize=10, label='p < 0.05 (more significant)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.9) # frame alpha is to make it transparent

ax.grid(True, axis='x', linestyle=':', alpha=0.4)
plt.tight_layout()
plt.savefig('figures/forest_plot_summary_13.png', dpi=300, bbox_inches='tight')
