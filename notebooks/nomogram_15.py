import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = pd.read_csv('data/TCGA/oligo_final_combined_model.csv')

# Only making it two variable as age/grade/methylation does not contribute to much, espcially as the nomogram is supposed to be a quick tool used by clinicians
covariates = ['has_secondary_og', 'rnaseq_risk_score']
covariate_labels = ['Secondary Oncogene\n(NOTCH1/PIK3CA)',
                     '9-Gene Transcriptomic\nRisk Score']

model_cols = covariates + ['PFI.time', 'PFI']
data = data.dropna(subset=model_cols).copy()
cph = CoxPHFitter(penalizer=0.1)
cph.fit(data[model_cols], duration_col='PFI.time', event_col='PFI')
coefs = cph.params_

var_ranges = {}
for cov in covariates:
    if data[cov].nunique() <= 2: # This is for binary- the secondary oncogene mutation is only 0 or 1 
        lo, hi = 0.0, 1.0
    else:
        lo, hi = data[cov].quantile(0.05), data[cov].quantile(0.95) # Because here I trim off the extrmeely far edges which is common for the outliers/nomogram procedures
    var_ranges[cov] = (lo, hi)

raw_spans = {}
risk_direction = {}
for cov, (lo, hi) in var_ranges.items():
    coef = coefs[cov]
    raw_spans[cov] = abs(coef) * (hi - lo) # This is to find the maximum shift that each variable spans to do the nomogram/ticks
    risk_direction[cov] = 1 if coef >= 0 else -1


max_span = max(raw_spans.values()) # finds the one with biggest impact
points_max = {cov: raw_spans[cov] / max_span * 100.0 for cov in covariates} # In nomograms, the most impacting variable is assigned 100 points, the others are assigned relative to the length of that one.

total_points_max = sum(points_max.values()) # So in this case, after running we see that the max points is ~115, 100 for transcriptomic, 15 for sec onc

prob_axis_points = np.linspace(0, total_points_max, 6) # Gemini snippit, helped with marks for the nomogram (like ticks)

# We will be measuring 3,5,10 years.
timepoints = [3.0, 5.0, 10.0]
prob_axis_values = {t: [] for t in timepoints} 

for target_points in prob_axis_points:
    frac = target_points / total_points_max if total_points_max > 0 else 0  # just sees percentage of points we have
    profile = {}
    for cov in covariates:
        lo, hi = var_ranges[cov]
        if risk_direction[cov] == 1:
            profile[cov] = lo + frac * (hi - lo) # if it increases risk, start at low and add by percentage then multiply by 'weight'
        else:
            profile[cov] = hi - frac * (hi - lo) # decerases, starts high and subtracts percentage multiplied by 'weight'
    profile_df = pd.DataFrame([profile])[covariates]
    
    for t in timepoints:
        surv_prob = cph.predict_survival_function(profile_df, times=[t]).iloc[0, 0] # predicting survival at that time
        prob_axis_values[t].append(surv_prob)


# onto actual figure generation
n_rows = len(covariates) + 2 + len(timepoints)  # the 2 is for the point metric system in nomograms, one for group points one for total
fig, ax = plt.subplots(figsize=(11, 1.3 * n_rows)) # Just scaling height, all this is in variables incase I want to add or detract # of timestamps/variables
y = n_rows

def draw_axis(y_pos, x_left, x_right, label, ticks_x, ticks_label, color='black'):
    ax.plot([x_left, x_right], [y_pos, y_pos], color=color, linewidth=1.8, zorder=2) # These are for invidivual lines the nomogram has, putting them above grid lines
    for tx, tl in zip(ticks_x, ticks_label):
        ax.plot([tx, tx], [y_pos - 0.08, y_pos + 0.08], color=color, linewidth=1.3, zorder=2)
        ax.text(tx, y_pos + 0.18, tl, ha='center', va='bottom', fontsize=9)
    ax.text(x_left - 8, y_pos, label, ha='right', va='center', fontsize=11, fontweight='bold')

for pts in range(0, 101, 20):
    ax.axvline(x=pts, color='gray', linestyle=':', linewidth=0.5, alpha=0.4, zorder=1)

y -= 1.0
pts_ticks = np.arange(0, 101, 10)
draw_axis(y, 0, 100.0, "Points", pts_ticks, [str(t) for t in pts_ticks], color='#333333')

for cov, label in zip(covariates, covariate_labels):
    y -= 1.0
    lo, hi = var_ranges[cov]
    is_binary = data[cov].nunique() <= 2
    if is_binary:
        val_ticks = np.array([lo, hi])
        fmt = "{:.0f}"
    else:
        val_ticks = np.linspace(lo, hi, 5)
        fmt = "{:.1f}"
    if risk_direction[cov] == 1:
        pt_ticks = np.interp(val_ticks, [lo, hi], [0, points_max[cov]])
    else:
        pt_ticks = np.interp(val_ticks, [lo, hi], [points_max[cov], 0])
    x_ticks = pt_ticks  # <-- was: pt_ticks / total_points_max * 100.0
    tick_labels = [fmt.format(v) for v in val_ticks]
    draw_axis(y, min(x_ticks), max(x_ticks), label, x_ticks, tick_labels, color='#2c7bb6')
    
y -= 1.0
total_ticks = np.linspace(0, total_points_max, 6)
x_total_ticks = total_ticks / total_points_max * 100.0 if total_points_max > 0 else total_ticks
draw_axis(y, 0, 100.0, "Total Points", x_total_ticks, [f"{t:.0f}" for t in total_ticks], color='#d7191c')

# Three separate probability axes, one per timepoint
timepoint_colors = {3.0: '#1a9641', 5.0: '#fdae61', 10.0: '#d7191c'}
for t in timepoints:
    y -= 1.0
    x_prob_ticks = prob_axis_points / total_points_max * 100.0 if total_points_max > 0 else prob_axis_points
    draw_axis(y, 0, 100.0, f"Predicted {int(t)}-yr\nPFS Probability", x_prob_ticks,
              [f"{p:.2f}" for p in prob_axis_values[t]], color=timepoint_colors[t])

ax.set_xlim(-25, 105)
ax.set_ylim(0, n_rows)
ax.axis('off')
ax.set_title('Oligodendroglioma Prognostic Nomogram\n(Secondary Oncogene + Transcriptomic Risk Score)',
             fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('figures/16_nomogram.png', dpi=300, bbox_inches='tight')
print("Saved figures/16_nomogram.png")
