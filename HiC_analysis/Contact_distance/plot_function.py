import numpy as np
import pandas as pd
import os
import cooler
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import matplotlib.patches as patches
from statsmodels.stats.multitest import multipletests
from matplotlib.colors import TwoSlopeNorm

def ratio_by_groups_scatter(
    df,
    plot_dict,               # dict: {lineage: [(level_col, label, age_val), ...]}
    age_color_map,           # dict: {age_val: color}
    age_col='age_groups',
    min_cells_per_tuple=200,
    figsize=(10, 3), dpi=150,
    gap_size=1, ylim=None,
    show_group_labels=True,
# min_cells_per_age=10000
):
    """Plot log2(SE/LE) across multiple lineages with Case-level jitter points and ANOVA per lineage."""

    # ---- collect and label data ----
    all_rows, labels, groups, x_positions = [], [], [], []
    pos = 0
    for lineage, tuples in plot_dict.items():
        for (level_col, label, age) in tuples:
            mask = (df[level_col] == label) & (df[age_col] == age)
            sub = df.loc[mask, ['SE_over_LE_ratio','Interaction_count','Case']].copy()
            if not sub.empty:
                sub['log2_ratio'] = np.log2(np.clip(sub['SE_over_LE_ratio'], 1e-9, None))
                sub['log10_inter'] = np.log10(np.clip(sub['Interaction_count'], 1e-9, None))
                sub['tuple_id'], sub['group'], sub['age'] = len(labels), lineage, age
                all_rows.append(sub)
            labels.append(f"{label} | {age}")
            groups.append(lineage)
            x_positions.append(pos)
            pos += 1
        pos += gap_size  # spacing between lineages

    dat = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    if dat.empty:
        print("[WARN] No matching data found.")
        return None, pd.DataFrame()
    pair_counts = dat.groupby(['Case', 'tuple_id']).size()
    valid_pairs = pair_counts[pair_counts >= min_cells_per_tuple].index
    dat = dat.set_index(['Case', 'tuple_id'])
    dat = dat.loc[dat.index.isin(valid_pairs)].reset_index()
    # ---- filter small tuples ----
    # keep_ids = dat['tuple_id'].value_counts()
    # valid_ids = keep_ids[keep_ids >= min_cells_per_tuple].index
    # dat = dat[dat['tuple_id'].isin(valid_ids)]
    if dat['tuple_id'].nunique() < 2:
        print("[WARN] Too few valid tuples after filtering.")
        return None, pd.DataFrame()

    # ---- plot ----
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    data = [dat.loc[dat['tuple_id'] == i, 'log2_ratio'].values for i in range(len(labels))]
    colors = [age_color_map.get(t[2], 'gray')
              for lineage, tuples in plot_dict.items() for t in tuples] + ['gray'] * gap_size

    # --- boxplots ---
    bp = ax.boxplot(
        data, positions=x_positions, widths=0.6,
        showfliers=False, patch_artist=True, zorder=1
    )

    for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
        patch.set_zorder(1)   # box below scatter

    for med in bp['medians']:
        med.set_color('white')
        med.set_zorder(1)

    # --- scatter (drawn on top) ---
    for tid, xpos in enumerate(x_positions):
        sub = dat.loc[dat['tuple_id'] == tid]
        if sub.empty:
            continue
        case_mean = sub.groupby('Case')['log2_ratio'].median().reset_index()
        jitter = np.random.normal(0, 0.1, size=len(case_mean))
        ax.scatter(
            np.full(len(case_mean), xpos) + jitter,
            case_mean['log2_ratio'],
            s=15, alpha=1, edgecolor='black', linewidth=0.3,
            color='red',
            rasterized=True,
            zorder=3    # <- higher zorder = drawn later
        )

    # ---- axis & labels ----
    # x_positions = [i + 0.5 for i in x_positions]
    labels = [i.split(' | ')[1] for i in labels]
    ax.set_xticks(x_positions)

    ax.set_xticklabels(labels, rotation=90, ha='right', fontsize=8)
    ax.set_ylabel('log2(SE/LE)')
    ax.set_title(f"")

    # ---- lineage dividers & group labels ----
    cursor, centers = 0, []
    for lineage, tuples in plot_dict.items():
        n = len(tuples)
        centers.append((cursor + (n - 1)/2, lineage))
        cursor += n + gap_size
        if cursor < max(x_positions):
            ax.axvline(cursor - gap_size/2, color='black', ls='--', lw=0.8, alpha=0.6)

    if show_group_labels:
        if ylim is None:
            ylim = ax.get_ylim()
        for xc, lineage in centers:
            ax.text(xc, ylim[1] + 0.01*(ylim[1]-ylim[0]), lineage,
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
        ax.set_ylim(ylim[0], ylim[1] * 1.15)

    plt.tight_layout()

    # ---- ANOVA per lineage ----
    results = []
    for lineage in plot_dict.keys():
        sub = dat[dat['group'] == lineage]
        if sub['tuple_id'].nunique() > 1:
            fit = smf.ols("log2_ratio ~ C(tuple_id) + log10_inter", data=sub).fit()
            aov = anova_lm(fit)
            p = aov.loc["C(tuple_id)", "PR(>F)"] if "C(tuple_id)" in aov.index else np.nan
            results.append({'lineage': lineage, 'anova_p': p, 'n_cells': len(sub)})

    res = pd.DataFrame(results)
    if not res.empty:
        _, q, _, _ = multipletests(res['anova_p'].fillna(1), method='fdr_bh')
        res['anova_q'] = q

    return fig, res    

import os
import numpy as np
import cooler
import matplotlib.pyplot as plt
from matplotlib import patches


def heatmap_by_L2_peak_balanced(
    normalized_weight, meta, 
    L1 = None,
    L2=None, L3=None,
    age_order=('2T','3T','1m','4-7m','adult'),
    cmap='viridis',
    split_by_region=True,
    min_cells_per_age=200,
    n_per_age=None,
    balance=True,
    figsize = None,
    seed=0
):
    rng = np.random.default_rng(seed)

    # --- subset meta for the chosen cluster ---
    if L1 is not None:
        meta_sub = meta
        cluster_label = f"L1={L1}"

    elif L2 is None:
        meta_sub = meta.loc[normalized_weight.index, ['L3','age_groups']].dropna()
        meta_sub = meta_sub[(meta_sub['L3'] == L3) & (meta_sub['age_groups'].isin(age_order))]
        cluster_label = f"L3={L3}"
    else:
        meta_sub = meta.loc[normalized_weight.index, ['L2','age_groups']].dropna()
        meta_sub = meta_sub[(meta_sub['L2'] == L2) & (meta_sub['age_groups'].isin(age_order))]
        cluster_label = f"L2={L2}"
    if meta_sub.empty:
        raise ValueError(f"No cells for {cluster_label}")

    # --- per-age available counts ---
    avail = {age: int((meta_sub['age_groups'] == age).sum()) for age in age_order}
    keep_ages = [a for a in age_order if avail[a] >= min_cells_per_age]
    if len(keep_ages) < 2:
        raise ValueError(f"Need >=2 ages with >={min_cells_per_age} cells; got {avail}")

    # --- decide number per age (only if balancing) ---
    if balance:
        if n_per_age is None:
            n_per_age = min(avail[a] for a in keep_ages)
        n_per_age = int(n_per_age)
        if any(avail[a] < n_per_age for a in keep_ages):
            n_per_age = min(avail[a] for a in keep_ages)
    else:
        n_per_age = None  # flag for using all cells

    # --- sample + order within each age (balance-aware) ---
    cell_order, cuts = [], []
    k = 0
    for age in keep_ages:
        idx_age = meta_sub.index[meta_sub['age_groups'] == age]
        sub_all = normalized_weight.loc[idx_age]
        peak_bin_all = sub_all.values.argmax(axis=1)
        age_meta = meta_sub.loc[idx_age].copy()
        age_meta['peak_bin'] = peak_bin_all

        chosen = []
        if balance:
            # same balancing logic as before
            if split_by_region:
                regs, counts = np.unique(age_meta['merge_regions'].values, return_counts=True)
                props = counts / counts.sum()
                alloc = np.floor(props * n_per_age).astype(int)
                remainder = n_per_age - alloc.sum()
                frac = props * n_per_age - np.floor(props * n_per_age)
                for r in np.argsort(frac)[::-1][:remainder]:
                    alloc[r] += 1

                for reg, take in zip(regs, alloc):
                    idx_reg = age_meta.index[age_meta['merge_regions'] == reg]
                    if take == 0 or len(idx_reg) == 0:
                        continue
                    pick = rng.choice(idx_reg, size=min(take, len(idx_reg)), replace=False)
                    order = np.argsort(age_meta.loc[pick, 'peak_bin'])
                    chosen.extend(list(np.array(pick)[order]))
            else:
                pick = rng.choice(age_meta.index, size=min(n_per_age, len(idx_age)), replace=False)
                order = np.argsort(age_meta.loc[pick, 'peak_bin'])
                chosen = list(np.array(pick)[order])
        else:
            # ---- no balancing: use all cells ----
            order = np.argsort(age_meta['peak_bin'])
            chosen = list(age_meta.index[order])

            # compute region cuts even when not balancing
            region_cuts = []
            if split_by_region:
                labels_seq = age_meta.loc[chosen, 'merge_regions'].tolist()
                last = labels_seq[0] if labels_seq else None
                for j, lab in enumerate(labels_seq):
                    if j > 0 and lab != last:
                        region_cuts.append((k - len(chosen) + j, None))
                        last = lab

        cell_order.extend(chosen)
        k += len(chosen)

        # compute region cuts (same as before)
        region_cuts = []
        if split_by_region:
            labels_seq = age_meta.loc[chosen, 'merge_regions'].tolist()
            last = labels_seq[0] if labels_seq else None
            for j, lab in enumerate(labels_seq):
                if j > 0 and lab != last:
                    region_cuts.append((k - len(chosen) + j, None))
                    last = lab
        cuts.append((k, age, region_cuts))

    # --- rest of the code unchanged (plotting, dashed lines, etc.) ---
    mat_df = normalized_weight.loc[cell_order]
    mat = mat_df.T.values
    # vmin, vmax = 0, np.nanpercentile(mat, 99)
    vmin, vmax = 0, 0.015
    # vmax = 0.02
    if vmax == 0: vmax = 1e-9

    fig, ax = plt.subplots(figsize=(6, 1.6) if figsize is None else  figsize , dpi=120)
    im = ax.imshow(mat, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True, interpolation='none')

    prev_age_edge = 0
    for cut_age, age, region_cuts in cuts:
        if split_by_region:
            for rc, _ in region_cuts:
                ax.axvline(rc - 0.5, color='white', ls='-', lw=0.5, alpha=0.9)
        if prev_age_edge > 0:
            ax.axvline(prev_age_edge - 0.5, color='black', ls='--', lw=1.2)
        xc = (prev_age_edge + cut_age - 1) / 2
        ax.text(xc, mat.shape[0] + 0.5, age, ha='center', va='bottom', fontsize=9)
        prev_age_edge = cut_age

    yt = [np.log2(xx/2500.0)/0.125 for xx in [5_000, 50_000, 200_000, 2_000_000, 20_000_000, 100_000_000]]
    ax.set_yticks(yt)
    ax.set_yticklabels(['5k','50k','200k','2M','20M','100M'],size = 6)
    ax.set_ylabel('')
    ax.set_xticks([])
    # ax.set_xlabel('Cells (balanced per age{})'.format(' and region' if split_by_region else ''))
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label='Normalized weight')
    plt.tight_layout()
    return fig, ax