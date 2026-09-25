## Env: allcools
import xarray as xr
import pandas as pd
import numpy as np
import os 
import pyranges as pr
# import umap
from scipy import ndimage as nd
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
# import scanpy as sc
import seaborn as sns
import anndata as ad
from sklearn.preprocessing import normalize
from scipy.signal import find_peaks
from scipy.stats import norm
from joblib import Parallel, delayed
from pathlib import Path
from collections import Counter
import glob
import math
import cooler
import re
import os
from joblib import Parallel, delayed
from scipy import sparse
import h5py
import anndata
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42 
matplotlib.rcParams['font.family'] = 'Arial'
import schicluster
PACKAGE_DIR = schicluster.__path__[0]


from statsmodels.sandbox.stats.multicomp import multipletests as FDR
from scipy.stats import chi2_contingency

chrom_size_path = '/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hg38.chrom_1-22.sizes'
chrom_sizes = cooler.read_chromsizes(chrom_size_path, all_names=True)

def diff_bound(bound_count_ct, cell_count_ct):
    """
    bound_count_ct: array-like, shape (n_groups, n_bins)
    cell_count_ct : array-like, shape (n_groups,)
    returns:
      chi2_stats: (n_bins,) chi-square stat per bin
      fdr:        (n_bins,) BH-FDR across bins
    """
    import numpy as np
    from scipy.stats import chi2_contingency
    from statsmodels.stats.multitest import multipletests

    # Force plain ndarrays (avoid numpy.matrix / pandas objects)
    B = np.asarray(bound_count_ct, dtype=np.int64)
    N = np.asarray(cell_count_ct, dtype=np.int64).ravel()  # (n_groups,)

    # Sanity check
    if B.shape[0] != N.size:
        raise ValueError(f"Groups mismatch: B has {B.shape[0]} rows, N has {N.size} elements")

    # Non-boundary counts for each group/bin
    T = N[:, None] - B  # (n_groups, n_bins)

    # Build a 1-D boolean mask of bins that have at least one boundary and one non-boundary across groups
    b_any = (B.sum(axis=0) > 0).ravel()
    t_any = (T.sum(axis=0) > 0).ravel()
    binfilter = (b_any & t_any)  # (n_bins,)

    n_bins = B.shape[1]
    stats = np.full(n_bins, np.nan, dtype=np.float64)
    pvals = np.ones(n_bins, dtype=np.float64)

    # Per-bin chi-square
    for i in range(n_bins):
        if not binfilter[i]:
            continue
        cont = np.vstack([B[:, i], T[:, i]])  # 2 x n_groups
        # If all groups identical, chi2 can fail; guard with try/except
        try:
            chi2, p, _, _ = chi2_contingency(cont, correction=False)
        except Exception:
            chi2, p = np.nan, 1.0
        stats[i] = chi2
        pvals[i] = p

    # FDR across bins
    try:
        _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    except Exception:
        fdr = np.full_like(pvals, np.nan, dtype=np.float64)

    return stats, fdr


def diff_bound_bulk(ins_count):
    stats = np.zeros(ins_count.shape[2])
    pv = np.ones(ins_count.shape[2])
    binfilter = (ins_count.min(axis=(0,1))>0)
    for i in range(ins_count.shape[2]):
        if binfilter[i]:
            stats[i], pv[i], _, _ = chi2_contingency(ins_count[:,:,i])
    fdr = FDR(pv, 0.001, 'fdr_bh')[1]
    return stats, pv

def domain_anova(
    leg,
    ins_count,
    bound_count_ct,
    bound_prob_ct,
    cell_count_ct,
    binall,
    chrom_sizes=chrom_sizes,
    outdir=None,
    fname = None,
    save=True
):
    if os.path.exists(f'{outdir}{fname}_bin_stats.hdf'):
        binall_computed = pd.read_hdf(f'{outdir}{fname}_bin_stats.hdf', key='data')

        return binall_computed

    binall = binall.copy()
    os.makedirs(outdir, exist_ok=True)

    # --- 1. Insulation matrix subset ---
    ins = ins_count['ratio'].to_pandas().reindex(leg).dropna(how="all")

    # --- 2. Chi-square ANOVA test ---
    chi2sc, fdr_sc = diff_bound(
        bound_count_ct.loc[leg].values,
        cell_count_ct.loc[leg].values
    )

    # Standardize chi-square values (Z-score)
    chi2_nonzero = chi2sc[chi2sc > 0]
    ave = np.mean(chi2_nonzero) if len(chi2_nonzero) > 0 else 0
    stdev = np.std(chi2_nonzero) if len(chi2_nonzero) > 0 else 1
    binall['chi2filter'] = ((chi2sc - ave) / stdev) > norm.isf(0.025)

    # --- 3. Detect insulation local minima per group ---
    binall['ins_lm'] = 0
    for xx in leg:
        sel_all = []
        for c in chrom_sizes.index:
            idx = np.where(binall['chrom'] == c)[0]
            if len(idx) > 0:
                data = -ins.loc[xx, idx]
                peaks, _ = find_peaks(data, distance=5)
                if len(peaks) > 0:
                    sel_all.append(idx.min() + peaks)
        if len(sel_all) > 0:
            sel_all = np.concatenate(sel_all)
            binall.loc[binall.index[sel_all], 'ins_lm'] = 1

    # --- 4. Boundary probability difference ---
    # Ensure bound_prob_ct matches group order
    bp = bound_prob_ct.loc[leg]
    binall['probdiff'] = (bp.max(axis=0) - bp.min(axis=0)).values
    if len(leg) == 2:
        binall['probdiff_signed'] = np.sign(bp.iloc[0] - bp.iloc[1])

    # --- 5. Store chi-square + insulation fold change ---
    binall['chi2_sc'] = chi2sc
    binall['insfc'] = ((ins.max(axis=0) + 0.01) / (ins.min(axis=0) + 0.01)).values

    # --- 6. Detect significantly differential domain bins ---
    sig_mask = fdr_sc < 1e-2
    thres = np.min(chi2sc[sig_mask])

    sel_bins = []
    for c in chrom_sizes.index:
        idx = np.where(binall['chrom'] == c)[0]
        if len(idx) > 0:
            data = chi2sc[idx]
            peaks, _ = find_peaks(data, height=thres, distance=5)
            if len(peaks) > 0:
                sel_bins.append(idx.min() + peaks)
    sel_bins = np.concatenate(sel_bins) if len(sel_bins) > 0 else np.array([], dtype=int)

    binall['diff_sc'] = 0
    if len(sel_bins) > 0:
        binall.loc[binall.index[sel_bins], 'diff_sc'] = 1

    # --- 7. Clean categorical columns ---
    for col in binall.select_dtypes(['category']).columns:
        binall[col] = binall[col].astype(str)

    # --- 8. Save results ---
    if save:
        binall.to_hdf(f'{outdir}{fname}_bin_stats.hdf', key='data')
    return binall

def plot_domain_window(
    bin_id,
    leg,
    leg_binall,
    binall,
    bound_prob_ct,
    ins_count,
    cooler_root="/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hic_pseudobulk_imputed/age_L3_25kb",
    resl=25_000,
    lslop=3_000_000,
    rslop=3_000_000,
    vmax=0.01,
    cmap="afmhot_r",
    marker_track="ins_lm",   # could be 'ins_lm', 'diff_sc', etc. (must exist in binall)
    dpi=110,
    figsize_per_row=0.5,
    annotation=None,
):
    """
    Plot a domain-centered panel:
      Top: chi-square per-bin (+ difference heatmap if len(leg)==2)
      Then, for each group in `leg`: rotated 25kb contact map + boundary probability track

    Parameters
    ----------
    bin_id : str
        Row key in binall, e.g. 'chr1_1000' that defines the focal domain.
    leg : list[str]
        List of group names (must match rows of bound_prob_ct and cooler folders).
    binall : pd.DataFrame
        Bin-level table with columns ['chrom','start','end','chi2_sc', marker_track, ...].
    bound_prob_ct : pd.DataFrame
        Boundary probability per group (rows) x bin (columns). Must include groups in `leg`.
    cooler_root : str
        Root folder where each group's {group}.cool lives under {cooler_root}/{group}/.
    resl : int
        Bin size (default 25kb).
    lslop, rslop : int
        Left/right flanks (bp) around the focal domain to plot.
    vmax, cmap : float, str
        Heatmap display settings.
    marker_track : str
        Name of a boolean/int column in `binall` to scatter as markers (e.g. 'ins_lm' or 'diff_sc').
    dpi, figsize_per_row : int, float
        Figure rendering controls.
    annotation : dict, optional
        Custom tick labels as {label: bp_position}.
    """

    # --- Resolve the focal window ---
    info = binall.loc[bin_id]
    chrom = info['chrom']
    diff_domain_l, diff_domain_r = int(info['start']), int(info['end'])

    ll, rr = (diff_domain_r - lslop), (diff_domain_l + rslop)
    loopl = max(0, ll // resl)
    loopr = max(loopl + 1, rr // resl)  # ensure at least one bin

    # --- Fetch contact maps, rotate 45° ---
    dstall = []
    for group in leg:
        base = os.path.join(cooler_root, group)
        if os.path.exists(base + ".cool"):
            cool_path = base + ".cool"
        elif os.path.exists(base + ".Q.cool"):
            cool_path = base + ".Q.cool"
        else:
            raise FileNotFoundError(f"No .cool or .Q.cool found for {group}")
        if not os.path.exists(cool_path):
            raise FileNotFoundError(f"Cool file not found: {cool_path}")
        c = cooler.Cooler(cool_path)
        # Ensure chrom exists in this cooler
        if chrom not in c.chromnames:
            raise ValueError(f"Chrom {chrom} not found in {cool_path}")
        Q = c.matrix(balance=False, sparse=True).fetch(chrom).tocsr()
        sub = Q[loopl:loopr, loopl:loopr].toarray()
        dst = nd.rotate(sub, 45, order=0, reshape=True, prefilter=False, cval=0)
        dstall.append(dst)

    # ---- Figure layout: per group we add [map, prob] ----
    # If pairwise comparison (len(leg)==2), add difference heatmap panel
    n_groups = len(leg)
    height_pattern = [3.5, 0.7]  # map, prob
    
    # if len(leg) == 2:
    #     # Add difference heatmap panel after chi2
    #     height_ratios = [1.2, 3.5] + np.tile(height_pattern, n_groups).tolist()
    #     n_rows = 2 + 2 * n_groups
    # else:
    height_ratios = [1.2] + np.tile(height_pattern, n_groups).tolist()
    n_rows = 1 + 2 * n_groups

    fig_h = max(4, np.sum(height_ratios) * figsize_per_row)
    fig, axes = plt.subplots(
        n_rows, 1,
        figsize=(6, fig_h),
        gridspec_kw={'height_ratios': height_ratios, 'hspace':0.1},
        dpi=dpi,
        sharex='all'
    )
    # --- X axis in "diamond" coordinates ---
    sel = (binall['chrom'].astype(str) == str(chrom)) & (binall['start'] >= ll) & (binall['start'] < rr)
    xpos = ((binall.loc[sel, 'start'] // resl) - loopl) * np.sqrt(2)

    # --- Top panel: chi-square track (if present) ---
    row = 0
    ax = axes[row]; row += 1
    if 'chi2_sc' in binall.columns:
        ax.plot(xpos, binall.loc[sel, 'chi2_sc'].values, c='C0', alpha=0.8)
        ax.set_title(f"{chrom}:{diff_domain_l}-{diff_domain_r}  |  Chi-square", fontsize=10)
    else:
        ax.set_title(f"{chrom}:{diff_domain_l}-{diff_domain_r}", fontsize=10)
    ax.set_xlim([0, max(1, (loopr - loopl - 1) * np.sqrt(2))])

    # --- Difference heatmap panel (only for pairwise comparisons) ---
    # if len(leg) == 2:
    #     ax_diff = axes[row]; row += 1
    #     ax_diff.set_title(f"Difference: {leg[0]} - {leg[1]}", fontsize=10)
    #     for sp in ("right", "top", "bottom", "left"):
    #         ax_diff.spines[sp].set_visible(False)
        
    #     # Compute difference
    #     diff_map = np.log2(dstall[0] / dstall[1])
        
    #     # Symmetric color scale around zero
    #     vmax_diff = np.percentile(np.abs(diff_map[diff_map != 0]), 98)
    #     img_diff = ax_diff.imshow(
    #         diff_map,
    #         cmap="RdBu",  # Red for positive (group1 > group2), Blue for negative
    #         vmin=-vmax_diff,
    #         vmax=vmax_diff
    #     )
        
    #     h = diff_map.shape[0]
    #     ax_diff.set_ylim([0.5*h, 0.35*h])
    #     ax_diff.set_xlim([0, h])
    #     ax_diff.set_yticks([])
    #     ax_diff.set_xticks([])
        
    #     # Add colorbar for difference
    #     cax_diff = fig.add_axes([0.99, 0.35, 0.02, 0.1])
    #     cbar_diff = fig.colorbar(img_diff, cax=cax_diff)
    #     cbar_diff.set_label("Contact diff.", fontsize=8)
    #     cbar_diff.ax.tick_params(labelsize=7)
        
    # --- For each group: contact map + boundary probability track ---
    ins_df = None
    if ins_count is not None:
        ins_df = ins_count['ratio'].to_pandas()
        ins_df.columns = binall.index
    
    # --- Calculate global y-axis limits for all groups ---
    # Boundary probability limits
    all_prob_vals = []
    for group in leg:
        y_prob = bound_prob_ct.loc[group, sel].values
        all_prob_vals.extend(y_prob[np.isfinite(y_prob)])
    
    if len(all_prob_vals) > 0:
        prob_ymax = min(0.3, 1.2 * np.nanmax(all_prob_vals))
    else:
        prob_ymax = 0.3
    prob_ymax = 0.2
    # Insulation limits
    all_ins_vals = []
    if ins_df is not None:
        for group in leg:
            if group in ins_df.index:
                y_ins = ins_df.loc[group, sel].values
                all_ins_vals.extend(y_ins[np.isfinite(y_ins)])
    
    if len(all_ins_vals) > 0:
        ins_lo = np.percentile(all_ins_vals, 2)
        ins_hi = np.percentile(all_ins_vals, 99)
        if ins_lo == ins_hi:
            ins_lo, ins_hi = np.min(all_ins_vals), np.max(all_ins_vals) + 1e-6
    else:
        ins_lo, ins_hi = 0.05, 0.4
    # ins_hi = 1
    
    for i, group in enumerate(leg):
        # Contact map
        ax_map = axes[row]; row += 1
        ax_map.set_title(group, fontsize=10)
        for sp in ("right", "top", "bottom", "left"):
            ax_map.spines[sp].set_visible(False)
        img = ax_map.imshow(dstall[i], cmap=cmap, vmin=0, vmax=vmax)
        h = dstall[i].shape[0]
        ax_map.set_ylim([0.5*h, 0.35*h])
        ax_map.set_xlim([0, h])
        ax_map.set_yticks([]); ax_map.set_xticks([])
        
        # Add colorbar for the first group only, positioned outside
        if i == 0:
            # create a new axis outside the main plot
            cax = fig.add_axes([0.99, 0.5, 0.02, 0.1])
            cbar = fig.colorbar(img, cax=cax)
            cbar.set_label("Contact frequency", fontsize=8)
            cbar.ax.tick_params(labelsize=7)

        # Boundary probability (+ insulation on twin y)
        ax_prob = axes[row]; row += 1
        y_prob = bound_prob_ct.loc[group, sel].values
        ax_prob.plot(xpos, y_prob, c='C0', alpha=0.9, lw=1.2, label='Boundary prob.')
        if marker_track in binall.columns:            
            marker_sel = sel & (binall[marker_track].astype(int) > 0) & (leg_binall['chi2filter']) & (leg_binall['ins_lm']) & (leg_binall['probdiff'] > 0.05)
            tmpd = (binall.loc[marker_sel, 'start'] // resl) - loopl
            ax_prob.scatter(tmpd * np.sqrt(2), np.zeros(len(tmpd)) + 0.05, color='r', s=6, alpha=0.7)
        
        # Use global y-axis limits for boundary probability
        ax_prob.set_ylim([0, prob_ymax])
        # Add y-ticks for boundary probability
        prob_ticks = np.linspace(0, prob_ymax, 2)
        ax_prob.set_yticks(prob_ticks)
        ax_prob.set_yticklabels([f"{y:.2f}" for y in prob_ticks], fontsize=8)
        ax_prob.set_ylabel("Boundary \n prob.", fontsize=9)

        # plot insulation on twin axis (optional)
        if ins_df is not None and group in ins_df.index:
            ax_ins = ax_prob.twinx()
            y_ins = ins_df.loc[group, sel].values
            ax_ins.plot(xpos, y_ins, c='C1', alpha=0.8, lw=1.0, label='Insulation')
            
            # Use global y-axis limits for insulation
            ax_ins.set_ylim(ins_lo, ins_hi)
            # Add y-ticks for insulation score
            ins_ticks = np.linspace(ins_lo, ins_hi, 2)
            ax_ins.set_yticks(ins_ticks)
            ax_ins.set_yticklabels([f"{y:.2f}" for y in ins_ticks], fontsize=8)
            ax_ins.set_ylabel("Insulation", fontsize=9)
    # --- Shared X formatting: tick labels in Mb ---
    # ax.set_xlim([0, max(1, (loopr - loopl - 1) * np.sqrt(2))])
    # xticks_base = np.arange(0, loopr - loopl + 1, 100)  # every 100 bins ~ 2.5 Mb at 25kb
    # xticks = np.sqrt(2) * xticks_base
    # plt.setp(axes, xticks=xticks)
    # # Mb labels
    # xticklabels = [f'{((xx + loopl) * resl) / 1e6:.2f}M' for xx in xticks_base]
    # plt.setp(axes, xticklabels=xticklabels)

    xmax = max(1, (loopr - loopl - 1) * np.sqrt(2))
    ax.set_xlim([0, xmax])
    # default ticks every 100 bins (~2.5 Mb for 25kb bins)
    xticks_base = np.arange(2, loopr - loopl + 1, 100).tolist()
    xticks = (np.sqrt(2) * np.array(xticks_base)).tolist()
    xticklabels = [f'{((xx + loopl) * resl) / 1e6:.2f}M' for xx in xticks_base]

    # --- Add custom ticks from annotation ---
    if annotation is not None:
        for label, bp in annotation.items():
            if ll <= bp < rr:
                # convert bp coordinate → diamond coordinate
                tick_pos = ((bp // resl) - loopl) * np.sqrt(2)
                # add tick + label
                xticks.append(tick_pos)
                xticklabels.append(label)

    # apply ticks to ALL shared axes
    plt.setp(axes, xticks=xticks)
    plt.setp(axes, xticklabels=xticklabels)

    fig.tight_layout()
    return fig, axes
