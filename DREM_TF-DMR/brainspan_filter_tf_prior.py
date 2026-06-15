#!/usr/bin/env python3
"""
Filter a TF-gene prior file by BrainSpan bulk RNA-seq expression.

Keeps only TFs whose mean RPKM across all BrainSpan samples exceeds the
threshold. Works on both collapsed priors (first column = TF name) and
annotated priors (first column = TF|cell|SRX); the TF name is extracted
as the portion before the first pipe.

Outputs the filtered prior alongside diagnostic files:
    - tf_mean_rpkm.tsv              per-TF BrainSpan mean RPKM
    - tfs_kept.txt                  TFs that passed the filter
    - tfs_dropped_below_thresh.txt  TFs present in BrainSpan but below threshold
    - tfs_dropped_not_in_brainspan  TFs absent from BrainSpan (e.g. chemical
                                    antibody targets like 8-Hydroxydeoxyguanosine)

Example invocation on Hoffman2:

    module load python/3.9.6
    python3 brainspan_filter_tf_prior.py \\
        --expr /u/project/cluo/chongyua/brain_dev_snm3C/2022/analysis/brainspan/expression_matrix.csv \\
        --meta /u/project/cluo/chongyua/brain_dev_snm3C/2022/analysis/brainspan/rows_metadata.csv \\
        --prior ~/epigenome/bican/TF_gene_prior_bican_top3_new.3col.txt \\
        --outdir ~/epigenome/bican/brainspan_filtered \\
        --mean-thresh 1.0
"""

import argparse
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--expr", required=True,
                    help="BrainSpan expression matrix CSV (first col = row_num)")
    ap.add_argument("--meta", required=True,
                    help="BrainSpan rows_metadata.csv (row_num, gene_symbol, ...)")
    ap.add_argument("--prior", required=True,
                    help="TF-gene prior file in 3-col format (TF, Gene, Input)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--mean-thresh", type=float, default=1.0,
                    help="Mean RPKM threshold (default: 1.0)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[1/6] Loading BrainSpan expression matrix from {args.expr}")
    expr = pd.read_csv(args.expr, low_memory=False, header=None)
    expr = expr.rename(columns={0: "row_num"})
    expr["row_num"] = pd.to_numeric(expr["row_num"], errors="coerce")
    expr = expr.dropna(subset=["row_num"]).astype({"row_num": int}).set_index("row_num")
    n_samples = expr.shape[1]
    print(f"       {expr.shape[0]} rows x {n_samples} samples")

    print(f"[2/6] Loading BrainSpan metadata from {args.meta}")
    meta = pd.read_csv(args.meta)
    meta = meta[["row_num", "gene_symbol"]].dropna()
    meta["row_num"] = meta["row_num"].astype(int)
    meta["gene_symbol"] = meta["gene_symbol"].astype(str).str.strip()
    print(f"       {len(meta)} gene rows")

    print(f"[3/6] Computing per-gene mean RPKM across {n_samples} samples")
    gene_mean = expr.mean(axis=1).rename("mean_rpkm")
    # Join with gene symbol; some genes have multiple transcript rows,
    # so take the max mean per symbol (a TF counts as expressed if any of
    # its transcript rows meets the threshold).
    joined = gene_mean.to_frame().join(meta.set_index("row_num"), how="inner")
    gene_mean_by_symbol = joined.groupby("gene_symbol")["mean_rpkm"].max()
    print(f"       {len(gene_mean_by_symbol)} unique gene symbols with expression")

    print(f"[4/6] Loading TF-gene prior from {args.prior}")
    prior = pd.read_csv(args.prior, sep="\t", dtype=str)
    print(f"       {len(prior)} rows, "
          f"{prior.iloc[:, 0].nunique()} unique TF-column entries")

    # Extract the "TF name" from the first column. If the column is
    # annotated as "TF|cell|SRX", take the part before the first pipe.
    # If it's already collapsed (just "TF"), the split is a no-op.
    tf_col = prior.columns[0]
    tf_symbol = prior[tf_col].str.split("|").str[0].str.strip()
    unique_tfs = tf_symbol.unique()
    print(f"       {len(unique_tfs)} unique TF symbols after collapsing")

    print(f"[5/6] Applying filter: mean RPKM > {args.mean_thresh}")
    tf_expr = pd.Series(
        {tf: gene_mean_by_symbol.get(tf, float("nan")) for tf in unique_tfs}
    )
    keep_mask = tf_expr > args.mean_thresh
    keep_tfs = set(tf_expr[keep_mask].index)
    missing_tfs = set(tf_expr[tf_expr.isna()].index)
    below_thresh_tfs = set(tf_expr[~keep_mask].index) - keep_tfs - missing_tfs

    print(f"       Kept:                  {len(keep_tfs)} TFs")
    print(f"       Dropped (below thresh): {len(below_thresh_tfs)} TFs")
    print(f"       Dropped (not in BS):   {len(missing_tfs)} TFs")

    print(f"[6/6] Writing outputs to {args.outdir}")
    prior_filtered = prior[tf_symbol.isin(keep_tfs)].copy()

    filt_path = os.path.join(
        args.outdir,
        f"{os.path.basename(args.prior).rsplit('.', 1)[0]}"
        f".mean_gt_{args.mean_thresh}.tsv",
    )
    prior_filtered.to_csv(filt_path, sep="\t", index=False)
    print(f"       Filtered prior:          {filt_path}  ({len(prior_filtered)} rows)")

    expr_path = os.path.join(args.outdir, "tf_mean_rpkm.tsv")
    tf_expr.to_frame("mean_rpkm").to_csv(expr_path, sep="\t", index_label="TF")
    print(f"       TF mean RPKM table:      {expr_path}")

    for name, tfs in [("kept", keep_tfs),
                      ("dropped_below_thresh", below_thresh_tfs),
                      ("dropped_not_in_brainspan", missing_tfs)]:
        p = os.path.join(args.outdir, f"tfs_{name}.txt")
        with open(p, "w") as f:
            for tf in sorted(tfs):
                rpkm = tf_expr.get(tf, float("nan"))
                line = f"{tf}\tNA\n" if pd.isna(rpkm) else f"{tf}\t{rpkm}\n"
                f.write(line)
        print(f"       {name:27s}: {p}")


if __name__ == "__main__":
    main()
