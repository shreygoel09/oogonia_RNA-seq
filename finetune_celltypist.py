import celltypist
import scanpy as sc
import pandas as pd
import torch
import argparse
import sys

from datetime import datetime


def pre_process(adata):
    # Combine immune and germ cell anndatas
    adata = sc.read_h5ad(adata)

    # Normalize anndatas
    sc.pp.log1p(adata)
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # obtain expression matrix from Anndata object
    combined_counts = pd.DataFrame.sparse.from_spmatrix(adata.X)
    combined_counts.columns = adata.var_names # assign gene names as column headers
    combined_counts.index = [cell for cell in adata.obs['cell_type']] # assign cell names as indices
    combined_counts = combined_counts.drop(combined_counts[combined_counts.index=='pre_spermatogonia'].index) # drop sperm cells
    combined_counts

    # check if normalization was done correctly
    sums = combined_counts.sum(axis=1).to_dict()
    un_normed = {}
    for k, v in sums.items():
        if abs(v - 10000) > 1:
            un_normed[k] = v
    print(f"Un-normalized cells: {un_normed}")
    sys.stdout.flush()

    # remove cells that were not normalized
    combined_counts.drop(combined_counts.index.intersection(un_normed.keys()), inplace=True)
    print(combined_counts)
    sys.stdout.flush()

    # save results
    combined_counts.to_csv("/workspace/a03-sgoel/Oogonia/CellTypist/Data/germ_immune_somatic_subsampled_10k.csv")


if __name__ == "__main__":
    print("processing args...")
    sys.stdout.flush()
    # Initialize argument parser
    parser = argparse.ArgumentParser(description="Process the args.")
    parser.add_argument('--anndata_path', type=str, help='cell by gene counts matrix')

    # Parse the arguments
    args = parser.parse_args()
    anndata = args.anndata_path

    # Check file paths
    print(f"Anndata path for pre-processing: {args.anndata_path}")
    sys.stdout.flush()

    # Check GPU
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    sys.stdout.flush()

    # pre-process data
    print("Pre-processing now...")
    sys.stdout.flush()
    pre_process(anndata)
    print("Pre-processing done!")
    sys.stdout.flush()
    counts = pd.read_csv("/workspace/a03-sgoel/Oogonia/CellTypist/Data/germ_immune_somatic_subsampled_10k.csv", index_col=0)
    print("Processed counts:")
    sys.stdout.flush()
    print(counts)
    sys.stdout.flush()

    # train model
    print("Starting training...")
    sys.stdout.flush()
    germ_immune_model = celltypist.train(counts,
                                   labels=counts.index.tolist(),
                                   genes=counts.columns.tolist(),
                                   check_expression=False,
                                   n_jobs=10,
                                   use_GPU=True,
                                   feature_selection=True)

    print("Finished training and saving model...")
    sys.stdout.flush()

    current_datetime = datetime.now()
    datetime_str = current_datetime.strftime('%m-%d-%Y-%H:%M:%S')
    germ_immune_model.write(f"/workspace/a03-sgoel/Oogonia/CellTypist/{datetime_str}_germ_immune_somatic_subsampled_10k_celltypist_model.pkl")
    print("Saved model!")
    sys.stdout.flush()