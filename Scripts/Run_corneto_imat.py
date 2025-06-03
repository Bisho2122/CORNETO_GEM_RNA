import corneto as cn
import pandas as pd
import numpy as np
import argparse
import os
from pandas.errors import SettingWithCopyWarning


import Functions as custom

import warnings

from corneto.methods.metabolism import evaluate_gpr_expression
from corneto._data import Data
from corneto.methods.future.imat import MultiSampleIMAT

warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", SettingWithCopyWarning)

parser = argparse.ArgumentParser(description="Run multisample iMAT with CORNETO on RNA-Seq data.")

parser.add_argument("-rna", "--rna_norm", type=str, help="Path to the RNA-Seq Normalized matrix")
parser.add_argument("-meta", "--metadata", type=str, help="Path to RNA-Seq metadata file")
parser.add_argument("-gem", "--model", type=str, help="Path to GEM model")
parser.add_argument("-disc", "--discr_method", type=str, default= "gmm", help="Discretization method [quantile or gmm] (default: 'gmm')")
parser.add_argument("-o", "--out_dir", type=str, help="Path to the output directory")

args = parser.parse_args()


G = cn.io.import_cobra_model(args.model)
df_expr = pd.read_csv(args.rna_norm)
metadata = pd.read_csv(args.metadata)

dataset = custom.average_expression_by_condition(df_expr, metadata)

discretized = custom.discretize_expression(df = dataset, method=args.discr_method, pool=False)

if(args.discr_method == "quantile"):
    conds = discretized['Condition'].unique()
    data = {}
    for cond in conds:
        cond_df = discretized[discretized['Condition'] == cond]
        cond_df.score = cond_df.score.astype(int)
        values = cond_df[["gene", "score"]].set_index("gene").to_dict()["score"]
        e = evaluate_gpr_expression(G.get_attr_from_edges("GPR"), values)

        feats = dict()
        for edge, val in zip(G.get_attr_from_edges("id"), e):
            feats[edge] = {"value": val, "mapping": "edge"}
        data[cond] = feats
else:
    conds = discretized['updated_df']['Condition'].unique()
    data = {}
    for cond in conds:
        cond_df = discretized['updated_df'][discretized['updated_df']['Condition'] == cond]
        cond_df.score = cond_df.score.astype(int)
        values = cond_df[["gene", "score"]].set_index("gene").to_dict()["score"]
        e = evaluate_gpr_expression(G.get_attr_from_edges("GPR"), values)

        feats = dict()
        for edge, val in zip(G.get_attr_from_edges("id"), e):
            feats[edge] = {"value": val, "mapping": "edge"}
        data[cond] = feats

data = Data.from_cdict(data)

# eps is the min absolute flux value that we will consider
# for reactions with positive scores, i.e., if a reaction
# with a positive score is included in the context specific
# metabolic network, it has to be able to carry an absolute
# metabolic flux value above eps.
eps = 1e-2

# We will also add a small penalty to make networks sparser.
# When we have more than 1 sample, this parameter controls
# the regularization across samples, blocking entire groups
# of reactions across samples. Here, we are just penalizing
# the inclusion of reactions not needed to fit the reaction
# scores
m = MultiSampleIMAT(lambda_reg=1e-3, eps=eps, use_bigm_constraints=False)
preprocessed_data = m.preprocess(G, data)
P = m.build(preprocessed_data[0], preprocessed_data[1])

P.solve(solver="scipy", verbosity=0)

df_sol = pd.DataFrame(
    P.expr.flow.value, index=G.get_attr_from_edges("id"), columns=[i + '_flux' for i in list(data.samples.keys())]
)

df_sol.to_csv(os.path.join(args.out_dir, "imat_fluxes.csv"))

print("iMAT fluxes saved to:", os.path.join(args.out_dir, "imat_fluxes.csv"))