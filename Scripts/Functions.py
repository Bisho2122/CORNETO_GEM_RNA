import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from scipy.stats import norm

def average_expression_by_condition(expr_df, metadata_df):
    """
    Averages expression values per condition for each gene.

    Parameters:
    - expr_df: pd.DataFrame, with 'gene' column and sample columns (e.g., S1_Replic_1, etc.)
    - metadata_df: pd.DataFrame, with columns ['Sample_ID', 'Condition']

    Returns:
    - pd.DataFrame with 'gene' and one column per condition (average expression)
    """
    expr_long = expr_df.melt(id_vars='gene', var_name='Sample_ID', value_name='expression')
    merged = expr_long.merge(metadata_df, on='Sample_ID')
    grouped = merged.groupby(['gene', 'Condition'])['expression'].mean().reset_index()
    result = grouped.pivot(index='gene', columns='Condition', values='expression').reset_index()
    result.columns.name = None

    return result

def discretize_expression_gmm(df,
                              gene_col='gene', n_components=3, labels=['low', 'neutral', 'high']):
    """
    Discretizes gene expression using a Gaussian Mixture Model across all genes pooled.

    Parameters:
    - df: DataFrame with expression values and a gene column.
    - gene_col: Name of the column containing gene names.
    - n_components: Number of GMM components (typically 3).
    - labels: Labels corresponding to expression levels.

    Returns:
    - df_discrete: DataFrame with same shape but discretized values.
    - gmm: The fitted GaussianMixture model.
    """
    # Flatten values across all genes
    expr_values = df.drop(columns=gene_col).values.flatten().reshape(-1, 1)

    # Fit GMM
    gmm = GaussianMixture(n_components=n_components, random_state=42)
    gmm.fit(expr_values)

    # Get component assignments for all values
    component_assignments = gmm.predict(expr_values)

    # Map components to labels based on means (ascending → low to high)
    sorted_indices = np.argsort(gmm.means_.flatten())
    component_to_label = {comp: labels[i] for i, comp in enumerate(sorted_indices)}

    label_values = np.vectorize(component_to_label.get)(component_assignments)

    # Reshape to match original dataframe structure (excluding gene column)
    label_matrix = label_values.reshape(df.drop(columns=gene_col).shape)

    # Construct final discretized dataframe
    df_discrete = df.copy()
    df_discrete.update(pd.DataFrame(label_matrix, columns=df.columns[1:], index=df.index))

    return (df_discrete, gmm)

def plot_gmm_fit(gmm, data, bins=50):
    """
    Plot histogram of data with GMM component densities and total density overlay.

    Parameters:
    - gmm: fitted GaussianMixture object
    - data: 1D numpy array of data points (expression values)
    - bins: number of bins for histogram
    """
    x = np.linspace(data.min(), data.max(), 1000)
    logprob = gmm.score_samples(x.reshape(-1, 1))
    pdf = np.exp(logprob)

    # Plot histogram of data
    plt.hist(data, bins=bins, density=True, alpha=0.5, color='gray', label='Data histogram')

    # Plot total GMM density
    plt.plot(x, pdf, '-k', label='GMM total density')

    # Plot each component density
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_).flatten()
    weights = gmm.weights_.flatten()

    for mean, std, weight in zip(means, stds, weights):
        plt.plot(x, weight * norm.pdf(x, mean, std), '--', label=f'Component (mean={mean:.2f})')

    plt.xlabel('Expression')
    plt.ylabel('Density')
    plt.title('GMM fit over expression distribution')
    plt.legend()
    plt.show()



def plot_quantile_discretization(data,n_components=3, bins=50, labels=['low','neutral','high']):
    data = np.asarray(data).flatten()

    plt.hist(data, bins=bins, density=True, alpha=0.5, color='gray', label='Expression histogram')
    # Use pd.qcut to get bins
    _, bins = pd.qcut(data, q=n_components, retbins=True, labels=labels)
    thresholds = list(bins[1:-1])  # exclude min and max

    for thr in thresholds:
        plt.axvline(thr, color='blue', linestyle='--', label='Quantile threshold' if thr == thresholds[0] else None)

    plt.xlabel('Expression')
    plt.ylabel('Density')
    plt.title('Histogram with Quantile thresholds (pd.qcut)')
    plt.legend()
    plt.show()


def discretize_expression(df, method = "quantile", pool = True,
                          gene_col='gene', condition_col = "Condition", n_components=3, labels=['low', 'neutral', 'high']):
    df_long = df.melt(id_vars="gene", var_name="Condition", value_name="expression")
    if method == "quantile":
        if pool:
            expr_vals = df.drop(columns=gene_col).values.flatten()
            _,thresholds = pd.qcut(expr_vals, q=n_components, retbins=True, labels=labels)
            df_long["high"] = (df_long["expression"] >= thresholds[2]).astype(int)
            df_long["low"] = (df_long["expression"] <= thresholds[1]).astype(int)
            df_long["score"] = df_long["high"] - df_long["low"]
        else:
            all_conds = df_long["Condition"].unique()
            for cond in all_conds:
                cond_df = df_long[df_long["Condition"] == cond]
                expr_vals = cond_df.drop(columns=["gene", "Condition"]).values.flatten()
                _, thresholds = pd.qcut(expr_vals, q=n_components, retbins=True, labels=labels)
                df_long.loc[df_long["Condition"] == cond, 'high'] = (df_long.loc[df_long["Condition"] == cond, 'expression'] >= thresholds[2]).astype(int)
                df_long.loc[df_long["Condition"] == cond, 'low'] = (df_long.loc[df_long["Condition"] == cond, 'expression'] <= thresholds[1]).astype(int)
                df_long.loc[df_long["Condition"] == cond, 'score'] = df_long.loc[df_long["Condition"] == cond, 'high'] - df_long.loc[df_long["Condition"] == cond, 'low']
        #plot_quantile_discretization(expr_vals, n_components=n_components, bins=50, labels=labels)
        return df_long
    elif method == "gmm":
        if pool:
            gmm_res = discretize_expression_gmm(df = df, gene_col=gene_col, n_components=n_components, labels=labels)
            df_discrete = (
                gmm_res[0]
                .melt(id_vars="gene", var_name="Condition", value_name="label")
                .pivot(index=["gene", "Condition"], columns="label", values="label")
                .drop(columns="neutral")
                .fillna(0)
                .reset_index()
            )
            df_discrete['high'][df_discrete['high'] == 'high'] = 1
            df_discrete['low'][df_discrete['low'] == 'low'] = 1
            df_discrete["score"] = df_discrete["high"] - df_discrete["low"]
            df_long = df_long.merge(df_discrete, on=["gene", "Condition"], how="left")
            final_results = {"updated_df": df_long, "gmm_results": gmm_res}
            return final_results
        else:
            all_conds = df_long["Condition"].unique()
            df_long['high'] = 0
            df_long['low'] = 0
            df_long['score'] = 0
            all_gmm_res = {}
            for cond in all_conds:
                cond_df = df_long[df_long["Condition"] == cond].drop(columns = ["Condition", "high", "low", "score"])
                gmm_res = discretize_expression_gmm(df=cond_df, gene_col=gene_col, n_components=n_components, labels=labels)
                df_discrete = (
                    gmm_res[0]
                    .melt(id_vars="gene", var_name="Condition", value_name="label")
                    .pivot(index=["gene", "Condition"], columns="label", values="label")
                    .drop(columns="neutral")
                    .fillna(0)
                    .reset_index()
                )
                df_discrete['high'][df_discrete['high'] == 'high'] = 1
                df_discrete['low'][df_discrete['low'] == 'low'] = 1
                df_discrete["score"] = df_discrete["high"] - df_discrete["low"]

                df_long.loc[df_long["Condition"] == cond, ['high', 'low', 'score']] = df_discrete[['high', 'low', 'score']].values
                all_gmm_res[cond] = gmm_res
            final_results = {"updated_df": df_long, "gmm_results": all_gmm_res}
            return final_results
    else:
        raise ValueError("Method must be either 'quantile' or 'gmm'.")