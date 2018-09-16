from statsmodels.discrete.discrete_model import Logit
import numpy as np


def check_is_diplo(result):
    return result.total == result.A_count+result.B_count


def simple_model(motif_results_A, non_results_A, motif_results_B, non_results_B):
    all_results = motif_results_A + non_results_A + motif_results_B + non_results_B
    is_diplo = np.array([check_is_diplo(result) for result in all_results], dtype="int")
    total_gpc = (len(motif_results_A)+len(non_results_A))
    is_gpc = np.zeros(len(all_results), dtype="int")
    is_gpc[:total_gpc] = 1
    motif = np.zeros(len(all_results), dtype="int")
    motif[:len(motif_results_A)] = 1
    motif[total_gpc:total_gpc+len(motif_results_B)] = 1
    X = np.hstack((np.ones_like(is_gpc)[:, None], is_gpc[:, None], is_diplo[:, None], (is_gpc*is_diplo)[:, None]))
    print(np.sum(X*motif[:, None], axis=0)/np.sum(X, axis=0))
    print(np.sum((1-X)*motif[:, None], axis=0)/np.sum((1-X), axis=0))
    y = motif
    model = Logit(y, X)
    result = model.fit()
    print(result.summary())
