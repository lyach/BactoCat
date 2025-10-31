import cobra
import numpy as np
from .fvfa import FVAProblem


def cobra_to_fva_problem(model, mu):
    from cobra.util.array import create_stoichiometric_matrix

    S = create_stoichiometric_matrix(model)
    v_l = np.array([rxn.lower_bound for rxn in model.reactions])
    v_u = np.array([rxn.upper_bound for rxn in model.reactions])
    objective_rxn = list(model.objective.variables)[0].name
    c = np.array([1 if rxn.id == objective_rxn else 0 for rxn in model.reactions])

    return FVAProblem(S=S, v_l=v_l, v_u=v_u, c=c, mu=mu)
