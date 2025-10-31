from enum import Enum

import numpy
import numpy as np
from dataclasses import dataclass
from typing import Optional, List, Tuple

from docplex.mp.model import Model
from docplex.mp.constants import OptimizationStatus


@dataclass
class FVAProblem:
    r"""
    This is the class object that defines FVA problems in this package

    .. math:: \min c^Tv

    Subject to:
        Sv = 0
        v_l <= v <= v_u
        v \in R^n
    """
    S: numpy.array
    v_u: numpy.array
    v_l: numpy.array
    c: numpy.array
    mu: float

    def __post_init__(self):
        """Some clean-up, and error checking is done here"""

        # assure that v_u and v_l are flattened arrays
        self.v_l = self.v_l.flatten()
        self.v_u = self.v_u.flatten()

        if not self.is_valid_dims():
            RuntimeError('Dimensions of vertices not correct! Please check the dims of the problems!')
        if not self.is_valid_param():
            RuntimeError('The mu parameter must be between [0,1]')

    def is_valid_dims(self) -> bool:
        """
        Checks the dimensions of the FVA problem matrices.

        :return: False if not consistent dimensions, otherwise true
        """

        num_lb = numpy.size(self.v_l)
        num_ub = numpy.size(self.v_u)
        num_c = numpy.size(self.c)

        # check upper bounds and lower bounds
        if num_lb != num_ub:
            return False

        # check the objective dims
        if (num_lb != num_c) or (num_ub != num_c):
            return False

        # check the dims of the metabolic network matrix
        if self.S.shape[1] != num_lb:
            return False

        return True

    def is_valid_param(self) -> bool:
        return 0 <= self.mu <= 1

    def num_v(self) -> int:
        return numpy.size(self.c)


@dataclass
class FVASolution:
    """
    The solution to FVAProblem instance. Contains the bounds of the fluxs, the best objective value, the number of
    LPs needed to solve and the FVAProblem instance.
    """
    lower_bound: numpy.array
    upper_bound: numpy.array
    Z: float
    number_LPs: int
    problem: FVAProblem


def build_fva_lp(problem: FVAProblem) -> Tuple[Model, List]:
    """
    Build initial CPLEX model via docplex.
    Returns (model, var_list)
    """
    # create model (silenced)
    model = Model(name="fva_model")
    model.context.solver.suppress_output = True

    n = problem.num_v()
    lb = list(problem.v_l.flatten())
    ub = list(problem.v_u.flatten())

    # create continuous vars with bounds
    v = model.continuous_var_list(n, lb=lb, ub=ub, name="v")

    # add stoichiometric constraints: S @ v == 0
    S = problem.S
    # S is numpy array rows = metabolites, cols = reactions
    for i in range(S.shape[0]):
        row = S[i, :]
        # build linear expression only for nonzero entries
        expr = model.sum(float(row[j]) * v[j] for j in range(n) if row[j] != 0)
        model.add_constraint(expr == 0, ctname=f"flux_network_{i}")

    # set objective maximize c^T v
    c = problem.c.flatten()
    obj_expr = model.sum(float(c[j]) * v[j] for j in range(n))
    model.maximize(obj_expr)

    return model, v


def setup_initial_fva_problem_solve(problem: FVAProblem):
    """
    Build model, solve initial optimization to obtain Z and add percent-optimality constraint.
    Returns (model, var_list, Z, v_vals_initial)
    """
    fva_model, flux_vars = build_fva_lp(problem)

    sol = fva_model.solve()
    if sol is None:
        raise RuntimeError("Initial optimization infeasible or solver failed")

    # get objective value
    Z = float(fva_model.objective_value)

    # fetch the solution variable values as a numpy array
    # fva_model.solution is the solution object after solve()
    sol_obj = fva_model.solution
    v_vals_initial = np.array([float(sol_obj.get_value(var)) for var in flux_vars])

    # Add percent optimality constraint: c^T v >= mu * Z
    c = problem.c.flatten()
    n = problem.num_v()
    expr = fva_model.sum(float(c[j]) * flux_vars[j] for j in range(n))
    fva_model.add_constraint(expr >= float(problem.mu * Z), ctname="Percent_Optimality")

    return fva_model, flux_vars, float(Z), v_vals_initial


def _get_flux_values_from_model(model: Model, flux_vars: List) -> numpy.ndarray:
    """
    Helper to extract the current solution flux values from a model.
    Returns numpy array of floats; if no solution available returns array of nans.
    """
    sol = model.solution
    if sol is None:
        return numpy.full(len(flux_vars), numpy.nan)
    # docplex solution.get_value(var) works for Var objects
    vals = numpy.array([float(sol.get_value(v)) for v in flux_vars], dtype=float)
    return vals


def find_variable_range(model: Model, flux_vars, var_index: int, sense: str):
    """
    Replaces the objective with flux_vars[var_index] and resolves.
    Returns tuple: (objective_value, v_vals) where v_vals is numpy array of current flux values.
    sense: 'MAX' or 'MIN'
    """
    target = flux_vars[var_index]
    if sense.upper() == "MAX":
        model.maximize(target)
    else:
        model.minimize(target)

    sol = model.solve()
    if sol is None:
        # return NaN on infeasible, and NaN vector
        n = len(flux_vars)
        return float("nan"), np.full(n, np.nan)

    sol_obj = model.solution
    obj_val = float(model.objective_value)
    v_vals = np.array([float(sol_obj.get_value(var)) for var in flux_vars])

    return obj_val, v_vals

def fva_solve_basic(problem: FVAProblem) -> Optional[FVASolution]:
    """
    Solve the FVA problem with the standard algorithm

    :param problem: A FVAProblem type
    :return: Optionally None if the problem is not feasible, otherwise a FVASolution object
    """

    # build to solve the first step of the FVA problem
    fva_model, flux_vars, Z = setup_initial_fva_problem_solve(problem)

    # instantiate the flux bounds and the LP counting
    lower_bound = numpy.zeros(problem.num_v())
    upper_bound = numpy.zeros(problem.num_v())
    number_lps = 1 + 2 * problem.num_v()

    # solve the FVA sub problems
    for index in range(problem.num_v()):
        upper_bound[index] = find_variable_range(fva_model, flux_vars, index, "MAX")
        lower_bound[index] = find_variable_range(fva_model, flux_vars, index, "MIN")

    return FVASolution(lower_bound, upper_bound, Z, number_lps, problem)


def fva_solve_basic_parallel(problem: FVAProblem) -> Optional[FVASolution]:
    """
    Solve the FVA problem with the standard algorithm (parallel-friendly structure, sequential here)

    :param problem: A FVAProblem type
    :return: Optionally None if the problem is not feasible, otherwise a FVASolution object
    """

    # build to solve the first step of the FVA problem
    fva_model, flux_vars, Z = setup_initial_fva_problem_solve(problem)

    # instantiate the flux bounds and the LP counting
    lower_bound = numpy.zeros(problem.num_v())
    upper_bound = numpy.zeros(problem.num_v())
    number_lps = 1 + 2 * problem.num_v()

    # solve the FVA sub problems
    for index in range(problem.num_v()):
        upper_bound[index] = find_variable_range(fva_model, flux_vars, index, "MAX")
        lower_bound[index] = find_variable_range(fva_model, flux_vars, index, "MIN")

    return FVASolution(lower_bound, upper_bound, Z, number_lps, problem)


def fva_solve_faster(problem: FVAProblem) -> Optional[FVASolution]:
    # build the solve the first step of the FVA problem
    fva_model, flux_vars, Z, initial_v_vals = setup_initial_fva_problem_solve(problem)

    # set up problems to be solved
    upper_bound_problems = set(range(problem.num_v()))
    lower_bound_problems = set(range(problem.num_v()))

    # instantiate the flux bounds and the LP counting
    lower_bound = np.empty(problem.num_v())
    lower_bound.fill(np.nan)

    upper_bound = np.empty(problem.num_v())
    upper_bound.fill(np.nan)

    number_lps = 1

    def remove_bound_problems(v_val_array, atol=1e-9):
        """Helper function to remove problems to solve from the lower and upper bound sets.
           v_val_array is a numpy array of variable values from the solver.
        """
        # indices where variable equals lower bound (within tolerance)
        lb_mask = np.isclose(v_val_array, problem.v_l, atol=atol)
        ub_mask = np.isclose(v_val_array, problem.v_u, atol=atol)

        lb_idx = np.where(lb_mask)[0]
        ub_idx = np.where(ub_mask)[0]

        lower_bound_problems.difference_update(set(lb_idx))
        upper_bound_problems.difference_update(set(ub_idx))

    # initial pruning using the result from the initial solve
    remove_bound_problems(initial_v_vals)

    # solve the upper bound problems
    while len(upper_bound_problems) != 0:
        ub_problem = upper_bound_problems.pop()
        obj_val, v_vals = find_variable_range(fva_model, flux_vars, ub_problem, "MAX")
        upper_bound[ub_problem] = obj_val
        number_lps += 1

        # prune based on this solve
        remove_bound_problems(v_vals)

    # solve the lower bound problems
    while len(lower_bound_problems) != 0:
        lb_problem = lower_bound_problems.pop()
        obj_val, v_vals = find_variable_range(fva_model, flux_vars, lb_problem, "MIN")
        lower_bound[lb_problem] = obj_val
        number_lps += 1

        # prune based on this solve
        remove_bound_problems(v_vals)

    # fill remaining NaNs with problem bounds
    for index in range(problem.num_v()):
        if np.isnan(upper_bound[index]):
            upper_bound[index] = problem.v_u[index]
        if np.isnan(lower_bound[index]):
            lower_bound[index] = problem.v_l[index]

    return FVASolution(lower_bound, upper_bound, Z, number_lps, problem)
