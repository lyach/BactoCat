"""
Flux Reconciliation Module

This module performs flux reconciliation using genome-scale metabolic models
and experimental flux measurements.
"""

import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition
import pandas as pd
from typing import Dict, List, Tuple, Optional
import warnings
import cobra
from cobra.io import read_sbml_model


def load_gem_model(model_path: str) -> cobra.Model:
    """
    Load a genome-scale metabolic model from file path.
    
    Parameters:
    -----------
    model_path : str
        Path to the model file (SBML format)
        
    Returns:
    --------
    cobra.Model
        The loaded COBRA model
    """
    try:
        model = read_sbml_model(model_path)
        print(f"Model loaded successfully from {model_path}")
        print(f"Model contains {len(model.reactions)} reactions and {len(model.metabolites)} metabolites")
        return model
    except Exception as e:
        raise ValueError(f"Failed to load model from {model_path}: {str(e)}")


def extract_model_info(model: cobra.Model) -> Tuple[np.ndarray, List[str], List[str], int, int, Dict]:
    """
    Extract stoichiometric matrix and other information from COBRA model.
    
    Parameters:
    -----------
    model : cobra.Model
        The loaded COBRA model
        
    Returns:
    --------
    tuple
        (S_matrix, metabolite_names, reaction_names, n_metabolites, n_reactions, flux_bounds)
    """
    # Get reaction and metabolite names
    reaction_names = [rxn.id for rxn in model.reactions]
    metabolite_names = [met.id for met in model.metabolites]
    
    n_reactions = len(reaction_names)
    n_metabolites = len(metabolite_names)
    
    # Create stoichiometric matrix
    S = np.zeros((n_metabolites, n_reactions))
    
    for j, reaction in enumerate(model.reactions):
        for metabolite, coefficient in reaction.metabolites.items():
            i = metabolite_names.index(metabolite.id)
            S[i, j] = coefficient
    
    # Extract flux bounds
    flux_bounds = {}
    for i, reaction in enumerate(model.reactions):
        flux_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
    
    # Debug prints
    print("\n=== Model Information ===")
    print(f"Number of reactions: {n_reactions}")
    print(f"Number of metabolites: {n_metabolites}")
    print(f"S matrix shape: {S.shape}")
    print(f"Number of non-zero elements in S: {np.count_nonzero(S)}")
    print("Sample of flux bounds:")
    for rxn_id in list(flux_bounds.keys())[:5]:
        print(f"{rxn_id}: {flux_bounds[rxn_id]}")
    
    return S, metabolite_names, reaction_names, n_metabolites, n_reactions, flux_bounds


def process_measured_fluxes_df(measured_fluxes_df: pd.DataFrame, reaction_names: List[str]) -> Dict:
    """
    Process measured fluxes dataframe to extract measurements.
    
    Parameters:
    -----------
    measured_fluxes_df : pd.DataFrame
        DataFrame with columns: ['rxn_id', 'exp_flux']
    reaction_names : List[str]
        List of reaction names from the model
        
    Returns:
    --------
    dict
        Dictionary of measurements
    """
    measurements = {}
    missing_reactions = []
    
    for _, row in measured_fluxes_df.iterrows():
        rxn_id = row['rxn_id']
        if rxn_id in reaction_names:
            measurements[rxn_id] = row['exp_flux']
        else:
            missing_reactions.append(rxn_id)
            warnings.warn(f"Reaction {rxn_id} not found in model reactions")
    
    # Debug prints
    print("\n=== Measured Fluxes Information ===")
    print(f"Total measured reactions: {len(measured_fluxes_df)}")
    print(f"Matched measurements: {len(measurements)}")
    print(f"Missing reactions: {len(missing_reactions)}")
    if missing_reactions:
        print("First few missing reactions:", missing_reactions[:5])
    print("\nSample of measurements:")
    for rxn_id in list(measurements.keys())[:5]:
        print(f"{rxn_id}: flux = {measurements[rxn_id]}")
    
    return measurements


def process_fixed_fluxes_df(fixed_fluxes_df: pd.DataFrame, reaction_names: List[str], measurements: Optional[Dict] = None) -> Dict:
    """
    Process fixed fluxes dataframe to extract fixed flux values.
    
    Parameters:
    -----------
    fixed_fluxes_df : pd.DataFrame
        DataFrame with columns: ['rxn_id', 'exp_fixed_flux']
    reaction_names : List[str]
        List of reaction names from the model
    measurements : Optional[Dict]
        Dictionary of measured fluxes to check for consistency
        
    Returns:
    --------
    dict
        Dictionary of fixed fluxes
    """
    fixed_fluxes = {}
    missing_reactions = []
    conflicts = []
    
    # Critical exchange reactions that shouldn't be fixed to zero
    critical_exchanges = ['EX_glc__D_e']
    
    for _, row in fixed_fluxes_df.iterrows():
        rxn_id = row['rxn_id']
        flux_value = row['exp_fixed_flux']
        
        # Skip if reaction not in model
        if rxn_id not in reaction_names:
            missing_reactions.append(rxn_id)
            continue
            
        # Check if fixing to zero for critical exchanges
        if rxn_id in critical_exchanges and abs(flux_value) < 1e-10:
            print(f"\nWarning: {rxn_id} has zero flux. Might lead to infeasibility.")
            continue
            
        # Check for conflicts with measurements
        if measurements and rxn_id in measurements:
            measured_value = measurements[rxn_id]
            if abs(measured_value - flux_value) > 1e-6:  # Allow small differences
                conflicts.append((rxn_id, measured_value, flux_value))
                print(f"\nWarning: Conflict for {rxn_id}:")
                print(f"  Measured value: {measured_value}")
                print(f"  Fixed value: {flux_value}")
                print("Using fixed value and not measured value.")
                
        fixed_fluxes[rxn_id] = flux_value
    
    # Debug prints
    print("\n=== Fixed Fluxes Information ===")
    print(f"Total fixed reactions: {len(fixed_fluxes_df)}")
    print(f"Matched fixed fluxes: {len(fixed_fluxes)}")
    print(f"Missing reactions: {len(missing_reactions)}")
    if missing_reactions:
        print("Missing reactions:", missing_reactions)
    if conflicts:
        print("\nConflicting reactions:")
        for rxn, meas, fix in conflicts:
            print(f"{rxn}: measured={meas}, fixed={fix}")
    print("\nAll fixed fluxes:")
    for rxn_id, value in fixed_fluxes.items():
        print(f"{rxn_id}: {value}")
    
    return fixed_fluxes


def build_flux_reconciliation_model(S: np.ndarray, 
                                  obj_rxn: str,
                                  metabolite_names: List[str],
                                  reaction_names: List[str],
                                  measurements: Dict,
                                  fixed_fluxes: Dict,
                                  flux_bounds: Dict) -> pyo.ConcreteModel:
    """
    Build the Pyomo optimization model for flux reconciliation.
    
    Parameters:
    -----------
    S : np.ndarray
        Stoichiometric matrix
    metabolite_names : List[str]
        List of metabolite names
    reaction_names : List[str]
        List of reaction names
    measurements : Dict
        Dictionary of measured fluxes
    fixed_fluxes : Dict
        Dictionary of fixed fluxes
    flux_bounds : Dict
        Dictionary of flux bounds
        
    Returns:
    --------
    pyo.ConcreteModel
        The constructed Pyomo model
    """
    n_metabolites = len(metabolite_names)
    n_reactions = len(reaction_names)
    
    # Set flux bounds
    bounds = {}
    default_lower = -1000.0
    default_upper = 1000.0

    for i, reaction_name in enumerate(reaction_names):
        if reaction_name in flux_bounds:
            bounds[i] = flux_bounds[reaction_name]
        else:
            bounds[i] = (default_lower, default_upper)
    
    # Debug prints for bounds
    print("\n=== Optimization Model Setup ===")
    print("Sample of bounds:")
    for i in range(min(5, len(reaction_names))):
        print(f"{reaction_names[i]}: {bounds[i]}")
    
    # Initialize Pyomo model
    model = pyo.ConcreteModel()
            
    # Sets
    model.metabolites = pyo.Set(initialize=range(n_metabolites))
    model.reactions = pyo.Set(initialize=range(n_reactions))
    
    # Dual variables for optimiality conditions
    model.alpha_L = pyo.Var(model.reactions, bounds=(0, 1000))    # Lower bound multipliers
    model.alpha_U = pyo.Var(model.reactions, bounds=(0, 1000))    # Upper bound multipliers
    
    # Split measured reactions into zero and non-zero measurements
    zero_measurements = [
        j for j, name in enumerate(reaction_names)
        if name in measurements and abs(measurements[name]) < 1e-10
    ]
    nonzero_measurements = [
        j for j, name in enumerate(reaction_names)
        if name in measurements and abs(measurements[name]) >= 1e-10
    ]
    
    model.zero_measured_reactions = pyo.Set(initialize=zero_measurements)
    model.nonzero_measured_reactions = pyo.Set(initialize=nonzero_measurements)
    
    print("\n=== Measurement Information ===")
    print(f"Number of zero measurements: {len(zero_measurements)}")
    print(f"Number of non-zero measurements: {len(nonzero_measurements)}")

    # Variables
    model.v = pyo.Var(model.reactions, bounds=(-1000, 1000))

    # Apply flux bounds to reaction variables
    for i in model.reactions:
        model.v[i].setlb(bounds[i][0])
        model.v[i].setub(bounds[i][1])
        
        
    # Fixed fluxes (discarded)
    # if fixed_fluxes:
    #     print("\n=== Fixed Fluxes in Model ===")
    #     for reaction_name, flux_value in fixed_fluxes.items():
    #         if reaction_name in reaction_names:
    #             reaction_idx = reaction_names.index(reaction_name)
    #             model.v[reaction_idx].fix(flux_value)
    #             print(f"Fixed {reaction_name} to {flux_value}")
        
    # Fixed fluxes as bounds
    if fixed_fluxes:
        print("\n=== Setting Flux Bounds in Model ===")
    for reaction_name, flux_value in fixed_fluxes.items():
        if reaction_name in reaction_names:
            reaction_idx = reaction_names.index(reaction_name)
            v = model.v[reaction_idx]
            if flux_value < 0:
                v.setlb(flux_value)  # uptake
                print(f"Set lower bound of {reaction_name} to {flux_value}")
            else:
                v.setub(flux_value)  # secretion
                print(f"Set upper bound of {reaction_name} to {flux_value}")
                
    # Debug prints for model structure
    print("\n=== Model Structure ===")
    print(f"Number of variables: {len(model.v)}")
    print(f"Number of fixed fluxes: {len(fixed_fluxes)}")
                
    # Objective function: combining relative error for non-zero measurements
    # and absolute error for zero measurements
    def objective_rule(m):
        # For non-zero measurements, use relative error
        relative_error = sum(
            ((m.v[j] - measurements[reaction_names[j]])/measurements[reaction_names[j]])**2  
            for j in m.nonzero_measured_reactions
        )
        
        # For zero measurements, use absolute error
        absolute_error = sum(
            (m.v[j])**2  # Since measurement is zero
            for j in m.zero_measured_reactions
        )

        return relative_error + absolute_error
    model.obj = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

    # Mass balance constraints (S matrix)
    def mass_balance_rule(m, i):
        return sum(S[i, j] * m.v[j] for j in m.reactions) == 0
    model.mass_balance = pyo.Constraint(model.metabolites, rule=mass_balance_rule)
    print(f"Number of mass balance constraints: {len(model.mass_balance)}")
    
    # Parameters
    barrier_param = 1e-6
    
    # Vector of initial weights (zeros). Each entry corresponds to a reaction
    objective_weights = np.zeros(n_reactions)

    # Assign weight of 1.0 to reaction of objective function
    objective_reaction = obj_rxn
    reaction_idx = reaction_names.index(objective_reaction)
    objective_weights[reaction_idx] = 1.0

    # Dual feasibility constraint
    def stationarity_rule(m, j):
        return (objective_weights[j] + m.alpha_L[j] - m.alpha_U[j] == 0)
    model.stationarity = pyo.Constraint(model.reactions, rule=stationarity_rule)
        
    # Complementarity constraints
    def comp_lower_rule(m, j):
        return (m.v[j] - bounds[j][0]) * m.alpha_L[j] <= barrier_param
    model.comp_lower = pyo.Constraint(model.reactions, rule=comp_lower_rule)
            
    def comp_upper_rule(m, j):
        return (bounds[j][1] - m.v[j]) * m.alpha_U[j] <= barrier_param
    model.comp_upper = pyo.Constraint(model.reactions, rule=comp_upper_rule)
    
    return model


def solve_flux_reconciliation_model(model: pyo.ConcreteModel,
                                  reaction_names: List[str],
                                  measurements: Dict,
                                  solver_name: str = 'ipopt',
                                  solver_options: Optional[Dict] = None) -> Dict:
    """
    Solve the flux reconciliation optimization model.
    
    Parameters:
    -----------
    model : pyo.ConcreteModel
        The Pyomo model to solve
    reaction_names : List[str]
        List of reaction names
    measurements : Dict
        Dictionary of measured fluxes
    solver_name : str
        Name of the solver to use
    solver_options : Dict
        Dictionary of solver options
        
    Returns:
    --------
    Dict
        Results dictionary
    """
    if solver_options is None:
        solver_options = {
            'max_iter': 3000,
            'tol': 1e-8,
            'acceptable_tol': 1e-6
        }

    solver = pyo.SolverFactory(solver_name)

    for key, value in solver_options.items():
        solver.options[key] = value

    try:
        print("\n=== Starting Optimization ===")
        results = solver.solve(model, tee=True)
        
        print("\n=== Optimization Results ===")
        print(f"Solver Status: {results.solver.status}")
        print(f"Termination Condition: {results.solver.termination_condition}")
        
        if (results.solver.termination_condition == TerminationCondition.optimal or
            results.solver.termination_condition == TerminationCondition.locallyOptimal):
            
            # Extract results
            reconciled_fluxes = {
                reaction_names[j]: pyo.value(model.v[j])
                for j in model.reactions
            }
            
            objective_value = pyo.value(model.obj)
            
            results_dict = {
                'reconciled_reaction_fluxes': reconciled_fluxes,
                'objective_value': objective_value,
                'solver_status': 'optimal',
                'original_measurements': measurements
            }
            
            print(f"\nObjective value: {objective_value}")
            
        else:
            print(f"\nSolver terminated with condition: {results.solver.termination_condition}")
            results_dict = {
                'solver_status': 'failed', 
                'termination_condition': str(results.solver.termination_condition)
            }
            
    except Exception as e:
        print(f"\nError during solving: {str(e)}")
        results_dict = {'solver_status': 'error', 'error_message': str(e)}
    
    return results_dict


def flux_reconciliation(model_path: str,
                       measured_fluxes_df: pd.DataFrame,
                       fixed_fluxes_df: pd.DataFrame,
                       obj_rxn: str,
                       solver_name: Optional[str] = 'ipopt',
                       solver_options: Optional[Dict] = None) -> pd.DataFrame:
    """
    Main function to perform flux reconciliation.
    
    Parameters:
    -----------
    model_path : str
        Path to the genome-scale model file
    measured_fluxes_df : pd.DataFrame
        DataFrame with columns: ['rxn_id', 'exp_flux']
    fixed_fluxes_df : pd.DataFrame
        DataFrame with columns: ['rxn_id', 'exp_fixed_flux']
    obj_rxn : str
        Reaction ID of the objective function
    solver_name : str
        Name of the solver to use
    solver_options : Dict
        Dictionary of solver options
        
    Returns:
    --------
    pd.DataFrame
        Results dataframe with reconciled fluxes and original measurements
    """
    # Load model
    model = load_gem_model(model_path)
    
    # Extract model information
    S, metabolite_names, reaction_names, n_metabolites, n_reactions, flux_bounds = extract_model_info(model)
    
    # Validate stoichiometric matrix
    print("\n=== Validating Stoichiometric Matrix ===")
    zero_rows = np.where(~S.any(axis=1))[0]
    if len(zero_rows) > 0:
        print(f"Warning: Found {len(zero_rows)} metabolites with no reactions:")
        for idx in zero_rows[:5]:  # Show first 5
            print(f"  {metabolite_names[idx]}")
            
    zero_cols = np.where(~S.any(axis=0))[0]
    if len(zero_cols) > 0:
        print(f"Warning: Found {len(zero_cols)} reactions with no metabolites:")
        for idx in zero_cols[:5]:  # Show first 5
            print(f"  {reaction_names[idx]}")
    
    # Process input dataframes
    measurements = process_measured_fluxes_df(measured_fluxes_df, reaction_names)
    fixed_fluxes = process_fixed_fluxes_df(fixed_fluxes_df, reaction_names, measurements)
    
    # Validate measurements against bounds
    print("\n=== Validating Measurements Against Bounds ===")
    bound_violations = []
    for rxn_id, measured_value in measurements.items():
        if rxn_id in flux_bounds:
            lb, ub = flux_bounds[rxn_id]
            if measured_value < lb - 1e-6 or measured_value > ub + 1e-6:
                bound_violations.append((rxn_id, measured_value, lb, ub))
                
    if bound_violations:
        print("\nWarning: Found measurements violating bounds:")
        for rxn, val, lb, ub in bound_violations:
            print(f"  {rxn}: value={val}, bounds=[{lb}, {ub}]")
    
    # Build optimization model
    optimization_model = build_flux_reconciliation_model(
        S, obj_rxn, metabolite_names, reaction_names, measurements,
        fixed_fluxes, flux_bounds
    )
    
    # Solve model
    results = solve_flux_reconciliation_model(
        optimization_model, reaction_names, measurements,
        solver_name, solver_options
    )
    
    # Create results dataframe
    if results['solver_status'] == 'optimal':
        results_data = []
        
        for rxn_id in reaction_names:
            row_data = {
                'rxn_id': rxn_id,
                'reconciled_flux': results['reconciled_reaction_fluxes'][rxn_id],
                'original_measured_flux': measurements.get(rxn_id, None)
            }
            results_data.append(row_data)
        
        results_df = pd.DataFrame(results_data)
        return results_df
    
    else:
        raise RuntimeError(f"Optimization failed: {results.get('error_message', 'Unknown error')}")
