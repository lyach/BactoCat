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
            measurements[rxn_id] = row.iloc[1]
        else:
            missing_reactions.append(rxn_id)
            warnings.warn(f"Reaction {rxn_id} not found in model reactions")
    
    # Debug prints
    print("\n=== Measured Fluxes Information ===")
    print(f"Total measured reactions: {len(measured_fluxes_df)}")
    print(f"Matched measurements: {len(measurements)}")
    print(f"Missing reactions: {len(missing_reactions)}")
    if missing_reactions:
        print("Missing reactions:", missing_reactions)
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
                print("Using fixed value as the new measured value.")
                
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
                                  flux_bounds: Dict,
                                  barrier_param: float = 1e-6) -> pyo.ConcreteModel:
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
    barrier_param : float
        Barrier parameter for complementarity constraints
        
    Returns:
    --------
    pyo.ConcreteModel
        The constructed Pyomo model
    """
    n_metabolites = len(metabolite_names)
    n_reactions = len(reaction_names)
    
    # Find biomass reaction index
    if obj_rxn not in reaction_names:
        raise ValueError(f"Objective reaction {obj_rxn} not found in model")
    biomass_idx = reaction_names.index(obj_rxn)
    
    # Set flux bounds
    bounds = {}
    default_lower = -1000.0
    default_upper = 1000.0

    for i, reaction_name in enumerate(reaction_names):
        if reaction_name in flux_bounds:
            bounds[i] = flux_bounds[reaction_name]
        else:
            bounds[i] = (default_lower, default_upper)
    
    # Initialize Pyomo model
    model = pyo.ConcreteModel()
            
    # Sets
    model.metabolites = pyo.Set(initialize=range(n_metabolites))
    model.reactions = pyo.Set(initialize=range(n_reactions))
    model.measured_reactions = pyo.Set(initialize=[
        i for i, name in enumerate(reaction_names) if name in measurements
    ])

    # Variables
    # Flux variables (decision variables from outer problem)
    model.v = pyo.Var(model.reactions, bounds=(-1000, 1000))
    
    # Dual variables (inner problem) 
    model.lambda_mass = pyo.Var(model.metabolites, bounds=(-1000, 1000))  # Mass balance duals
    model.alpha_L = pyo.Var(model.reactions, bounds=(0, 1000))  # Lower bound duals
    model.alpha_U = pyo.Var(model.reactions, bounds=(0, 1000))  # Upper bound duals
    
    # Apply flux bounds to reaction variables
    for i in model.reactions:
        model.v[i].setlb(bounds[i][0])
        model.v[i].setub(bounds[i][1])
    
    # Objective function: minimize weighted sum of squared deviations
    def objective_rule(m):
        return sum(
            ((m.v[i] - measurements[reaction_names[i]]) / measurements[reaction_names[i]])**2  
            for i in m.measured_reactions
        )
    model.obj = pyo.Objective(rule=objective_rule, sense=pyo.minimize)
    
    # Constraints from inner problem optimality conditions
    
    # 1. Mass balance constraints (primal feasibility)
    def mass_balance_rule(m, i):
        return sum(S[i, j] * m.v[j] for j in m.reactions) == 0
    model.mass_balance = pyo.Constraint(model.metabolites, rule=mass_balance_rule)
    
    # 2. Dual feasibility (stationarity condition)
    # For biomass reaction: 1 + S^T * lambda + alpha_L - alpha_U = 0
    def dual_biomass_rule(m):
        return (1 + sum(S[i, biomass_idx] * m.lambda_mass[i] for i in m.metabolites) + 
                m.alpha_L[biomass_idx] - m.alpha_U[biomass_idx] == 0)
    model.dual_biomass = pyo.Constraint(rule=dual_biomass_rule)
    
    # For all other reactions: S^T * lambda + alpha_L - alpha_U = 0
    def dual_other_rule(m, j):
        if j == biomass_idx:
            return pyo.Constraint.Skip
        return (sum(S[i, j] * m.lambda_mass[i] for i in m.metabolites) + 
                m.alpha_L[j] - m.alpha_U[j] == 0)
    model.dual_other = pyo.Constraint(model.reactions, rule=dual_other_rule)
    
    # 3. Complementarity constraints (relaxed with barrier parameter)
    def comp_lower_rule(m, j):
        return (m.v[j] - bounds[j][0]) * m.alpha_L[j] <= barrier_param
    model.comp_lower = pyo.Constraint(model.reactions, rule=comp_lower_rule)
    
    def comp_upper_rule(m, j):
        return (bounds[j][1] - m.v[j]) * m.alpha_U[j] <= barrier_param
    model.comp_upper = pyo.Constraint(model.reactions, rule=comp_upper_rule)
        
    print(f"Number of mass balance constraints: {len(model.mass_balance)}")
    print(f"Number of dual feasibility constraints: {len(model.dual_other) + 1}")
    print(f"Number of complementarity constraints: {len(model.comp_lower) + len(model.comp_upper)}")
    print(f"Barrier parameter: {barrier_param}")
    
    return model


def solve_flux_reconciliation_model(model: pyo.ConcreteModel,
                                  reaction_names: List[str],
                                  measurements: Dict,
                                  solver_name: Optional[str] = 'ipopt',
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
            
            # Extract dual variables for analysis
            dual_mass_balance = {
                j: pyo.value(model.lambda_mass[j])
                for j in model.metabolites
            }
            
            dual_lower = {
                reaction_names[j]: pyo.value(model.alpha_L[j])
                for j in model.reactions
            }
            
            dual_upper = {
                reaction_names[j]: pyo.value(model.alpha_U[j])
                for j in model.reactions
            }
            
            objective_value = pyo.value(model.obj)
            
            results_dict = {
                'reconciled_reaction_fluxes': reconciled_fluxes,
                'dual_mass_balance': dual_mass_balance,
                'dual_lower_bounds': dual_lower,
                'dual_upper_bounds': dual_upper,
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
                       obj_rxn: str,
                       solver_name: Optional[str] = 'ipopt',
                       solver_options: Optional[Dict] = None,
                       barrier_param: Optional[float] = 1e-6) -> pd.DataFrame:
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
        flux_bounds, barrier_param
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
                'original_measured_flux': measurements.get(rxn_id, None),
                'dual_lower_bound': results['dual_lower_bounds'][rxn_id],
                'dual_upper_bound': results['dual_upper_bounds'][rxn_id]
            }
            results_data.append(row_data)
        
        results_df = pd.DataFrame(results_data)
        
        return results_df
    
    else:
        raise RuntimeError(f"Optimization failed: {results.get('error_message', 'Unknown error')}")



def build_simple_flux_reconciliation(S: np.ndarray, 
                                   metabolite_names: List[str],
                                   reaction_names: List[str],
                                   measurements: Dict,
                                   flux_bounds: Dict,
                                   biomass_reaction: str,
                                   measurement_weights: Optional[Dict] = None,
                                   biomass_weight: float = 1.0) -> pyo.ConcreteModel:
    """
    Build a simplified single-level flux reconciliation model.
    
    This formulation minimizes:
    1. Weighted deviations from measured fluxes
    2. Negative biomass production (to encourage growth)
    
    Subject to:
    - Mass balance constraints
    - Flux bounds
    """
    
    n_metabolites = len(metabolite_names)
    n_reactions = len(reaction_names)
    
    # Find biomass reaction index
    if biomass_reaction not in reaction_names:
        raise ValueError(f"Biomass reaction {biomass_reaction} not found in model")
    biomass_idx = reaction_names.index(biomass_reaction)
    
    # Default weights
    if measurement_weights is None:
        measurement_weights = {rxn_id: 1.0 for rxn_id in measurements.keys()}
    
    # Initialize model
    model = pyo.ConcreteModel()
    
    # Sets
    model.metabolites = pyo.Set(initialize=range(n_metabolites))
    model.reactions = pyo.Set(initialize=range(n_reactions))
    model.measured_reactions = pyo.Set(initialize=[
        i for i, name in enumerate(reaction_names) if name in measurements
    ])
    
    # Variables
    model.v = pyo.Var(model.reactions, bounds=(-1000, 1000))  # Flux variables
    model.dev_pos = pyo.Var(model.measured_reactions, bounds=(0, None))  # Positive deviations
    model.dev_neg = pyo.Var(model.measured_reactions, bounds=(0, None))  # Negative deviations
    
    # Apply flux bounds
    for i in model.reactions:
        rxn_name = reaction_names[i]
        if rxn_name in flux_bounds:
            lb, ub = flux_bounds[rxn_name]
            model.v[i].setlb(lb)
            model.v[i].setub(ub)
    
    # Objective: minimize weighted deviations + encourage biomass
    def objective_rule(m):
        deviation_term = sum(
            measurement_weights.get(reaction_names[i], 1.0) * (m.dev_pos[i] + m.dev_neg[i])
            for i in m.measured_reactions
        )
        # Encourage biomass production (negative coefficient)
        biomass_term = -biomass_weight * m.v[biomass_idx]
        
        return deviation_term + biomass_term
    
    model.obj = pyo.Objective(rule=objective_rule, sense=pyo.minimize)
    
    # Mass balance constraints
    def mass_balance_rule(m, i):
        return sum(S[i, j] * m.v[j] for j in m.reactions) == 0
    
    model.mass_balance = pyo.Constraint(model.metabolites, rule=mass_balance_rule)
    
    # Deviation constraints
    def deviation_rule(m, i):
        rxn_name = reaction_names[i]
        measured_value = measurements[rxn_name]
        return m.v[i] - measured_value == m.dev_pos[i] - m.dev_neg[i]
    
    model.deviation_constraint = pyo.Constraint(model.measured_reactions, rule=deviation_rule)
    
    # Minimum biomass constraint 
    model.min_biomass = pyo.Constraint(expr=model.v[biomass_idx] >= 0.01)
    
    print(f"Model built with {len(model.mass_balance)} mass balance constraints")
    print(f"Number of measured reactions: {len(model.measured_reactions)}")
    print(f"Biomass reaction: {biomass_reaction} (index {biomass_idx})")
    
    return model


def solve_simple_flux_reconciliation(model: pyo.ConcreteModel,
                                   reaction_names: List[str],
                                   measurements: Dict,
                                   solver_name: str = 'ipopt') -> Dict:
    """
    Solve the simplified flux reconciliation model.
    """
    
    solver = pyo.SolverFactory(solver_name)
    
    # Solver options for better numerical stability
    solver.options['tol'] = 1e-6
    solver.options['max_iter'] = 3000
    solver.options['print_level'] = 5
    solver.options['hessian_approximation'] = 'limited-memory'
    
    try:
        print("\n=== Solving Simplified Model ===")
        results = solver.solve(model, tee=True)
        
        print(f"Solver Status: {results.solver.status}")
        print(f"Termination Condition: {results.solver.termination_condition}")
        
        if (results.solver.termination_condition == TerminationCondition.optimal or
            results.solver.termination_condition == TerminationCondition.locallyOptimal):
            
            # Extract results
            reconciled_fluxes = {
                reaction_names[j]: pyo.value(model.v[j])
                for j in model.reactions
            }
            
            # Calculate deviations
            deviations = {}
            for i in model.measured_reactions:
                rxn_name = reaction_names[i]
                measured = measurements[rxn_name]
                reconciled = reconciled_fluxes[rxn_name]
                deviations[rxn_name] = reconciled - measured
            
            objective_value = pyo.value(model.obj)
            
            results_dict = {
                'reconciled_fluxes': reconciled_fluxes,
                'deviations': deviations,
                'objective_value': objective_value,
                'solver_status': 'optimal',
                'measurements': measurements
            }
            
            print(f"\nObjective value: {objective_value:.6f}")
            print("\nMeasurement vs Reconciled:")
            for rxn_name in measurements.keys():
                measured = measurements[rxn_name]
                reconciled = reconciled_fluxes[rxn_name]
                deviation = deviations[rxn_name]
                print(f"  {rxn_name}: measured={measured:.3f}, reconciled={reconciled:.3f}, dev={deviation:.3f}")
            
            return results_dict
            
        else:
            print(f"Optimization failed: {results.solver.termination_condition}")
            return {'solver_status': 'failed', 'termination_condition': str(results.solver.termination_condition)}
            
    except Exception as e:
        print(f"Error during optimization: {e}")
        return {'solver_status': 'error', 'error_message': str(e)}


def simple_flux_reconciliation_main(model_path: str,
                                  measured_fluxes_df: pd.DataFrame,
                                  biomass_reaction: str = 'BIOMASS_Ec_iML1515_core_75p37M',
                                  solver_name: str = 'ipopt') -> pd.DataFrame:
    """
    Main function for simplified flux reconciliation.
    """
    
    # Load model (reuse your existing function)
    model = read_sbml_model(model_path)
    
    # Extract model information (reuse your existing function)
    # You can copy the extract_model_info function from your original code
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
    for reaction in model.reactions:
        flux_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
    
    # Process measurements
    measurements = {}
    for _, row in measured_fluxes_df.iterrows():
        rxn_id = row['rxn_id']
        if rxn_id in reaction_names:
            measurements[rxn_id] = row.iloc[1]
    
    # Build and solve model
    opt_model = build_simple_flux_reconciliation(
        S, metabolite_names, reaction_names, measurements,
        flux_bounds, biomass_reaction
    )
    
    results = solve_simple_flux_reconciliation(opt_model, reaction_names, measurements, solver_name)
    
    # Create results DataFrame
    if results['solver_status'] == 'optimal':
        results_data = []
        for rxn_id in reaction_names:
            results_data.append({
                'rxn_id': rxn_id,
                'reconciled_flux': results['reconciled_fluxes'][rxn_id],
                'measured_flux': measurements.get(rxn_id, None),
                'deviation': results['deviations'].get(rxn_id, None)
            })
        
        return pd.DataFrame(results_data)
    else:
        raise RuntimeError(f"Optimization failed: {results.get('error_message', 'Unknown error')}")