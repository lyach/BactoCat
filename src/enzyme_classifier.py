"""
Module Description: enzyme_classifier.py

A comprehensive toolkit for analyzing and classifying enzyme types from 
Gene-Protein-Reaction (GPR) rules in genome-scale metabolic models (GEMs).

This module processes COBRA metabolic models to extract and classify 
enzymes based on their genetic architecture as encoded in GPR rules. It automatically 
categorizes enzymes into three main types based on their gene organization patterns.

3 Types of Enzyme Classification:
    - Homomeric enzymes: Single gene products that function alone
    - Enzyme complexes: Multi-subunit enzymes requiring multiple genes (connected by AND logic)
    - Isoenzymes: Alternative enzyme forms that can perform the same reaction (connected by OR logic)
"""


import ast
import pandas as pd
import cobra

from cobra.core.gene import GPR
from itertools import product


def classify_gpr_type(gpr_rule):
    """
    Classify the GPR rule as 'simple', 'or_only', 'and_only', or 'complex'.
    """
    gpr = str(gpr_rule).lower()
    has_and = 'and' in gpr
    has_or = 'or' in gpr
    if not has_and and not has_or:
        return 'simple'
    elif has_or and not has_and:
        return 'or_only'
    elif has_and and not has_or:
        return 'and_only'
    else:
        return 'complex'


def _ast_to_dnf(node):
    """
    Recursively expand a GPR abstract-syntax-tree node into disjunctive normal
    form (DNF), i.e. an OR of AND-groups.

    Parameters
    ----------
    node : ast.AST
        A node from ``cobra.core.gene.GPR.body`` (``ast.Name`` for a single
        gene or ``ast.BoolOp`` for AND/OR combinations).

    Returns
    -------
    list[frozenset[str]]
        Each frozenset is one OR-branch containing the gene IDs joined by AND.
    """
    if node is None:
        return []
    if isinstance(node, ast.Name):
        return [frozenset({node.id})]
    if isinstance(node, ast.BoolOp):
        children = [_ast_to_dnf(value) for value in node.values]
        if isinstance(node.op, ast.Or):
            # OR: concatenate the branches of every child
            branches = []
            for child in children:
                branches.extend(child)
            return branches
        # AND: cartesian product, unioning one branch from each child
        return [frozenset().union(*combo) for combo in product(*children)]
    raise ValueError(f"Unsupported GPR AST node: {type(node).__name__}")


def gpr_to_dnf(gpr):
    """
    Convert a GPR rule into disjunctive normal form (an OR of AND-groups).

    Each AND-group represents one alternative enzyme (a single gene for a
    homomeric/isozyme alternative, or several genes for a complex). The list of
    groups represents the isozyme alternatives for the reaction.

    Parameters
    ----------
    gpr : cobra.Reaction, cobra.core.gene.GPR, or str
        A reaction, a GPR object, or a GPR rule string.

    Returns
    -------
    list[frozenset[str]]
        DNF representation. Empty list when the GPR is empty.

    Examples
    --------
    ``(b1 and b2) or b3`` -> ``[{b1, b2}, {b3}]``
    """
    # Resolve to a GPR object
    if isinstance(gpr, cobra.Reaction):
        gpr_obj = gpr.gpr
    elif isinstance(gpr, GPR):
        gpr_obj = gpr
    elif isinstance(gpr, str):
        rule = gpr.strip()
        if not rule or rule.lower() == 'none':
            return []
        gpr_obj = GPR.from_string(rule)
    else:
        raise TypeError(f"Unsupported GPR input type: {type(gpr).__name__}")

    body = getattr(gpr_obj, 'body', None)
    return _ast_to_dnf(body)


def classify_enzyme_from_dnf(dnf):
    """
    Assign a reaction-level enzyme class from a DNF representation.

    Parameters
    ----------
    dnf : list[frozenset[str]]
        DNF as returned by :func:`gpr_to_dnf`.

    Returns
    -------
    str
        One of ``'homomeric'``, ``'isoenzyme'``, ``'complex'``, ``'mixed'`` or
        ``'none'`` (empty GPR).

        - homomeric : a single gene catalyses the reaction.
        - complex   : a single AND-group of multiple genes (one heteromer).
        - isoenzyme : several OR alternatives, each a single gene.
        - mixed     : several OR alternatives where at least one is a complex
                      (i.e. isozymes that are themselves heteromers).
    """
    if not dnf:
        return 'none'
    n_branches = len(dnf)
    max_branch_size = max(len(branch) for branch in dnf)
    all_genes = frozenset().union(*dnf)
    if len(all_genes) == 1:
        return 'homomeric'
    if n_branches == 1:
        # single AND-group with >1 gene
        return 'complex'
    if max_branch_size == 1:
        # several single-gene alternatives
        return 'isoenzyme'
    return 'mixed'


def build_enzyme_id(rxn_id, all_genes):
    """
    Build a stable identifier for the integrated enzyme of a reaction.

    For all enzyme classes the integrated enzyme is keyed by the reaction and
    its sorted gene set, so every gene row of a given reaction shares one id.
    """
    genes_token = '_'.join(sorted(str(g) for g in all_genes))
    return f"{rxn_id}__{genes_token}"


def create_gpr_dataframe(model):
    """
    Create a dataframe that represents gene-protein-reaction rules from a COBRA model
    
    Parameters:
        model (cobra.Model): COBRA genome-scale metabolic model
    
    Returns:
        pd.DataFrame: DataFrame with columns
        ['gene', 'type', 'rxn', 'subsystem', 'subunit', 'GPR', 'enzyme_ID',
         'gpr_class', 'enzyme_class', 'enzyme_id', 'n_components']

        - ``gpr_class`` : raw GPR structure ('simple', 'or_only', 'and_only',
          'complex').
        - ``enzyme_class`` : reaction-level biological class ('homomeric',
          'isoenzyme', 'complex', 'mixed') derived from the GPR DNF.
        - ``enzyme_id`` : stable identifier shared by all gene rows of the same
          integrated enzyme (reaction + sorted gene set).
        - ``n_components`` : number of distinct genes in the integrated enzyme.
    """
    
    rows = []
    
    # Iterate through all reactions in the model
    for reaction in model.reactions:
        gpr_rule = str(reaction.gene_reaction_rule)
        
        # Skip reactions without GPR rules
        if not gpr_rule or gpr_rule.strip() == '' or gpr_rule == 'None':
            continue
        
        # Get genes involved in this reaction
        genes_in_rule = [gene.id for gene in reaction.genes]
        
        if not genes_in_rule:
            continue
        
        # Classify GPR type for this reaction
        gpr_class = classify_gpr_type(gpr_rule)

        # Robustly decompose the GPR into isozyme alternatives (OR of AND-groups)
        dnf = gpr_to_dnf(reaction.gpr)
        enzyme_class = classify_enzyme_from_dnf(dnf)
        all_genes = sorted(frozenset().union(*dnf)) if dnf else sorted(genes_in_rule)
        enzyme_id = build_enzyme_id(reaction.id, all_genes)
        
        # Process each gene in the GPR rule
        for gene in genes_in_rule:
            # Determine per-gene enzyme type and the subunits it is bound to
            enzyme_type, subunits, legacy_id = determine_enzyme_properties(
                gene, dnf, reaction.id, all_genes
            )
            row = {
                'gene': gene,
                'type': enzyme_type,
                'rxn': reaction.id,
                'subsystem': reaction.subsystem,
                'subunit': subunits if subunits else '-',
                'GPR': gpr_rule,
                'enzyme_ID': legacy_id,
                'gpr_class': gpr_class,
                'enzyme_class': enzyme_class,
                'enzyme_id': enzyme_id,
                'n_components': len(all_genes),
            }
            rows.append(row)
    
    return pd.DataFrame(rows)


def determine_enzyme_properties(gene, dnf, rxn_id, all_genes):
    """
    Determine the per-gene enzyme type, its complex partners, and a legacy
    enzyme ID, using the DNF decomposition of the GPR.

    Parameters
    ----------
    gene : str
        Current gene ID.
    dnf : list[frozenset[str]]
        DNF of the reaction GPR (OR of AND-groups), from :func:`gpr_to_dnf`.
    rxn_id : str
        Reaction ID.
    all_genes : list[str]
        All genes in the GPR rule.

    Returns
    -------
    tuple(str, str | None, str)
        (enzyme_type, subunits, legacy_enzyme_id)
    """
    # Single gene reaction (homomeric enzyme)
    if len(all_genes) == 1:
        return 'homomeric', None, f"{gene}_h_{rxn_id}"

    # Branches of the DNF that contain this gene
    gene_branches = [branch for branch in dnf if gene in branch]

    # Genes that this gene is bound to via AND (across every branch it appears in)
    complex_partners = set()
    in_complex = False
    for branch in gene_branches:
        if len(branch) > 1:
            in_complex = True
            complex_partners.update(branch)

    if in_complex:
        subunits = ','.join(sorted(complex_partners))
        complex_genes = ''.join(sorted(complex_partners))
        return 'complex', subunits, f"{complex_genes}_c_{rxn_id}"

    # Gene only ever appears as a stand-alone OR alternative -> isoenzyme
    return 'isoenzyme', None, f"{gene}_i_{rxn_id}"


def analyze_model_gprs(model):
    """
    Analyze GPR rules in the model and provide summary statistics
    
    Parameters:
        model (cobra.Model): COBRA genome-scale metabolic model
    
    Returns:
        dict: Summary statistics
    """
    
    total_reactions = len(model.reactions)
    reactions_with_gpr = sum(1 for r in model.reactions if r.gene_reaction_rule and str(r.gene_reaction_rule) != '')
    total_genes = len(model.genes)
    
    gpr_types = {'simple': 0, 'or_only': 0, 'and_only': 0, 'complex': 0}
    
    for reaction in model.reactions:
        gpr = str(reaction.gene_reaction_rule).lower()
        if not gpr or gpr == 'none':
            continue
            
        has_and = 'and' in gpr
        has_or = 'or' in gpr
        
        if not has_and and not has_or:
            gpr_types['simple'] += 1
        elif has_or and not has_and:
            gpr_types['or_only'] += 1
        elif has_and and not has_or:
            gpr_types['and_only'] += 1
        else:
            gpr_types['complex'] += 1
    
    return {
        'total_reactions': total_reactions,
        'reactions_with_gpr': reactions_with_gpr,
        'total_genes': total_genes,
        'gpr_complexity': gpr_types
    }


def create_gem_rxns_df(model):
    '''
    Extract GEM reaction information from a model.
    
    Parameters:
        model (cobra.Model): COBRA genome-scale metabolic model
    
    Returns:
        pd.DataFrame: DataFrame with columns 
        ['gem_rxn_id', 'gem_rxn', 'gem_bigg', 'gem_rxn_name']
    '''
    
    data = []
    for rxn in model.reactions:
        bigg_id = rxn.annotation.get('bigg.reaction', None)
        if isinstance(bigg_id, list):
            bigg_id = ', '.join(bigg_id) # Ensure it's a list for uniformity
        data.append({
            'gem_rxn_id': rxn.id,
            'gem_rxn': rxn.reaction,
            'gem_bigg': bigg_id,
            'gem_rxn_name': rxn.name
        })
    return pd.DataFrame(data)
