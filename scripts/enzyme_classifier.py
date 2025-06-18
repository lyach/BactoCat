"""
Module Description: enzyme_classifier.py

Purpose: A comprehensive toolkit for analyzing and classifying enzyme types from 
Gene-Protein-Reaction (GPR) rules in genome-scale metabolic models (GEMs).

Overview: This module processes COBRA metabolic models to extract and classify 
enzymes based on their genetic architecture as encoded in GPR rules. It automatically 
categorizes enzymes into three main types based on their gene organization patterns.

3 Types of Enzyme Classification:
    - Homomeric enzymes: Single gene products that function alone
    - Enzyme complexes: Multi-subunit enzymes requiring multiple genes (connected by AND logic)
    - Isoenzymes: Alternative enzyme forms that can perform the same reaction (connected by OR logic)
"""


import pandas as pd
import cobra
import re


def create_gpr_dataframe(model):
    """
    Create a dataframe that represents gene-protein-reaction rules from a COBRA model
    
    Parameters:
        model (cobra.Model): COBRA genome-scale metabolic model
    
    Returns:
        pd.DataFrame: DataFrame with columns 
        ['gene', 'type', 'rxn', 'subunit', 'GPR', 'enzyme_ID']
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
            
        # Process each gene in the GPR rule
        for gene in genes_in_rule:
            # Determine enzyme type based on GPR complexity
            enzyme_type, subunits, enzyme_id = determine_enzyme_properties(
                gene, gpr_rule, reaction.id, genes_in_rule
            )
            
            row = {
                'gene': gene,
                'type': enzyme_type,
                'rxn': reaction.id,
                'subunit': subunits if subunits else '-',
                'GPR': gpr_rule,
                'enzyme_ID': enzyme_id
            }
            rows.append(row)
    
    return pd.DataFrame(rows)


def determine_enzyme_properties(gene, gpr_rule, rxn_id, all_genes):
    """
    Determine enzyme type, subunits, and enzyme ID based on GPR rule
    
    Parameters:
        gene (str): Current gene ID
        gpr_rule (str): GPR rule string
        rxn_id (str): Reaction ID
        all_genes (list): All genes in the GPR rule
    """
    
    # Clean the GPR rule for analysis
    gpr_clean = gpr_rule.replace(' ', '').lower()
    gene_clean = gene.lower()
    
    # Single gene reaction (homomeric enzyme)
    if len(all_genes) == 1:
        return 'homomeric', None, f"{gene}_h_{rxn_id}"
    
    # Multiple genes - check for AND/OR relationships
    has_and = 'and' in gpr_clean
    has_or = 'or' in gpr_clean
    
    if has_and and has_or:
        # Complex case with both AND and OR
        # Check if this gene is part of a complex (connected by AND)
        if is_gene_in_complex(gene, gpr_rule):
            complex_partners = get_complex_partners(gene, gpr_rule, all_genes)
            subunits = ','.join(sorted(complex_partners))
            complex_genes = ''.join(sorted(complex_partners))
            return 'complex', subunits, f"{complex_genes}_c_{rxn_id}"
        else:
            # This gene is an isoenzyme alternative
            return 'isoenzyme', None, f"{gene}_i_{rxn_id}"
    
    elif has_and and not has_or:
        # Pure complex (all genes required)
        subunits = ','.join(sorted(all_genes))
        complex_genes = ''.join(sorted(all_genes))
        return 'complex', subunits, f"{complex_genes}_c_{rxn_id}"
    
    elif has_or and not has_and:
        # Pure isoenzymes (alternative genes)
        return 'isoenzyme', None, f"{gene}_i_{rxn_id}"
    
    else:
        # Fallback case
        return 'isoenzyme', None, f"{gene}_i_{rxn_id}"


def is_gene_in_complex(gene, gpr_rule):
    """
    Check if a gene is part of a protein complex (connected by AND)
    """
    
    # Simple heuristic: if gene appears near an AND, it's likely in a complex
    gpr_parts = re.split(r'\s+or\s+', gpr_rule.lower())
    
    for part in gpr_parts:
        if gene.lower() in part and 'and' in part:
            return True
    return False


def get_complex_partners(gene, gpr_rule, all_genes):
    """
    Get the genes that form a complex with the given gene
    """
    
    # Find the part of GPR rule containing this gene and connected by AND
    gpr_parts = re.split(r'\s+or\s+', gpr_rule.lower())
    gene_lower = gene.lower()
    
    for part in gpr_parts:
        if gene_lower in part and 'and' in part:
            # Extract genes from this part
            genes_in_part = []
            for g in all_genes:
                if g.lower() in part:
                    genes_in_part.append(g)
            if len(genes_in_part) > 1:
                return genes_in_part
    
    # Fallback: return all genes if complex structure unclear
    return all_genes


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
