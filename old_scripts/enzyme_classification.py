import cobra
import pandas as pd
import os
import re

'''
Script for: classifying enzymes from a GEM using GPR rules.

'''


def load_model(model_id):
    """
    Load a metabolic model using COBRA.
    
    Parameters:
    -----------
    model_id : str
        ID of the model (e.g., 'iML1515') to load using COBRA's model repository,
        or path to the model file if not using the cobra repository.
        
    Returns:
    --------
    cobra.Model or None
        The loaded COBRA model, or None if loading failed.
    """
    try:
        # From the cobra repository
        model = cobra.io.load_model(model_id)
        print(f"Model {model_id} loaded successfully from COBRA repository.")
        return model
    except Exception:
        # Try alternate loading method (model ID is a file path)
        print(f"Loading model as path...")
        try:
            model = cobra.io.read_sbml_model(model_id)
            print(f"Model {model_id} loaded successfully from path.")
            return model
        except Exception as e2:
            print(f"Failed to load model: {e2}")
            return None


def count_enzyme_types(model): 
    """     
    Returns an enzyme dictionary with each gene's information and a summary count.
    
    Classification types:
      • 1: homomeric enzymes    (single gene or geneX AND geneX) 
      • 2: isoenzymes           (gene1 OR gene2 …) - also called promiscuous in process_gpr
      • 3: heteromeric complexes (geneA AND geneB …, ≥2 different genes) 
    
    Parameters 
    ---------- 
    model : cobra.Model  (already loaded with cobrapy)
    uniprot_mapping : dict, optional
        A dictionary mapping model gene IDs to UniProt information including:
        - uniprot_id
        - sequence
        - ec_number
    
    Returns 
    ------- 
    enzyme_dict : dict
        A dictionary with gene IDs as keys, containing enzyme information:
        - uniprot_id: UniProt identifier (if mapping provided)
        - protein_seq: Protein sequence (if mapping provided)
        - ec_number: EC number (if mapping provided)
        - reactions: List of reaction IDs catalyzed by this enzyme
        - type: Enzyme type (1=homomeric, 2=promiscuous/isoenzyme, 3=homomeric complex, 4=heteromeric complex)
        - components: List of gene IDs for complex subunits (only for protein complexes)
    
    count_summary : dict
        A dictionary with counts of each enzyme type
    """ 
    import re
    from sympy import symbols, sympify
    from sympy.logic.boolalg import to_dnf
    
    # Initialize the output dictionaries
    enzyme_dict = {}
    counts = {1: 0, 2: 0, 3: 0, 4: 0, 'error': 0}  # Use numeric keys to match original
    
    # Default empty mapping if none provided
    if uniprot_mapping is None:
        uniprot_mapping = {}
    
    # Process each reaction in the model
    for rxn in model.reactions:
        gpr = rxn.gene_reaction_rule.strip()
        if not gpr:  # Skip reactions without a GPR
            continue
            
        # Get genes involved in this reaction
        reaction_genes = [gene.id for gene in rxn.genes]
        if not reaction_genes:
            continue
            
        # --- 1. First normalize all logic operators --------------------------- 
        # Convert all operator variations to & and | with proper spacing
        expr = gpr.replace(" and ", " & ").replace(" or ", " | ") 
        expr = re.sub(r'\b(?:and|AND)\b', ' & ', expr) 
        expr = re.sub(r'\b(?:or|OR)\b', ' | ', expr) 
        expr = re.sub(r'\s+', ' ', expr).strip()  # Clean whitespace
        
        # --- 2. Find actual gene IDs (excluding operators) ------------------- 
        genes = []
        for token in re.findall(r'[A-Za-z0-9_.-]+', expr):
            # Skip tokens that are actually operators
            if token.lower() not in ('and', 'or', '&', '|'):
                genes.append(token)
        
        genes = sorted(set(genes))
        if not genes: 
            continue 
            
        # Create symbol mapping for actual genes only - use underscores to avoid SymPy conflicts
        sym_map = {g: f"gene_{idx}" for idx, g in enumerate(genes)} 
        
        # --- 3. Replace gene IDs with symbols in the expression --------------
        for gene, sym in sym_map.items(): 
            expr = re.sub(rf'\b{re.escape(gene)}\b', sym, expr) 
        
        # --- 4. Clean up and normalize the expression -----------------------
        # Ensure all operators have spaces on both sides
        expr = re.sub(r'(\s*)&(\s*)', r' & ', expr)
        expr = re.sub(r'(\s*)\|(\s*)', r' | ', expr)
        
        # Fix implicit AND operations
        expr = re.sub(r'(\bgene_\d+\b)\s+(\bgene_\d+\b)', r'\1 & \2', expr)  # gene_0 gene_1 -> gene_0 & gene_1
        expr = re.sub(r'(\bgene_\d+\b)\s*\(', r'\1 & (', expr)          # gene_0(... -> gene_0 & (...
        expr = re.sub(r'\)\s*(\bgene_\d+\b)', r') & \1', expr)          # )gene_1 -> ) & gene_1
        expr = re.sub(r'\)\s*\(', r') & (', expr)                   # )( -> ) & (
        
        # Final spacing cleanup
        expr = re.sub(r'\s+', ' ', expr).strip()
        
        # --- 5. Process and classify the Boolean expression ----------------
        enzyme_type = 0  # Default, will be updated
        try:
            locals_dict = {v: symbols(v) for v in sym_map.values()} 
            
            # Important: set simplify=False to prevent over-simplification
            sympy_expr = sympify(expr, locals=locals_dict)
            dnf = to_dnf(sympy_expr, simplify=False)
            
            # Classify the expression
            is_isoenzyme = False
            if "or" in gpr.lower() and "(" in gpr and ")" in gpr:
                # Isoenzymes with one option being a complex - e.g. "b3857 or (b3857 and b3856)
                enzyme_type = 2 
                is_isoenzyme = True
            elif dnf.func.__name__ == "Or":  
                # Standard for isoenzymes
                enzyme_type = 2 
                is_isoenzyme = True
            elif dnf.func.__name__ == "Symbol":  
                # Single gene
                enzyme_type = 1 
            elif dnf.func.__name__ == "And":
                # Complex
                genes_in_term = {str(arg) for arg in dnf.args} 
                if len(genes_in_term) == 1:  
                    # geneX & geneX → homomeric complex
                    enzyme_type = 3
                else:  
                    # geneA & geneB … → heteromeric complex
                    enzyme_type = 4 
            else:
                # For any other case, examine the original GPR for 'or'
                if " or " in gpr.lower() or " | " in expr:
                    enzyme_type = 2  # Isoenzyme
                    is_isoenzyme = True
                else:
                    enzyme_type = 1  # Default to homomeric enzyme
            
            # Update enzyme dictionary for each gene based on its type
            if is_isoenzyme:
                # Each gene is an independent enzyme in an OR relationship
                for gene_id in reaction_genes:
                    # Retrieve UniProt info if available
                    uniprot_id = ""
                    protein_seq = ""
                    ec_number = ""
                    if gene_id in uniprot_mapping:
                        uniprot_info = uniprot_mapping[gene_id]
                        uniprot_id = uniprot_info.get('uniprot_id', "")
                        protein_seq = uniprot_info.get('sequence', "")
                        ec_number = uniprot_info.get('ec_number', "")
                    
                    # Initialize or update enzyme entry
                    if gene_id not in enzyme_dict:
                        enzyme_dict[gene_id] = {
                            'uniprot_id': uniprot_id,
                            'protein_seq': protein_seq,
                            'ec_number': ec_number,
                            'reactions': [rxn.id],
                            'type': enzyme_type
                        }
                    else:
                        # Update type if needed (e.g., if gene already exists as another type)
                        enzyme_dict[gene_id]['type'] = max(enzyme_dict[gene_id]['type'], enzyme_type)
                        if rxn.id not in enzyme_dict[gene_id]['reactions']:
                            enzyme_dict[gene_id]['reactions'].append(rxn.id)
            
            elif enzyme_type in [3, 4]:  # Complex (homomeric or heteromeric)
                # For each gene in the complex
                for gene_id in reaction_genes:
                    # Retrieve UniProt info if available
                    uniprot_id = ""
                    protein_seq = ""
                    ec_number = ""
                    if gene_id in uniprot_mapping:
                        uniprot_info = uniprot_mapping[gene_id]
                        uniprot_id = uniprot_info.get('uniprot_id', "")
                        protein_seq = uniprot_info.get('sequence', "")
                        ec_number = uniprot_info.get('ec_number', "")
                    
                    # Initialize enzyme entry if it doesn't exist
                    if gene_id not in enzyme_dict:
                        enzyme_dict[gene_id] = {
                            'uniprot_id': uniprot_id,
                            'protein_seq': protein_seq,
                            'ec_number': ec_number,
                            'reactions': [rxn.id],
                            'type': enzyme_type,
                            'components': []
                        }
                    else:
                        # Update existing entry
                        if 'components' not in enzyme_dict[gene_id]:
                            enzyme_dict[gene_id]['components'] = []
                        # Update type if needed
                        enzyme_dict[gene_id]['type'] = max(enzyme_dict[gene_id]['type'], enzyme_type)
                        if rxn.id not in enzyme_dict[gene_id]['reactions']:
                            enzyme_dict[gene_id]['reactions'].append(rxn.id)
                    
                    # Add other genes as components
                    for other_gene_id in reaction_genes:
                        if other_gene_id != gene_id and other_gene_id not in enzyme_dict[gene_id]['components']:
                            enzyme_dict[gene_id]['components'].append(other_gene_id)
            
            else:  # Homomeric enzyme
                for gene_id in reaction_genes:
                    # Retrieve UniProt info if available
                    uniprot_id = ""
                    protein_seq = ""
                    ec_number = ""
                    if gene_id in uniprot_mapping:
                        uniprot_info = uniprot_mapping[gene_id]
                        uniprot_id = uniprot_info.get('uniprot_id', "")
                        protein_seq = uniprot_info.get('sequence', "")
                        ec_number = uniprot_info.get('ec_number', "")
                    
                    # Initialize or update enzyme entry
                    if gene_id not in enzyme_dict:
                        enzyme_dict[gene_id] = {
                            'uniprot_id': uniprot_id,
                            'protein_seq': protein_seq,
                            'ec_number': ec_number,
                            'reactions': [rxn.id],
                            'type': enzyme_type
                        }
                    else:
                        # Update type if needed
                        enzyme_dict[gene_id]['type'] = max(enzyme_dict[gene_id]['type'], enzyme_type)
                        if rxn.id not in enzyme_dict[gene_id]['reactions']:
                            enzyme_dict[gene_id]['reactions'].append(rxn.id)
            
            # Update the count for this enzyme type
            if enzyme_type > 0:
                counts[enzyme_type] += 1
                
        except Exception as e:
            counts["error"] += 1
            print(f"Error parsing expression: '{expr}' from original: '{gpr}'")
            print(f"Genes: {genes}")
            print(f"Symbol map: {sym_map}")
            print(f"Error: {str(e)}")
    
    return enzyme_dict, counts