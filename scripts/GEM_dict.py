import cobra
import pandas as pd
import os
import re

'''
WIP Script for: Extracting enzyme information from a GEM

Need to correct the enzyme type classification.
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



def process_gem(model_id, uniprot_mapping_file):
    """
    ####### WIP #######
    
    Process a GEM to extract enzyme information using a pre-generated UniProt mapping file.
    
    Parameters:
    -----------
    model_id : str
        ID of the model (e.g., 'iML1515') to load using COBRA's model repository.
        Path to the model file if not using the cobra repository.
        
    uniprot_mapping_file : str
        Path to a file mapping model gene IDs to UniProt IDs, EC numbers, and sequences.
        This file is mandatory and should contain columns: model_gene_id, uniprot_id, ec_number, sequence.
        
    Returns:
    --------
    pd.DataFrame
        A DataFrame with model gene IDs as index and enzyme information as columns.
    dict
        A summary dictionary with the count of each enzyme type.
    """    
    # Load the model
    model = load_model(model_id)
    if model is None:
        return None, None
    
    # Initialize the output dictionary
    enzyme_dict = {}
    
    # Load UniProt mapping
    if not os.path.exists(uniprot_mapping_file):
        raise FileNotFoundError(f"UniProt mapping file not found: {uniprot_mapping_file}")
    
    # Load the mapping data
    uniprot_df = pd.read_csv(uniprot_mapping_file)
    
    # Create a dictionary for quick lookups - mapping model gene IDs to UniProt info
    uniprot_mapping = {}
    for _, row in uniprot_df.iterrows():
        model_gene_id = row['model_gene_id']
        uniprot_mapping[model_gene_id] = {
            'uniprot_id': row['uniprot_id'],
            'ec_number': row['ec_number'],
            'sequence': row['sequence']
        }
    
    # Classify enzyme types
    classify_enzyme_types(enzyme_dict)
    
    # Count the number of each type of enzyme
    enzyme_type_summary = {1: 0, 2: 0, 3: 0, 4: 0}
    for info in enzyme_dict.values():
        enzyme_type_summary[info['type']] += 1
    
    # Create DataFrame from enzyme dictionary
    data = []
    for gene_id, info in enzyme_dict.items():
        # Get UniProt info from mapping
        if gene_id in uniprot_mapping:
            uniprot_info = uniprot_mapping[gene_id]
            row = {
                'model_gene_id': gene_id,
                'uniprot_id': uniprot_info['uniprot_id'],
                'ec_number': uniprot_info['ec_number'],
                'protein_seq': uniprot_info['sequence'],
                'num_reactions': len(info['reactions']),
                'enzyme_type': info['type'],
                'reactions': ';'.join(info['reactions']),
            }
            # Add components if it's a complex
            if 'components' in info:
                row['components'] = ';'.join(info['components'])
            else:
                row['components'] = ''
            
            data.append(row)
    
    # Create DataFrame
    df_enzyme_info = pd.DataFrame(data)
    
    # Set model_gene_id as the index
    df_enzyme_info.set_index('model_gene_id', inplace=True)
    
    return df_enzyme_info, enzyme_type_summary


def count_homomeric_enzymes(model):
    """
    Counts the number of homomeric enzymes in a model.
    
    Parameters:
    -----------
    model_id : str
        ID of the model (e.g., 'iML1515') to load using COBRA's model repository,
        or path to the model file if not using the cobra repository.
    
    Returns:
    --------
    enzyme_dict : dict
        A dictionary with enzyme type as keys:
        • homomeric enzymes: list of tuples (reaction_id, gene_id)
        • isoenzymes: list of tuples ([reaction_ids], gene_id)
        • complexes: list of tuples ([reaction_ids], [gene_ids])
    """
    # Initialize the output dictionary
    enzyme_dict = {}
    
    # Process each reaction in the model
    for rxn in model.reactions:
    


def count_enzyme_types(model, uniprot_mapping=None): 
    """ 
    ####### WIP #######
    
    Returns an enzyme dictionary with each gene's information and a summary count.
    Classification types:
      • 1: homomeric enzymes    (single gene or geneX AND geneX …) 
      • 2: isoenzymes           (gene1 OR gene2 …) - also called promiscuous in process_gpr
      • 3: homomeric complexes  (not explicitly classified in original function)
      • 4: heteromeric complexes (geneA AND geneB …, ≥2 different genes) 
    
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



def count_enzyme_types2(model, uniprot_mapping=None): 
    """ 
    ####### WIP #######
    Returns a DataFrame with enzyme information and a summary count.
    Classification types:
      • homomeric: single gene enzymes
      • isoenzyme: multiple alternative genes (OR relationship)
      • homo_complex: homomeric complexes (same gene repeated)
      • hetero_complex: heteromeric complexes (different genes in AND relationship)
    
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
    enzyme_df : pd.DataFrame
        A DataFrame with columns:
        - enzyme_id: unique identifier for each enzyme
        - type: enzyme type (homomeric, isoenzyme, homo_complex, hetero_complex)
        - genes: list of gene IDs associated with the enzyme
        - reactions: list of reaction IDs catalyzed by the enzyme
    
    count_summary : dict
        A dictionary with counts of each enzyme type
    """ 
    import re
    from sympy import symbols, sympify
    from sympy.logic.boolalg import to_dnf
    
    # Initialize the output lists
    enzyme_data = []
    enzyme_id_counter = 0
    counts = {'homomeric': 0, 'isoenzyme': 0, 'homo_complex': 0, 'hetero_complex': 0, 'error': 0}
    
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
            
        # Create symbol mapping for actual genes only
        sym_map = {g: f"g{idx}" for idx, g in enumerate(genes)} 
        
        # --- 3. Replace gene IDs with symbols in the expression --------------
        for gene, sym in sym_map.items(): 
            expr = re.sub(rf'\b{re.escape(gene)}\b', sym, expr) 
        
        # --- 4. Clean up and normalize the expression -----------------------
        # Ensure all operators have spaces on both sides
        expr = re.sub(r'(\s*)&(\s*)', r' & ', expr)
        expr = re.sub(r'(\s*)\|(\s*)', r' | ', expr)
        
        # Fix implicit AND operations
        expr = re.sub(r'(\bg\d+\b)\s+(\bg\d+\b)', r'\1 & \2', expr)  # g0 g1 -> g0 & g1
        expr = re.sub(r'(\bg\d+\b)\s*\(', r'\1 & (', expr)          # g0(... -> g0 & (...
        expr = re.sub(r'\)\s*(\bg\d+\b)', r') & \1', expr)          # )g1 -> ) & g1
        expr = re.sub(r'\)\s*\(', r') & (', expr)                   # )( -> ) & (
        
        # Final spacing cleanup
        expr = re.sub(r'\s+', ' ', expr).strip()
        
        # --- 5. Process and classify the Boolean expression ----------------
        try:
            locals_dict = {v: symbols(v) for v in sym_map.values()} 
            
            # Important: set simplify=False to prevent over-simplification
            sympy_expr = sympify(expr, locals=locals_dict)
            dnf = to_dnf(sympy_expr, simplify=False)
            
            # Classify and create enzyme entries
            if dnf.func.__name__ == "Or":  
                # Isoenzymes - each alternative is a separate enzyme
                enzyme_type = 'isoenzyme'
                
                # Parse OR terms to create individual enzyme entries
                or_terms = dnf.args if hasattr(dnf, 'args') else [dnf]
                for term in or_terms:
                    enzyme_id_counter += 1
                    
                    if term.func.__name__ == "Symbol":
                        # Single gene alternative
                        gene_symbol = str(term)
                        gene_id = [k for k, v in sym_map.items() if v == gene_symbol][0]
                        enzyme_id = gene_id  # For single gene, use gene ID as enzyme ID
                        enzyme_genes = [gene_id]
                    elif term.func.__name__ == "And":
                        # Complex alternative
                        gene_symbols = [str(arg) for arg in term.args]
                        enzyme_genes = [k for k, v in sym_map.items() if v in gene_symbols]
                        enzyme_id = f"complex_{enzyme_id_counter}"
                    else:
                        # Single term that's not a symbol or and
                        gene_symbol = str(term)
                        gene_id = [k for k, v in sym_map.items() if v == gene_symbol][0]
                        enzyme_id = gene_id
                        enzyme_genes = [gene_id]
                    
                    enzyme_data.append({
                        'enzyme_id': enzyme_id,
                        'type': enzyme_type,
                        'genes': enzyme_genes,
                        'reactions': [rxn.id]
                    })
                    counts[enzyme_type] += 1
                    
            elif dnf.func.__name__ == "Symbol":  
                # Single gene enzyme
                enzyme_type = 'homomeric'
                gene_symbol = str(dnf)
                gene_id = [k for k, v in sym_map.items() if v == gene_symbol][0]
                enzyme_id = gene_id  # For single gene, use gene ID as enzyme ID
                
                enzyme_data.append({
                    'enzyme_id': enzyme_id,
                    'type': enzyme_type,
                    'genes': [gene_id],
                    'reactions': [rxn.id]
                })
                counts[enzyme_type] += 1
                
            elif dnf.func.__name__ == "And": 
                # Complex
                gene_symbols = [str(arg) for arg in dnf.args]
                enzyme_genes = [k for k, v in sym_map.items() if v in gene_symbols]
                
                if len(set(enzyme_genes)) == 1:  
                    # Homomeric complex (same gene repeated)
                    enzyme_type = 'homo_complex'
                    enzyme_id = enzyme_genes[0]  # Use the gene ID since it's homomeric
                else:  
                    # Heteromeric complex (different genes)
                    enzyme_type = 'hetero_complex'
                    enzyme_id_counter += 1
                    enzyme_id = f"complex_{enzyme_id_counter}"
                
                enzyme_data.append({
                    'enzyme_id': enzyme_id,
                    'type': enzyme_type,
                    'genes': enzyme_genes,
                    'reactions': [rxn.id]
                })
                counts[enzyme_type] += 1
                
            else:
                # Fallback classification
                if " or " in gpr.lower() or " | " in expr:
                    enzyme_type = 'isoenzyme'
                else:
                    enzyme_type = 'homomeric'
                
                # Create a single enzyme entry for all genes
                if len(reaction_genes) == 1:
                    enzyme_id = reaction_genes[0]
                else:
                    enzyme_id_counter += 1
                    enzyme_id = f"enzyme_{enzyme_id_counter}"
                
                enzyme_data.append({
                    'enzyme_id': enzyme_id,
                    'type': enzyme_type,
                    'genes': reaction_genes,
                    'reactions': [rxn.id]
                })
                counts[enzyme_type] += 1
                
        except Exception as e:
            counts["error"] += 1
            print(f"Error parsing expression: '{expr}' from original: '{gpr}'")
            print(f"Genes: {genes}")
            print(f"Symbol map: {sym_map}")
            print(f"Error: {str(e)}")
            
            # Create fallback entry for error cases
            enzyme_id_counter += 1
            enzyme_data.append({
                'enzyme_id': f"error_enzyme_{enzyme_id_counter}",
                'type': 'error',
                'genes': reaction_genes,
                'reactions': [rxn.id]
            })
    
    # Create DataFrame from enzyme data
    enzyme_df = pd.DataFrame(enzyme_data)
    
    # Consolidate enzymes that appear in multiple reactions
    if not enzyme_df.empty:
        # Group by enzyme_id and consolidate reactions
        consolidated_data = []
        for enzyme_id, group in enzyme_df.groupby('enzyme_id'):
            consolidated_data.append({
                'enzyme_id': enzyme_id,
                'type': group['type'].iloc[0],  # All should be the same
                'genes': group['genes'].iloc[0],  # All should be the same
                'reactions': [rxn for rxn_list in group['reactions'] for rxn in rxn_list]  # Flatten and combine
            })
        
        enzyme_df = pd.DataFrame(consolidated_data)
        
        # Recalculate counts based on consolidated data
        counts = enzyme_df['type'].value_counts().to_dict()
        # Ensure all expected keys are present
        for key in ['homomeric', 'isoenzyme', 'homo_complex', 'hetero_complex', 'error']:
            if key not in counts:
                counts[key] = 0
    
    return enzyme_df, counts


def process_gpr(reaction, enzyme_dict, uniprot_mapping):
    """
    ####### WIP #######
    
    Process gene-protein-reaction (GPR) association for a single reaction in a metabolic model.
    
    Parameters:
    -----------
    reaction : cobra.core.Reaction
        A COBRA model reaction object containing GPR rules and gene associations
    
    enzyme_dict : dict
        A dictionary to store enzyme information, keyed by gene ID. Each entry contains:
        - uniprot_id: UniProt identifier
        - protein_seq: Protein sequence
        - ec_number: EC number
        - reactions: List of reaction IDs catalyzed by this enzyme
        - type: Enzyme type (1=homomeric, 2=promiscuous, 3=homomeric complex, 4=heteromeric complex)
        - components: List of gene IDs for complex subunits (only for protein complexes)
    
    uniprot_mapping : dict
        A dictionary mapping model gene IDs to UniProt information including:
        - uniprot_id
        - sequence
        - ec_number
    """
    # Extract GPR rule
    gpr = reaction.gene_reaction_rule.lower()
    
    # Simple case: single gene (homomeric enzyme)
    if len(reaction.genes) == 1:
        gene = list(reaction.genes)[0]
        gene_id = gene.id
        
        # Get UniProt info from mapping using model_gene_id
        if gene_id in uniprot_mapping:
            uniprot_info = uniprot_mapping[gene_id]
            uniprot_id = uniprot_info['uniprot_id']
            sequence = uniprot_info['sequence']
            ec_number = uniprot_info['ec_number']
        else:
            # If gene is not in mapping, use gene ID as fallback
            uniprot_id = gene_id
            sequence = ""
            ec_number = ""
        
        if gene_id not in enzyme_dict:
            enzyme_dict[gene_id] = {
                'uniprot_id': uniprot_id,
                'protein_seq': sequence,
                'ec_number': ec_number,
                'reactions': [],
                'type': 1  # Will be updated later if needed
            }
        
        # Add reaction information
        if reaction.id not in enzyme_dict[gene_id]['reactions']:
            enzyme_dict[gene_id]['reactions'].append(reaction.id)
    
    # Complex case: multiple genes (potentially heteromeric complex)
    else:
        # Extract genes involved in the reaction
        genes = list(reaction.genes)
        gene_ids = [gene.id for gene in genes]
        
        # Check if it's an OR relationship (isozymes) or AND relationship (complex)
        if ' or ' in gpr:
            # Isozymes case (promiscuous)
            for gene in genes:
                gene_id = gene.id
                
                # Get UniProt info from mapping using model_gene_id
                if gene_id in uniprot_mapping:
                    uniprot_info = uniprot_mapping[gene_id]
                    uniprot_id = uniprot_info['uniprot_id']
                    sequence = uniprot_info['sequence']
                    ec_number = uniprot_info['ec_number']
                else:
                    # If gene is not in mapping, use gene ID as fallback
                    uniprot_id = gene_id
                    sequence = ""
                    ec_number = ""
                
                if gene_id not in enzyme_dict:
                    enzyme_dict[gene_id] = {
                        'uniprot_id': uniprot_id,
                        'protein_seq': sequence,
                        'ec_number': ec_number,
                        'reactions': [],
                        'type': 2  # Promiscuous, will be updated later if needed
                    }
                
                # Add reaction information
                if reaction.id not in enzyme_dict[gene_id]['reactions']:
                    enzyme_dict[gene_id]['reactions'].append(reaction.id)
        
        elif ' and ' in gpr:
            # For each gene in the complex
            for gene in genes:
                gene_id = gene.id
                
                # Initialize the enzyme entry if it doesn't exist
                if gene_id not in enzyme_dict:
                    enzyme_dict[gene_id] = {
                        'uniprot_id': "",  # Complexes don't have a single UniProt ID
                        'protein_seq': "",  # Complexes don't have a single sequence
                        'ec_number': "",    # Complexes don't have a single EC number
                        'reactions': [],
                        'type': 4,  # Heteromeric complex by default, will be updated if homomeric
                        'components': []    # Store component genes
                    }
                # Ensure components list exists even if the entry already exists
                elif 'components' not in enzyme_dict[gene_id]:
                    enzyme_dict[gene_id]['components'] = []
                
                # Add all other genes in the complex as components
                for other_gene_id in gene_ids:
                    if other_gene_id != gene_id and other_gene_id not in enzyme_dict[gene_id]['components']:
                        enzyme_dict[gene_id]['components'].append(other_gene_id)
                
                # Add reaction information
                if reaction.id not in enzyme_dict[gene_id]['reactions']:
                    enzyme_dict[gene_id]['reactions'].append(reaction.id)
        
        else:
            # TO DO: Handle more complex GPR expressions
            # This would require a more sophisticated parser for complex GPR rules
            pass

def classify_enzyme_types(enzyme_dict):
    """
    Classify enzymes into types:
    1: homomeric (single enzyme, single reaction)
    2: promiscuous (single enzyme, multiple reactions)
    3: homomeric complex (complex with identical subunits)
    4: heteromeric complex (complex with different subunits)
    """
    for enzyme_id, info in enzyme_dict.items():
        # For complexes, check if they're homomeric or heteromeric
        if 'components' in info:
            # Check if all components are the same (homomeric complex)
            if len(set(info['components'])) == 1:
                info['type'] = 3  # homomeric complex
            else:
                info['type'] = 4  # heteromeric complex
        else:
            # Single enzyme
            # Check if it catalyzes multiple reactions (promiscuous)
            if len(info['reactions']) > 1:
                info['type'] = 2  # promiscuous
            else:
                info['type'] = 1  # homomeric
                

def collect_substrates(model_id, is_human_gem=False):
    """
    Collect metabolites from a model by iterating through all reactions.
    Returns a DataFrame with rows for each reaction-metabolite pair.
    
    Parameters:
    -----------
    model_id : str
        ID of the model to load using COBRA's model repository.
        Path to the model file if not using the cobra repository.
    is_human_gem : bool, optional
        If True, indicates that the model is the Human GEM and compartment tags
        should not be stripped if they don't have an underscore.
    """
    # Load the model
    model = load_model(model_id)
    if model is None:
        return None
        
    # Create list to store reaction-metabolite pairs
    data = []
    
    for rxn in model.reactions:
        subs = {s.id for s in rxn.reactants}
        if rxn.lower_bound < 0:   # check if it's a reversible reaction
            subs.update(s.id for s in rxn.products)

        # Strip compartment tag based on the model type
        if is_human_gem:
            # For Human GEM, strip the last character if it's a letter
            subs = {re.sub(r'[a-zA-Z]$', '', s) for s in subs}
        else:
            # For other models, strip compartment tag and remove duplicates due to compartments
            subs = {re.sub(r'_[a-z]$', '', s) for s in subs}

        for sub in sorted(subs):
            data.append({'reaction_id': rxn.id, 'metabolite_id': sub})
    
    # Create DataFrame
    reaction_substrates_df = pd.DataFrame(data)
    
    return reaction_substrates_df


def calculate_kcats(df_enzyme_info, reaction_substrates_df):
    """
    Create a DataFrame with expanded reactions and substrates for each model_gene_id.
    
    Parameters:
    -----------
    df_enzyme_info : pd.DataFrame
        DataFrame containing enzyme information with columns 'model_gene_id', 'uniprot_id', and 'reactions'.
        
    reaction_substrates_df : pd.DataFrame
        DataFrame containing reaction-substrate pairs with columns 'reaction_id' and 'metabolite_id'.
        
    Returns:
    --------
    pd.DataFrame
        A DataFrame with columns 'model_gene_id', 'uniprot_id', 'reaction_id', and 'metabolite_id'.
    """
    # List to store expanded data
    expanded_data = []
    
    # Iterate over each row in df_enzyme_info
    for model_gene_id, row in df_enzyme_info.iterrows():
        uniprot_id = row['uniprot_id']
        reactions = row['reactions'].split(';')
        
        # For each reaction, find corresponding substrates
        for reaction_id in reactions:
            substrates = reaction_substrates_df[reaction_substrates_df['reaction_id'] == reaction_id]['metabolite_id']
            for substrate in substrates:
                expanded_data.append({
                    'model_gene_id': model_gene_id,
                    'uniprot_id': uniprot_id,
                    'reaction_id': reaction_id,
                    'metabolite_id': substrate
                })
    
    # Create a new DataFrame from the expanded data
    kcat_df = pd.DataFrame(expanded_data)
    
    return kcat_df

def process_human_gem(model_id, uniprot_mapping_file):
    """
    Process a human GEM to extract enzyme information using a pre-generated UniProt mapping file
    that contains ENSP-ENSG mappings.
    
    Parameters:
    -----------
    model_id : str
        ID of the model (e.g., 'Human1') to load using COBRA's model repository.
        Path to the model file if not using the cobra repository.
        
    uniprot_mapping_file : str
        Path to a file mapping ENSP/ENSG IDs to UniProt IDs, EC numbers, and sequences.
        This file should contain columns: ENSP, ENSG, uniprot_id, ec_number, sequence.
        
    Returns:
    --------
    pd.DataFrame
        A DataFrame with enzyme information as columns.
    dict
        A summary dictionary with the count of each enzyme type.
    """
    
    # Load the model
    model = load_model(model_id)
    if model is None:
        return None, None
    
    # Initialize the output dictionary
    enzyme_dict = {}
    
    # Load UniProt mapping
    if not os.path.exists(uniprot_mapping_file):
        raise FileNotFoundError(f"UniProt mapping file not found: {uniprot_mapping_file}")
    
    # Load the mapping data
    uniprot_df = pd.read_csv(uniprot_mapping_file)
    
    # Create a mapping dictionary, keyed by ENSG (gene id in the model)
    uniprot_mapping = {}
    for _, row in uniprot_df.iterrows():
        # Skip rows with missing ENSG or ENSP
        if pd.isna(row.get('ENSG')) or pd.isna(row.get('ENSP')):
            continue
            
        ensg = row['ENSG']
        ensp = row['ENSP']
        
        # If this ENSG doesn't exist in the mapping yet, create a new entry
        if ensg not in uniprot_mapping:
            uniprot_mapping[ensg] = {
                'ENSP_list': [],
                'uniprot_ids': [],
                'ec_numbers': [],
                'sequences': []
            }
        
        # Add the ENSP and related data to the ENSG entry
        uniprot_mapping[ensg]['ENSP_list'].append(ensp)
        uniprot_mapping[ensg]['uniprot_ids'].append(row.get('uniprot_id', ''))
        uniprot_mapping[ensg]['ec_numbers'].append(row.get('ec_number', ''))
        uniprot_mapping[ensg]['sequences'].append(row.get('sequence', ''))
    
    # Process each reaction in the model
    for reaction in model.reactions:
        # Skip reactions without genes
        if not reaction.genes:
            continue
        
        # Process the gene-protein-reaction association
        process_human_gpr(reaction, enzyme_dict, uniprot_mapping)
    
    # Classify enzyme types
    classify_enzyme_types(enzyme_dict)
    
    # Count the number of each type of enzyme
    enzyme_type_summary = {1: 0, 2: 0, 3: 0, 4: 0}
    for info in enzyme_dict.values():
        enzyme_type_summary[info['type']] += 1
    
    # Create DataFrame from enzyme dictionary
    data = []
    for gene_id, info in enzyme_dict.items():
        # Get information for each ENSP associated with this gene
        ensp_list = uniprot_mapping.get(gene_id, {}).get('ENSP_list', [])
        
        if not ensp_list:
            # If no ENSP mapping found, create a single row with just the gene ID
            row = {
                'ENSG': gene_id,
                'ENSP': '',
                'uniprot_id': info.get('uniprot_id', ''),
                'ec_number': info.get('ec_number', ''),
                'protein_seq': info.get('protein_seq', ''),
                'num_reactions': len(info['reactions']),
                'enzyme_type': info['type'],
                'reactions': ';'.join(info['reactions']),
            }
            # Add components if it's a complex
            if 'components' in info:
                row['components'] = ';'.join(info['components'])
            else:
                row['components'] = ''
            
            data.append(row)
        else:
            # Create a row for each ENSP associated with this gene
            for i, ensp in enumerate(ensp_list):
                row = {
                    'ENSG': gene_id,
                    'ENSP': ensp,
                    'uniprot_id': uniprot_mapping[gene_id]['uniprot_ids'][i],
                    'ec_number': uniprot_mapping[gene_id]['ec_numbers'][i],
                    'protein_seq': uniprot_mapping[gene_id]['sequences'][i],
                    'num_reactions': len(info['reactions']),
                    'enzyme_type': info['type'],
                    'reactions': ';'.join(info['reactions']),
                }
                # Add components if it's a complex
                if 'components' in info:
                    row['components'] = ';'.join(info['components'])
                else:
                    row['components'] = ''
                
                data.append(row)
    
    # Create DataFrame
    df_enzyme_info = pd.DataFrame(data)
    
    return df_enzyme_info, enzyme_type_summary

def process_human_gpr(reaction, enzyme_dict, uniprot_mapping):
    """
    Process gene-protein-reaction association for a reaction in a human GEM.
    This function handles the ENSG-ENSP mapping for human genes.
    """
    # Extract GPR rule
    gpr = reaction.gene_reaction_rule.lower()
    
    # Simple case: single gene
    if len(reaction.genes) == 1:
        gene = list(reaction.genes)[0]
        gene_id = gene.id  # This is ENSG in Human1 model
        
        # Initialize enzyme entry if it doesn't exist
        if gene_id not in enzyme_dict:
            enzyme_dict[gene_id] = {
                'reactions': [],
                'type': 1  # Will be updated later if needed
            }
        
        # Add reaction information
        if reaction.id not in enzyme_dict[gene_id]['reactions']:
            enzyme_dict[gene_id]['reactions'].append(reaction.id)
    
    # Complex case: multiple genes (potentially heteromeric complex)
    else:
        # Extract genes involved in the reaction
        genes = list(reaction.genes)
        gene_ids = [gene.id for gene in genes]
        
        # Check if it's an OR relationship (isozymes) or AND relationship (complex)
        if ' or ' in gpr:
            # Isozymes case (promiscuous)
            for gene in genes:
                gene_id = gene.id
                
                # Initialize the enzyme entry if it doesn't exist
                if gene_id not in enzyme_dict:
                    enzyme_dict[gene_id] = {
                        'reactions': [],
                        'type': 2  # Promiscuous, will be updated later if needed
                    }
                
                # Add reaction information
                if reaction.id not in enzyme_dict[gene_id]['reactions']:
                    enzyme_dict[gene_id]['reactions'].append(reaction.id)
        
        elif ' and ' in gpr:
            # For each gene in the complex
            for gene in genes:
                gene_id = gene.id
                
                # Initialize the enzyme entry if it doesn't exist
                if gene_id not in enzyme_dict:
                    enzyme_dict[gene_id] = {
                        'reactions': [],
                        'type': 4,  # Heteromeric complex by default, will be updated if homomeric
                        'components': []    # Store component genes
                    }
                # Ensure components list exists even if the entry already exists
                elif 'components' not in enzyme_dict[gene_id]:
                    enzyme_dict[gene_id]['components'] = []
                
                # Add all other genes in the complex as components
                for other_gene_id in gene_ids:
                    if other_gene_id != gene_id and other_gene_id not in enzyme_dict[gene_id]['components']:
                        enzyme_dict[gene_id]['components'].append(other_gene_id)
                
                # Add reaction information
                if reaction.id not in enzyme_dict[gene_id]['reactions']:
                    enzyme_dict[gene_id]['reactions'].append(reaction.id)
        
        else:
            # Handle more complex GPR expressions (if needed)
            pass

def calculate_human_kcats(df_enzyme_info, reaction_substrates_df):
    """
    Create a DataFrame with expanded reactions and substrates for each ENSP-ENSG pair.
    
    Parameters:
    -----------
    df_enzyme_info : pd.DataFrame
        DataFrame containing enzyme information with ENSP, ENSG, uniprot_id, and reactions columns.
        
    reaction_substrates_df : pd.DataFrame
        DataFrame containing reaction-substrate pairs with columns 'reaction_id' and 'metabolite_id'.
        
    Returns:
    --------
    pd.DataFrame
        A DataFrame with columns 'ENSG', 'ENSP', 'uniprot_id', 'reaction_id', and 'metabolite_id'.
    """
    # List to store expanded data
    expanded_data = []
    
    # Iterate over each row in df_enzyme_info
    for _, row in df_enzyme_info.iterrows():
        ensg = row['ENSG']
        ensp = row['ENSP']
        uniprot_id = row['uniprot_id']
        reactions = row['reactions'].split(';')
        
        # For each reaction, find corresponding substrates
        for reaction_id in reactions:
            substrates = reaction_substrates_df[reaction_substrates_df['reaction_id'] == reaction_id]['metabolite_id']
            for substrate in substrates:
                expanded_data.append({
                    'ENSG': ensg,
                    'ENSP': ensp,
                    'uniprot_id': uniprot_id,
                    'reaction_id': reaction_id,
                    'metabolite_id': substrate
                })
    
    # Create a new DataFrame from the expanded data
    kcat_df = pd.DataFrame(expanded_data)
    
    return kcat_df