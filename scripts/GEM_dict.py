import cobra
import pandas as pd
import os

def process_gem(model, uniprot_mapping_file, organism="E coli"):
    """
    Process a GEM to extract enzyme information using a pre-generated UniProt mapping file.
    
    Parameters:
    -----------
    model : cobra.Model
        The genome-scale metabolic model.
    uniprot_mapping_file : str
        Path to a file mapping model gene IDs to UniProt IDs, EC numbers, and sequences.
        This file is mandatory and should contain columns: model_gene_id, uniprot_id, ec_number, sequence.
    organism : str, optional
        Organism name for reference (default: "E coli")
        
    Returns:
    --------
    dict
        A dictionary with enzyme IDs as keys and dictionaries with enzyme 
        information as values.
    """
    # Initialize the output dictionary
    enzyme_dict = {}
    
    # Load UniProt mapping (mandatory)
    if not os.path.exists(uniprot_mapping_file):
        raise FileNotFoundError(f"UniProt mapping file not found: {uniprot_mapping_file}")
    
    # Load the mapping data
    uniprot_df = pd.read_csv(uniprot_mapping_file)
    
    # Create a dictionary for quick lookups - mapping model gene IDs to UniProt info
    uniprot_mapping = {}
    for _, row in uniprot_df.iterrows():
        uniprot_mapping[row['model_gene_id']] = {
            'uniprot_id': row['uniprot_id'],
            'ec_number': row['ec_number'],
            'sequence': row['sequence']
        }
    
    # Process each reaction in the model
    for reaction in model.reactions:
        # Skip reactions without genes
        if not reaction.genes:
            continue
        
        # Process the gene-protein-reaction association
        process_gpr(reaction, enzyme_dict, uniprot_mapping)
    
    # Classify enzyme types
    classify_enzyme_types(enzyme_dict)
    
    return enzyme_dict

def process_gpr(reaction, enzyme_dict, uniprot_mapping):
    """
    Process gene-protein-reaction association for a reaction.
    """
    # Extract GPR rule
    gpr = reaction.gene_reaction_rule.lower()
    
    # Simple case: single gene (homomeric enzyme)
    if len(reaction.genes) == 1:
        gene = list(reaction.genes)[0]
        gene_id = gene.id
        
        # Get UniProt info from mapping
        if gene_id in uniprot_mapping:
            uniprot_id = uniprot_mapping[gene_id]['uniprot_id']
            sequence = uniprot_mapping[gene_id]['sequence']
            ec_number = uniprot_mapping[gene_id]['ec_number']
        else:
            # If gene is not in mapping, use gene ID as fallback
            uniprot_id = gene_id
            sequence = ""
            ec_number = ""
        
        if uniprot_id not in enzyme_dict:
            enzyme_dict[uniprot_id] = {
                'protein_seq': sequence,
                'ec_number': ec_number,
                'substrates': [],
                'reactions': [],
                'type': 1  # Will be updated later if needed
            }
        
        # Add reaction information
        if reaction.id not in enzyme_dict[uniprot_id]['reactions']:
            enzyme_dict[uniprot_id]['reactions'].append(reaction.id)
        
        # Add substrate information
        substrates = [m.id for m in reaction.reactants]
        for substrate in substrates:
            if substrate not in enzyme_dict[uniprot_id]['substrates']:
                enzyme_dict[uniprot_id]['substrates'].append(substrate)
    
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
                
                # Get UniProt info from mapping
                if gene_id in uniprot_mapping:
                    uniprot_id = uniprot_mapping[gene_id]['uniprot_id']
                    sequence = uniprot_mapping[gene_id]['sequence']
                    ec_number = uniprot_mapping[gene_id]['ec_number']
                else:
                    # If gene is not in mapping, use gene ID as fallback
                    uniprot_id = gene_id
                    sequence = ""
                    ec_number = ""
                
                if uniprot_id not in enzyme_dict:
                    enzyme_dict[uniprot_id] = {
                        'protein_seq': sequence,
                        'ec_number': ec_number,
                        'substrates': [],
                        'reactions': [],
                        'type': 2  # Promiscuous, will be updated later if needed
                    }
                
                # Add reaction information
                if reaction.id not in enzyme_dict[uniprot_id]['reactions']:
                    enzyme_dict[uniprot_id]['reactions'].append(reaction.id)
                
                # Add substrate information
                substrates = [m.id for m in reaction.reactants]
                for substrate in substrates:
                    if substrate not in enzyme_dict[uniprot_id]['substrates']:
                        enzyme_dict[uniprot_id]['substrates'].append(substrate)
        
        elif ' and ' in gpr:
            # Complex case - create an entry for the complex
            complex_id = '_'.join(sorted(gene_ids))
            
            if complex_id not in enzyme_dict:
                enzyme_dict[complex_id] = {
                    'protein_seq': "",  # Complexes don't have a single sequence
                    'ec_number': "",    # Complexes don't have a single EC number
                    'substrates': [],
                    'reactions': [],
                    'type': 4,  # Heteromeric complex by default, will be updated if homomeric
                    'components': []    # Store component genes
                }
            
            # Add component genes
            for gene_id in gene_ids:
                if gene_id in uniprot_mapping:
                    uniprot_id = uniprot_mapping[gene_id]['uniprot_id']
                else:
                    uniprot_id = gene_id
                    
                if uniprot_id not in enzyme_dict[complex_id]['components']:
                    enzyme_dict[complex_id]['components'].append(uniprot_id)
            
            # Add reaction information
            if reaction.id not in enzyme_dict[complex_id]['reactions']:
                enzyme_dict[complex_id]['reactions'].append(reaction.id)
            
            # Add substrate information
            substrates = [m.id for m in reaction.reactants]
            for substrate in substrates:
                if substrate not in enzyme_dict[complex_id]['substrates']:
                    enzyme_dict[complex_id]['substrates'].append(substrate)
        
        else:
            # Handle more complex GPR expressions
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

# Example usage
if __name__ == "__main__":
    import cobra
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python script.py <model_id> <uniprot_mapping_file>")
        sys.exit(1)
    
    model_id = sys.argv[1]
    uniprot_mapping_file = sys.argv[2]
    
    try:
        # Load the model
        model = cobra.io.load_model(model_id)
        print(f"Model {model_id} loaded successfully.")
    except Exception:
        try:
            # Try loading as file path
            model = cobra.io.read_sbml_model(model_id)
            print(f"Model loaded from file: {model_id}")
        except Exception as e:
            print(f"Failed to load model: {e}")
            sys.exit(1)
    
    # Process the model
    enzyme_info = process_gem(model, uniprot_mapping_file)
    
    # Save results
    output_file = f"{model_id}_enzyme_info.csv"
    
    # Convert enzyme dictionary to DataFrame for easy saving
    rows = []
    for enzyme_id, info in enzyme_info.items():
        row = {
            'enzyme_id': enzyme_id,
            'protein_seq': info['protein_seq'],
            'ec_number': info['ec_number'],
            'substrates': ','.join(info['substrates']),
            'reactions': ','.join(info['reactions']),
            'type': info['type']
        }
        if 'components' in info:
            row['components'] = ','.join(info['components'])
        rows.append(row)
    
    result_df = pd.DataFrame(rows)
    result_df.to_csv(output_file, index=False)
    print(f"Enzyme information saved to {output_file}")