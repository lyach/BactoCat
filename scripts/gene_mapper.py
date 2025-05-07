import cobra
import pandas as pd
import os
from bioservices import UniProt


def map_gem_genes_to_uniprot(model_id='iML1515', organism="E coli"):
    """
    Map gene IDs from a GEM to UniProt IDs.
    
    Parameters:
    -----------
    model_id : str
        ID of the model (e.g., 'iML1515') to load using COBRA's model repository.
    organism : str
        Organism name for UniProt queries, must match model organism
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the mapping between model gene IDs and UniProt IDs,
        with columns: model_gene_id, uniprot_id, ec_number, sequence
    """
    
    # Set default output file
    output_file = f"{model_id}_uniprot_mapping.csv"
    
    # Set taxon ID based on organism
    if organism.lower() == 'e coli':
        taxon_ID = '83333'
    elif organism.lower() == 'yeast':
        taxon_ID = '4932'
    else:
        # Default to E. coli if not specified
        taxon_ID = '83333'
        print(f"Using E. coli taxon ID (83333) as GEM wasn't specified. Organism is not specified.")
    
    try:
        # Load the model from the cobra repository
        model_read = str(model_id)
        model = cobra.io.load_model(model_read)
        print(f"Model {model_id} loaded successfully.")
    except Exception:
        print(f"Loading model as path...")
        try:
            # Try alternate loading method (model ID is a file path)
            model = cobra.io.read_sbml_model(model_id)
            print(f"Model loaded from file: {model_id}")
        except Exception as e2:
            print(f"Failed to load model: {e2}")
            return None
    
    # Extract all unique gene IDs from the model
    all_gene_ids = set()
    for reaction in model.reactions:
        for gene in reaction.genes:
            all_gene_ids.add(gene.id)
    
    print(f"Found {len(all_gene_ids)} unique gene IDs in model {model_id}")
    
    # Initialize UniProt service
    service = UniProt()
    
    # Map gene IDs to UniProt IDs
    mapping_data = []
    
    # Track progress
    total_genes = len(all_gene_ids)
    processed = 0
    
    print(f"Starting mapping of {total_genes} gene IDs to UniProt IDs...")
    
    for gene_id in all_gene_ids:
        processed += 1
        if processed % 50 == 0 or processed == total_genes:
            print(f"Progress: {processed}/{total_genes} genes processed")
        
        uniprot_id = None
        ec_number = None
        sequence = None
        
        try:
            # Try mapping based on gene name to get UniProt ID and EC number
            query = f"gene_exact:({gene_id}) AND taxonomy_id:({taxon_ID})"
            result = service.search(query, frmt="tsv", columns="id,ec")
            
            if result and "\n" in result:
                lines = result.strip().split("\n")
                if len(lines) > 1:  # Skip header
                    fields = lines[1].split("\t")
                    if len(fields) >= 1:
                        uniprot_id = fields[0]
                        ec_number = fields[1] if len(fields) > 1 else None
            
            # Get sequence using FASTA format
            if uniprot_id:
                try:
                    fasta_query = f"id:{uniprot_id}"
                    fasta_result = service.search(fasta_query, frmt="fasta")
                    if fasta_result:
                        # Extract sequence part from FASTA format
                        sequence = fasta_result.split("\n", 1)[-1]
                        sequence = sequence.replace("\n", "")
                        sequence = sequence.strip()
                except Exception as e:
                    print(f"Error fetching sequence for {uniprot_id}: {e}")
        except Exception as e:
            print(f"Error looking up {gene_id}: {e}")
        
        # If no mapping was found, use the original gene ID
        if not uniprot_id:
            uniprot_id = gene_id
            
        # Add the entry to our mapping data
        mapping_data.append({
            'model_gene_id': gene_id,
            'uniprot_id': uniprot_id,
            'ec_number': ec_number if ec_number else '',
            'sequence': sequence if sequence else ''
        })
    
    # Create DataFrame from the mapping data we've already collected
    mapping_df = pd.DataFrame(mapping_data)
    
    # Save the mapping file
    mapping_df.to_csv(output_file, index=False)
    print(f"Created mapping file {output_file}")
    
    return mapping_df


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Map gene IDs from a GEM to UniProt IDs')
    parser.add_argument('--model_id', type=str, default='iML1515', help='Model ID (e.g., iML1515) or SBML file path')
    parser.add_argument('--organism', type=str, default='E coli', help='Organism name (default: E coli)')    
    args = parser.parse_args()
    
   # Run the mapping function
    map_gem_genes_to_uniprot(args.model_id, args.organism)