"""
Module Description: gene_sequence_mapper.py

Purpose: 
Recover protein sequences using the gene IDs annotated in genome-scale models.

Overview: 
This module provides functions to:
1. Map gene IDs from a genome-scale model (GEM) to UniProt IDs, including retrieving EC numbers and protein sequences.
2. Map Ensembl Gene IDs (ENSG) to Ensembl Protein IDs (ENSP) and further to UniProt IDs.
3. Map PaxDB IDs to UniProt IDs and sequences.
4. Save the mapping results to CSV files for further analysis.

"""

import cobra
import pandas as pd
import os
from bioservices import UniProt
from mygene import MyGeneInfo


def map_gem_genes_to_uniprot(model_id, organism):
    """
    Map gene IDs from a GEM to UniProt IDs.
    
    Parameters:
        model_id : str
            ID of the model (e.g., 'iML1515') to load using COBRA's model repository.
            Path to the model file if not using the cobra repository.
        organism : str
            Organism name for UniProt queries, must match model organism
        
    Returns:
        pandas.DataFrame
            DataFrame containing the mapping between model gene IDs and UniProt IDs,
            with columns: model_gene_id, uniprot_id, ec_number, sequence
    """
        
    # Set taxon ID based on organism
    # TO DO: make a file mapping model IDs to taxon IDs
    if organism.lower() == 'e coli':
        taxon_ID = '83333'
        #taxon_ID = '511145'
        print(f'Using taxon ID: {taxon_ID} - E. coli')
    elif organism.lower() == 'yeast':
        taxon_ID = '4932'
        print(f'Using taxon ID: {taxon_ID} - Yeast')
    else:
        # Abort if the organism is not specified
        raise ValueError(f"Organism {organism} not supported.")
        
    # Load the model
    try:
        # From the cobra repository
        model_read = str(model_id)
        model = cobra.io.load_model(model_read)
        print(f"Model {model_id} loaded successfully from COBRA repository.")
    except Exception:
        # Try alternate loading method (model ID is a file path)
        print(f"Loading model as path...")
        try:
            model = cobra.io.read_sbml_model(model_id)
            print(f"Model {model_id} loaded successfully from path.")
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
    output_file = os.path.join('..', 'results', f"{model_id}_{organism}_{taxon_ID}_UniProt.csv")
    mapping_df.to_csv(output_file, index=False)
    print(f"Created mapping file {output_file}")
    
    return mapping_df


def ENSG_to_ENSP(model_id):
    """
    Map ENSG IDs to ENSP IDs for a given model.
    
    Parameters:
        model_id : str
            ID of the model (e.g., 'iML1515') to load using COBRA's model repository.
            Path to the model file if not using the cobra repository.
        
    Returns:
        pandas.DataFrame    
    """
    
    # Load the model
    try:
        # From the cobra repository
        model_read = str(model_id)
        model = cobra.io.load_model(model_read)
        print(f"Model {model_id} loaded successfully from COBRA repository.")
    except Exception:
        # Try alternate loading method (model ID is a file path)
        print(f"Loading model as path...")
        try:
            model = cobra.io.read_sbml_model(model_id)
            print(f"Model {model_id} loaded successfully from path.")
        except Exception as e2:
            print(f"Failed to load model: {e2}")
            return None
    
    # Extract all unique gene IDs from the model
    ensg_list = set()
    for reaction in model.reactions:
        for gene in reaction.genes:
            ensg_list.add(gene.id)
            
    mg = MyGeneInfo()
    results = mg.querymany(ensg_list, scopes="ensembl.gene", fields="ensembl.protein", species="human")
    
    records = []

    for res in results:
        ensg = res["query"]
        ensembl_data = res.get("ensembl", {})
        
        # Handle both dict and list cases
        if isinstance(ensembl_data, dict):
            proteins = ensembl_data.get("protein")
            if proteins:
                if isinstance(proteins, list):
                    for ensp in proteins:
                        records.append({"ENSP": ensp, "ENSG": ensg})
                else:
                    records.append({"ENSP": proteins, "ENSG": ensg})
        elif isinstance(ensembl_data, list):
            for entry in ensembl_data:
                protein = entry.get("protein")
                if protein:
                    if isinstance(protein, list):
                        for ensp in protein:
                            records.append({"ENSP": ensp, "ENSG": ensg})
                    else:
                        records.append({"ENSP": protein, "ENSG": ensg})
    
    ENSP_df = pd.DataFrame(records)
    
    # Save the mapping file
    output_file = os.path.join('..', 'results', f"{model_id}_ENSP_to_ENSG.csv")
    ENSP_df.to_csv(output_file, index=False)
    print(f"Created mapping file {output_file}")
    
    return ENSP_df



def map_ENSP_to_UniProt(ENSP_df, taxon_ID="9606"):
    """
    Map ENSP IDs to UniProt IDs and sequences for a given model.
    
    Parameters:
        ENSP_df : pandas.DataFrame
            DataFrame containing the mapping between model gene IDs and UniProt IDs,
            with columns: model_gene_id, uniprot_id, ec_number, sequence
         taxon_ID : str
            Taxon ID for UniProt queries, must match model organism
        
    Returns:
        pandas.DataFrame
    """
        
    # Initialize UniProt service
    service = UniProt()
    
    # Get unique ENSPs
    all_protein_ids = ENSP_df["ENSP"].dropna().unique()
    total_proteins = len(all_protein_ids)
    mapping_data = []

    print(f"Starting mapping of {total_proteins} ENSP IDs to UniProt entries...")

    for idx, protein_id in enumerate(all_protein_ids, 1):
        if idx % 50 == 0 or idx == total_proteins:
            print(f"Progress: {idx}/{total_proteins} proteins processed")

        uniprot_id = None
        ec_number = None
        sequence = None

        try:
            # Query by Ensembl Protein ID (ENSP)
            query = f"xref:ensembl-{protein_id} AND organism_id:{taxon_ID}"
            print(query)  # Print the query for debugging
            result = service.search(query, frmt="tsv", columns="id,ec")

            if result and "\n" in result:
                lines = result.strip().split("\n")
                if len(lines) > 1:  # Skip header
                    fields = lines[1].split("\t")
                    if len(fields) >= 1:
                        uniprot_id = fields[0]
                        ec_number = fields[1] if len(fields) > 1 else None

            # If first query failed, try alternative format
            if not uniprot_id:
                query = f"database:ensembl AND xref:{protein_id} AND organism_id:{taxon_ID}"
                print(query)  # Print the query for debugging
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
                        sequence = fasta_result.split("\n", 1)[-1].replace("\n", "").strip()
                except Exception as e:
                    print(f"Error fetching sequence for {uniprot_id}: {e}")
        except Exception as e:
            print(f"Error looking up {protein_id}: {e}")

        # If no mapping was found, use the original ENSP ID
        if not uniprot_id:
            uniprot_id = protein_id

        mapping_data.append({
            'ensp_id': protein_id,
            'uniprot_id': uniprot_id,
            'ec_number': ec_number if ec_number else '',
            'sequence': sequence if sequence else ''
        })

    # Create a DataFrame from the mapping data
    mapping_df = pd.DataFrame(mapping_data)

    # Merge the new mapping data with the original ENSP_df
    protein_mapping_df = pd.merge(ENSP_df, mapping_df, left_on='ENSP', right_on='ensp_id', how='left')

    # Save the mapping file
    output_file = os.path.join('..', 'data', f"Human1_{taxon_ID}_UniProt.csv")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    protein_mapping_df.to_csv(output_file, index=False)
    print(f"Created mapping file {output_file}")

    return protein_mapping_df


def map_PaxDB_to_UniProt(PaxDB_df, taxon_ID="9606"):
    """
    Map PaxDB IDs to UniProt IDs and sequences.
    
    Parameters:
        PaxDB_df : pandas.DataFrame
            DataFrame containing PaxDB IDs (must have a 'PaxDB_ID' column)
        taxon_ID : str
            Taxon ID for UniProt queries, must match model organism
        
    Returns:
        pandas.DataFrame
            DataFrame containing the mapping with all original columns plus UniProt data
    """
        
    # Initialize UniProt service
    service = UniProt()
    
    # Get unique PaxDB IDs
    all_protein_ids = PaxDB_df["PaxDB_ID"].dropna().unique()
    total_proteins = len(all_protein_ids)
    mapping_data = []

    print(f"Starting mapping of {total_proteins} PaxDB IDs to UniProt entries...")

    for idx, protein_id in enumerate(all_protein_ids, 1):
        if idx % 50 == 0 or idx == total_proteins:
            print(f"Progress: {idx}/{total_proteins} proteins processed")

        uniprot_id = None
        ec_number = None
        sequence = None

        try:
            # Try different query formats for PaxDB IDs
            # Option 1: Direct mapping if PaxDB ID is an Ensembl ID
            if protein_id.startswith("ENSP"):
                query = f"xref:ensembl-{protein_id} AND organism_id:{taxon_ID}"
            else:
                # Option 2: Try generic database cross-reference
                query = f"xref:{protein_id} AND organism_id:{taxon_ID}"
                
            result = service.search(query, frmt="tsv", columns="id,ec")

            if result and "\n" in result:
                lines = result.strip().split("\n")
                if len(lines) > 1:  # Skip header
                    fields = lines[1].split("\t")
                    if len(fields) >= 1:
                        uniprot_id = fields[0]
                        ec_number = fields[1] if len(fields) > 1 else None

            # If the above queries failed, try additional formats
            if not uniprot_id:
                # Try as a gene name
                query = f"gene:{protein_id} AND organism_id:{taxon_ID}"
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
                        sequence = fasta_result.split("\n", 1)[-1].replace("\n", "").strip()
                except Exception as e:
                    print(f"Error fetching sequence for {uniprot_id}: {e}")
        except Exception as e:
            print(f"Error looking up {protein_id}: {e}")

        # If no mapping was found, use the original PaxDB ID
        if not uniprot_id:
            uniprot_id = protein_id

        mapping_data.append({
            'PaxDB_ID': protein_id,
            'uniprot_id': uniprot_id,
            'ec_number': ec_number if ec_number else '',
            'sequence': sequence if sequence else ''
        })

    # Create a DataFrame from the mapping data
    mapping_df = pd.DataFrame(mapping_data)

    # Merge the new mapping data with the original PaxDB_df
    protein_mapping_df = pd.merge(PaxDB_df, mapping_df, on='PaxDB_ID', how='left')

    # Save the mapping file
    output_file = os.path.join('..', 'data', f"PaxDB_{taxon_ID}_UniProt.csv")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    protein_mapping_df.to_csv(output_file, index=False)
    print(f"Created mapping file {output_file}")

    return protein_mapping_df


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Map gene IDs from a GEM to UniProt IDs')
    parser.add_argument('--model_id', type=str, required=True, 
                       help='Model ID (e.g., iML1515) or SBML file path')
    parser.add_argument('--organism', type=str, required=True, 
                       help='Organism name (e.g., E coli)')    
    args = parser.parse_args()
    
    # Run the mapping function
    map_gem_genes_to_uniprot(args.model_id, args.organism)