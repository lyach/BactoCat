"""
gene_sequence_mapper.py

Purpose: 
Recover protein sequences from UniProt using the gene IDs annotated in genome-scale models.

Overview: 
This module provides functions to:
1. Map organism names to UniProt taxon IDs, including retrieving EC numbers and protein sequences.
2. Map gene IDs from a genome-scale model (GEM) to UniProt IDs, including retrieving EC numbers and protein sequences.
2. Map Ensembl Gene IDs (ENSG) to Ensembl Protein IDs (ENSP) and further to UniProt IDs.
3. Map PaxDB IDs to UniProt IDs and sequences.
4. Save the mapping results to CSV files for further analysis.

"""

import cobra
import pandas as pd
import os
from bioservices import UniProt
from mygene import MyGeneInfo
from tqdm import tqdm

# Dictionary mapping lowercase organism names to representative GEM IDs
DEFAULT_GEM_IDS = {
    'e coli': 'iML1515',      # Representative E. coli model
    'ecoli': 'iML1515',       # Allow single word for E. coli
    'b subtilis': 'iBsu1103',  # Representative B. subtilis model
    'bsubtilis': 'iBsu1103',   # Allow single word for B. subtilis
    'p putida': 'iJN746',      # Representative P. putida model
    'pputida': 'iJN746',       # Allow single word for P. putida
    'm tuberculosis': 'iNJ661', # Representative M. tuberculosis model
    'mtuberculosis': 'iNJ661',  # Allow single word for M. tuberculosis
    'p aeruginosa': 'iPAU1129', # Representative P. aeruginosa model
    'paeruginosa': 'iPAU1129',  # Allow single word for P. aeruginosa
    's aureus': 'iTCS900',     # Representative S. aureus model
    'saureus': 'iTCS900',      # Allow single word for S. aureus
}

def map_organism_to_uniprot(organism: str) -> pd.DataFrame:
    """
    Map gene IDs from a representative GEM (inferred from organism name) to UniProt IDs, 
    using the organism name to determine the correct NCBI Taxonomy ID for search grounding.
    
    Parameters:
        organism : str
            Organism name for UniProt queries (e.g., 'E coli', 'B subtilis', or 'ecoli', 'bsubtilis'). 
            Must be one of the supported organisms with a known default GEM ID.
        
    Returns:
        pandas.DataFrame
            DataFrame containing the mapping between model gene IDs and UniProt IDs,
            with columns: model_gene_id, uniprot_id, ec_number, sequence
    """
    
    organism_lower = organism.lower()
    
    # --- 1. Map Organism to Taxon ID & Default Model ID ---
    
    # Logic to handle different organism name formats and assign Taxon ID
    if organism_lower == 'e coli' or organism_lower == 'ecoli':
        taxon_ID = '83333'
    elif organism_lower == 'b subtilis' or organism_lower == 'bsubtilis':
        taxon_ID = '1423'
    elif organism_lower == 'p putida' or organism_lower == 'pputida':
        taxon_ID = '303'
    elif organism_lower == 'm tuberculosis' or organism_lower == 'mtuberculosis':
        taxon_ID = '1773'
    elif organism_lower == 'p aeruginosa' or organism_lower == 'paeruginosa':
        taxon_ID = '287'
    elif organism_lower == 's aureus' or organism_lower == 'saureus':
        taxon_ID = '1280'
    else:
        raise ValueError(f"Organism '{organism}' not supported for automatic Taxon ID lookup.")
    
    print(f'Using taxon ID: {taxon_ID} - {organism}')
    
    # --- Simplified logic to get the default model ID ---
    # Find the key in DEFAULT_GEM_IDS by normalizing the input (removing spaces)
    normalized_organism = organism_lower.replace(' ', '')
    model_id = DEFAULT_GEM_IDS.get(organism_lower) or DEFAULT_GEM_IDS.get(normalized_organism)
    
    if not model_id:
         # This should never be reached if the if/elif above passed, but kept for safety
         raise ValueError(f"Internal Error: Could not find a default GEM ID for organism '{organism}'.")
    
    # --- 2. Load the Model using the inferred model_id ---
    try:
        # Try loading from the cobra repository using the model_id
        model = cobra.io.load_model(model_id)
        print(f"Model {model_id} loaded successfully from COBRA repository (inferred from organism).")
    except Exception:
        # Try alternate loading method (assuming model_id is a file path)
        print(f"Failed to load inferred model {model_id} from repository. Attempting to load as path...")
        try:
            model = cobra.io.read_sbml_model(model_id)
            print(f"Model {model_id} loaded successfully from path.")
        except Exception as e2:
            print(f"Failed to load model {model_id}: {e2}")
            return None
    
    # --- 3. Extract Genes and Prepare UniProt Mapping ---
    
    # Extract all unique gene IDs from the model
    all_gene_ids = set()
    for reaction in model.reactions:
        for gene in reaction.genes:
            all_gene_ids.add(gene.id)
    
    print(f"Found {len(all_gene_ids)} unique gene IDs in inferred model {model_id}")
    
    # Initialize UniProt service
    service = UniProt()
    
    # Map gene IDs to UniProt IDs
    mapping_data = []
    
    # --- 4. Loop Through Genes and Query UniProt ---
    
    for gene_id in tqdm(all_gene_ids, desc="Mapping gene IDs to UniProt"):
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
    
    # NOTE: File saving is moved to the calling pipeline function (step 4)
    return mapping_df


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
    
    for gene_id in tqdm(all_gene_ids, desc="Mapping gene IDs to UniProt"):
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
    mapping_data = []

    for protein_id in tqdm(all_protein_ids, desc="Mapping ENSP IDs to UniProt"):
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
    mapping_data = []

    for protein_id in tqdm(all_protein_ids, desc="Mapping PaxDB IDs to UniProt"):
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