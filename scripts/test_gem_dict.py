from GEM_dict import process_gem
import pandas as pd

def main():
    # Test with iML1515 model
    model_id = 'iML1515'
    
    # Path to your UniProt mapping file
    uniprot_mapping_file = 'results/iML1515_E coli_83333_UniProt - Copy.csv'
    
    # Process the GEM
    df_enzyme_info, enzyme_type_summary = process_gem(model_id, uniprot_mapping_file)
    
    if df_enzyme_info is None:
        print("Failed to process the model")
        return
    
    # Print summary of enzyme types
    print("\nEnzyme Type Summary:")
    print(f"Type 1 (Homomeric): {enzyme_type_summary[1]}")
    print(f"Type 2 (Promiscuous): {enzyme_type_summary[2]}")
    print(f"Type 3 (Homomeric Complex): {enzyme_type_summary[3]}")
    print(f"Type 4 (Heteromeric Complex): {enzyme_type_summary[4]}")
    
    # Print first few rows of the DataFrame
    print("\nFirst few rows of the DataFrame:")
    print(df_enzyme_info.head())
    
    # Save DataFrame to CSV
    output_file = 'results/enzyme_info_Ecoli.csv'
    df_enzyme_info.to_csv(output_file)
    print(f"\nResults saved to {output_file}")
    
    # Print some basic statistics
    print("\nBasic Statistics:")
    print(f"Total number of enzymes: {len(df_enzyme_info)}")
    print(f"Number of enzymes with UniProt IDs: {df_enzyme_info['uniprot_id'].notna().sum()}")
    print(f"Number of enzymes with EC numbers: {df_enzyme_info['ec_number'].notna().sum()}")
    print(f"Number of enzymes with protein sequences: {df_enzyme_info['protein_seq'].notna().sum()}")

if __name__ == "__main__":
    main() 