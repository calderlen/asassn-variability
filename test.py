import pandas as pd
import os

# 1. The list of IDs extracted from your input text
target_ids = [
    'J183153-284827',
    'J070519+061219',
    'J081523-385923',
    'J085816-430955',
    'J114712-621037'
]

# 2. Load the manifest file
# Ensure 'lc_manifest.csv' is in the same directory as this script
try:
    df = pd.read_csv('lc_manifest.csv')
    print("Successfully loaded lc_manifest.csv")
    
    # 3. Identify the ID column
    # Adjust 'id_col_name' if the column in your CSV has a different name 
    # (common names: 'id', 'assassn_id', 'name')
    id_col_name = 'id' 
    if id_col_name not in df.columns:
        # Fallback: try to find a column that contains "id" or "name"
        potential_cols = [c for c in df.columns if 'id' in c.lower() or 'name' in c.lower()]
        if potential_cols:
            id_col_name = potential_cols[0]
            print(f"Assuming ID column is: '{id_col_name}'")
        else:
            print("Error: Could not find an ID column. Please edit the script to specify 'id_col_name'.")
            print("Available columns:", df.columns.tolist())
            exit()

    # 4. Filter for the target IDs
    # We use strip() to ensure no hidden spaces cause mismatches
    matches = df[df[id_col_name].astype(str).str.strip().isin(target_ids)]

    if matches.empty:
        print("\nNo matches found for these IDs.")
    else:
        print(f"\nFound {len(matches)} matches:\n")
        
        # 5. Print the location info
        # Checks for 'mag_bin' and 'lc_cal' columns to construct the path
        for index, row in matches.iterrows():
            star_id = row[id_col_name]
            
            # Try to grab path components if they exist
            mag_bin = row.get('mag_bin', 'N/A')
            lc_cal = row.get('lc_cal', 'N/A')
            
            # Construct the path
            path = os.path.join(str(mag_bin), str(lc_cal))
            
            print(f"ID: {star_id}")
            print(f"  -> Subdir Path: {path}")
            print("-" * 30)

except FileNotFoundError:
    print("Error: 'lc_manifest.csv' not found in the current directory.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")