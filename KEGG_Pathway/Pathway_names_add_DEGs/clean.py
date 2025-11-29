import pandas as pd

# Load CSV files with error handling
file1 = "up_DEGs.xlsx"
file2 = "KO_to_Pathway_Output.csv"

try:
    df1 = pd.read_excel(file1)
    df2 = pd.read_csv(file2)
except FileNotFoundError as e:
    print(f"Error: {e}")
    exit()
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    exit()

# Rename columns for clarity
df2.rename(columns={df2.columns[0]: "KEGG_ko", df2.columns[1]: "Pathway_name"}, inplace=True)

# Check if required columns exist
if "KEGG_ko" not in df1.columns or "KEGG_ko" not in df2.columns:
    print("Error: 'KEGG_ko' column is missing in one of the files.")
    exit()

# Merge on the matching column
merged_df = df1.merge(df2, on="KEGG_ko", how="left")
print("Merging completed. Saving merged output...")
merged_df.to_csv("merged_output_pathway.csv", index=False)

# Define the column to split
column_to_split = "Pathway_name"
if column_to_split not in merged_df.columns:
    print(f"Error: Column '{column_to_split}' not found in the merged dataframe.")
    exit()

# Split the column by semicolon and explode into multiple rows
merged_df[column_to_split] = merged_df[column_to_split].str.split(";")
exploded_df = merged_df.explode(column_to_split)
exploded_df.reset_index(drop=True, inplace=True)
print("Splitting and exploding completed. Saving final output...")
exploded_df.to_csv("final_splited_pathway_output.csv", index=False)

# Define columns to clean
column_to_clean_1 = "KEGG_ko"
column_to_clean_2 = "Pathway_name"

# Remove rows with invalid or unwanted values
filtered_df = exploded_df[~exploded_df[column_to_clean_1].str.strip().isin(["#NA", "-"])]
cleaned_df = filtered_df[
    (filtered_df[column_to_clean_2].str.strip().str.lower() != "no pathway found") &
    (filtered_df[column_to_clean_2].notna()) &
    (filtered_df[column_to_clean_2].str.strip() != "")
]
cleaned_df.reset_index(drop=True, inplace=True)
print("Cleaning completed. Saving cleaned output...")
cleaned_df.to_csv("cleaned_pathway_output.csv", index=False)

print("Cleaning and splitting completed successfully!")