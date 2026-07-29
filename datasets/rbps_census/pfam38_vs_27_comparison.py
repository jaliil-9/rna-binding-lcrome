import pandas as pd

# 1. File paths (Use raw strings to avoid escape sequence issues)
candidates_path = r"datasets/rbps_census/ensembl_v116/ensembl116_pfam38_selected_proteins.xlsx"
census_path = r"datasets/rbps_census/gerstberg/rna_binding_proteins.xls"
output_path = r"pfam38_census_comparison_v2.xlsx"

# 2. Load and combine both candidate sheets ("Strict" and "Borderline")
xls = pd.ExcelFile(candidates_path)
df_strict = pd.read_excel(xls, sheet_name="Strict")
df_borderline = pd.read_excel(xls, sheet_name="Borderline")

# (Optional) Tag them before combining so you know where they came from
df_strict['candidate_class'] = "Strict"
df_borderline['candidate_class'] = "Borderline"

df_combined = pd.concat([df_strict, df_borderline], ignore_index=True)

# 3. Create a clean mapping key by removing the version from the ENSP ID (e.g., ENSP000...7390.2 -> ENSP000...7390)
df_combined['ensp_base'] = df_combined['query_name'].astype(str).str.split('.').str[0]

# 4. Create Deduplicated subset (based on the base ENSP ID)
df_dedup = df_combined.drop_duplicates(subset=['ensp_base']).copy()

# 5. Load 2014 census and extract its ENSP IDs
df_census = pd.read_excel(census_path, sheet_name="RBP table")
# The 2014 census uses 'protein id' for its ENSP identifiers
census_ensps = set(df_census['protein id'].dropna().astype(str))

# 6. Keep only deduplicated rows whose ensp_base is NOT in the 2014 census set
df_not_in_2014 = df_dedup[~df_dedup['ensp_base'].isin(census_ensps)]

# 7. Export to three distinct sheets
with pd.ExcelWriter(output_path) as writer:
    df_combined.to_excel(writer, sheet_name="1_Combined_Original", index=False)
    df_dedup.to_excel(writer, sheet_name="2_Deduplicated", index=False)
    df_not_in_2014.to_excel(writer, sheet_name="3_Not_In_2014", index=False)

print(f"Export complete. File saved to {output_path}")
print(f"Total Combined: {len(df_combined)} | Deduplicated: {len(df_dedup)} | Not in 2014: {len(df_not_in_2014)}")