import pandas as pd

candidates_path = "datasets/rbps_census/ensembl_v116/ensembl116_pfam38_selected_proteins.xlsx"
census_path = "datasets/rbps_census/gerstberg/rna_binding_proteins.xls"
output_path = "pfam38_census_comparison.xlsx"

xls = pd.ExcelFile(candidates_path)
df_strict = pd.read_excel(xls, sheet_name="Strict")
df_borderline = pd.read_excel(xls, sheet_name="Borderline")

df_strict['candidate_class'] = "Strict"
df_borderline['candidate_class'] = "Borderline"

df_combined = pd.concat([df_strict, df_borderline], ignore_index=True)

df_combined['ensp_base'] = df_combined['query_name'].astype(str).str.split('.').str[0]

df_dedup = df_combined.drop_duplicates(subset=['ensp_base']).copy()

df_census = pd.read_excel(census_path, sheet_name="RBP table")
census_ensps = set(df_census['protein id'].dropna().astype(str))

df_not_in_2014 = df_dedup[~df_dedup['ensp_base'].isin(census_ensps)]

with pd.ExcelWriter(output_path) as writer:
    df_combined.to_excel(writer, sheet_name="1_Combined_Original", index=False)
    df_dedup.to_excel(writer, sheet_name="2_Deduplicated", index=False)
    df_not_in_2014.to_excel(writer, sheet_name="3_Not_In_2014", index=False)

print(f"Export complete. File saved to {output_path}")
print(f"Total Combined: {len(df_combined)} | Deduplicated: {len(df_dedup)} | Not in 2014: {len(df_not_in_2014)}")