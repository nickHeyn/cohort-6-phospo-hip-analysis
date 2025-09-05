import pandas as pd

PHOSPHO_FILE_PATH = "/home/nick/Documents/Repos/cohort-6-phospo-hip-analysis/input/female_Hipp_phospho_combinedDE(Proteins).csv"
TOTAL_FILE_PATH = "/home/nick/Documents/Repos/cohort-6-phospo-hip-analysis/input/female_Hipp_totalDE_p05(Proteins).csv"

# Columns
GENE_ID = "Gene ID"
PHOSPHO_ABUND_P_VAL_A_B = "Abundance Ratio P-Value: (Female, Acute) / (Female, Binge)"
PHOSPHO_ABUND_P_VAL_A_C = "Abundance Ratio P-Value: (Female, Acute) / (Female, Control)"
PHOSPHO_ABUND_P_VAL_B_C = "Abundance Ratio P-Value: (Female, Binge) / (Female, Control)"
TOTAL_ABUND_P_VAL_A_B = "Abundance Ratio P-Value: (Binge) / (Acute)"
TOTAL_ABUND_P_VAL_A_C = "Abundance Ratio P-Value: (Acute) / (Control)"
TOTAL_ABUND_P_VAL_B_C = "Abundance Ratio P-Value: (Binge) / (Control)"

DEFAULT_PHOSPHO_ABUNDANCE_COL_NAMES = [PHOSPHO_ABUND_P_VAL_A_B, PHOSPHO_ABUND_P_VAL_A_C, PHOSPHO_ABUND_P_VAL_B_C]
DEFAULT_COL_NAME_PAIRS = [(TOTAL_ABUND_P_VAL_A_B, PHOSPHO_ABUND_P_VAL_A_B), (TOTAL_ABUND_P_VAL_A_C, PHOSPHO_ABUND_P_VAL_A_C), (TOTAL_ABUND_P_VAL_B_C, PHOSPHO_ABUND_P_VAL_B_C)]
SIG_PHOSPHO_OUTPUT_COLUMN_NAMES = [GENE_ID] + DEFAULT_PHOSPHO_ABUNDANCE_COL_NAMES
DEFAULT_ABUNDANCE_THRESHOLD = 0.05

DEFAULT_ABUNDANCE_P_VAL_PERCENTAGE_DIFF_THRRESHOLD = 5

def get_percentage_diff(valA, valB):
    diff = abs(valA - valB)
    avg = (valA + valB / 2.0)
    return (diff / avg) * 100


def locate_sig_abundance_p_val(df, sig_threshold, col_names):
    conditions = []
    for col_name in col_names:
        conditions.append(df[col_name] <= sig_threshold)

    combined_condition = conditions[0]
    for cond in conditions[1:]:
        combined_condition = combined_condition | cond

    
    return df[combined_condition]

def locate_sig_diff_abundance_p_val(df_total, df_phospho, col_name_pairs, target_diff_percentage):
    gene_ids = df_total[GENE_ID].to_list()
    result_column_names = ['Gene ID', 'Diff percentage', 'Targeted Total Col Name', 'Targeted Phospho Col Name']
    result_df = pd.DataFrame(columns=result_column_names)
    result_gene_ids = set()
    missing_gene_ids = set()
    for gene_id in gene_ids:
        total_row = df_total[df_total[GENE_ID] == gene_id]
        phospho_row = df_phospho[df_phospho[GENE_ID] == gene_id]

        if len(phospho_row) > 0:
            for col_name_pair in col_name_pairs:
                total_val = pd.to_numeric(total_row[col_name_pair[0]]).item()
                phospho_val = pd.to_numeric(phospho_row[col_name_pair[1]]).item()

                percentage_diff = get_percentage_diff(total_val, phospho_val)
                if percentage_diff <= target_diff_percentage:
                    result_gene_ids.add(gene_id)
                    result_df_row = pd.DataFrame(data=[[gene_id, percentage_diff, col_name_pair[0], col_name_pair[1]]], columns=result_column_names)
                    result_df = pd.concat([result_df, result_df_row], ignore_index=True)
        else:
            missing_gene_ids.add(gene_id)
    return result_df

            

df_phospho = pd.read_csv(PHOSPHO_FILE_PATH)
df_total = pd.read_csv(TOTAL_FILE_PATH)

sig_diff_df = locate_sig_diff_abundance_p_val(df_total, df_phospho, DEFAULT_COL_NAME_PAIRS, DEFAULT_ABUNDANCE_P_VAL_PERCENTAGE_DIFF_THRRESHOLD)

sig_diff_df.to_csv("output/sig_diff_output.csv", index=False)

sig_phospho = locate_sig_abundance_p_val(df_phospho, DEFAULT_ABUNDANCE_THRESHOLD, DEFAULT_ABUNDANCE_COL_NAMES)
sig_phospho.to_csv("output/sig_phospho_output.csv", index=False)

