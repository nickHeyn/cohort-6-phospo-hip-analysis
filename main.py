import pandas as pd

PHOSPHO_FILE_PATH = "/home/nick/Documents/Repos/cohort-6-phospo-hip-analysis/input/female_Hipp_phospho_combinedDE(Proteins).csv"
TOTAL_FILE_PATH = "/home/nick/Documents/Repos/cohort-6-phospo-hip-analysis/input/female_Hipp_totalDE_p05(Proteins).csv"

DEFAULT_ABUNDANCE_COL_NAMES = ["Abundance Ratio P-Value: (Female, Acute) / (Female, Binge)", "Abundance Ratio P-Value: (Female, Acute) / (Female, Control)", "Abundance Ratio P-Value: (Female, Binge) / (Female, Control)"]
SIG_PHOSPHO_OUTPUT_COLUMN_NAMES = ["Gene ID"] + DEFAULT_ABUNDANCE_COL_NAMES
DEFAULT_ABUNDANCE_THRESHOLD = 0.05


def locate_sig_abundance_p_val(df, sig_threshold, col_names):
    conditions = []
    for col_name in col_names:
        conditions.append(df[col_name] <= sig_threshold)

    combined_condition = conditions[0]
    for cond in conditions[1:]:
        combined_condition = combined_condition | cond

    
    return df[combined_condition]

df_phospho = pd.read_csv(PHOSPHO_FILE_PATH)
df_total = pd.read_csv(TOTAL_FILE_PATH)

sig_phospho = locate_sig_abundance_p_val(df_phospho, DEFAULT_ABUNDANCE_THRESHOLD, DEFAULT_ABUNDANCE_COL_NAMES)
sig_phospho.to_csv("output/sig_phospho_output.csv", index=False)

