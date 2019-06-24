import math
import numpy as np
import pandas as pd
import statsmodels.sandbox.stats.multicomp as mc
import statsmodels.api as sm
import statsmodels.formula.api as smf
from  scipy.stats import f

# Calculates and returns the dispersion (scale) parameter to use in GLM.
def compute_scale(ms, vs):
	res = 0
	c = 0
	for i in range(0, len(ms)):
		r = ((vs[i] - ms[i]) / (ms[i]**2))
		if r < 0:
			r = 0
		if not np.isnan(r):
			res += r
			c += 1
	if (res==0) or (c==0):
		return 0
	return res/float(c)

def run(counts_table, model_design, model_columns=None, experimental_column=None, log=False):
	counts_table.sort_values(by=['unique_id'], inplace=True) # Sort using annotation information.

	sample_names = list(model_design.index) # pull out sample names from the model.
	
	if model_columns is not None: # if user specifies the model design file columns to use.
		variable_names = [c for c in model_design.columns if c in model_columns]
	else: # default is to use all columns in model design file.
		variable_names = list(model_design.columns)
	variable_indicators = [list(model_design[c]) for c in variable_names] # sample indicators for each parameter.

	columns_to_keep = list(counts_table[sample_names+['unique_id']].set_index('unique_id').dropna(how='all').index)
	counts_table = counts_table[counts_table['unique_id'].isin(columns_to_keep)] # remove unwanted columns from counts table.

	if experimental_column is None: # if the experimental column is not already specified
		experimental_column = variable_names[-1]

	assert experimental_column in variable_names

	experimental_formula = "counts ~ " + " + ".join(variable_names) + " + barcode"
	if (len(variable_names)-1) != 0:
		control_formula = "counts ~ " + " + ".join([v for v in variable_names if v != experimental_column]) + " + barcode"
	else:
		control_formula = "counts ~ barcode"

	result_list = []
	if log: 
		failed_df = pd.DataFrame(columns=counts_table.columns)
	
	index = 0
	while index < len(counts_table):
		row = counts_table.iloc[index, :] # Get current row.
		unique_id = row['unique_id'] # Get the annotated reference.

		if pd.isnull(unique_id): # If barcode cannot be paired to feature.
			unique_id = row.name
			barcode_count = 1
		else:
			barcode_count = ((counts_table['unique_id'] == unique_id)).sum()  # Number of barcodes with same unique_id.
		end_index = index + barcode_count # Increment by the number of barcodes processed during this iteration.

		row_subset = counts_table.iloc[index:end_index] # Take the subset of rows to model collectively.
		assert list(row_subset['unique_id'])[1:] == list(row_subset['unique_id'])[:-1]

		simple_barcodes = list(row_subset.index)
		counts = [c for r in row_subset[sample_names].values.tolist() for c in r]
		row_df = pd.DataFrame([v*barcode_count for v in variable_indicators]).T
		row_df.columns = variable_names

		bcol = []
		for b in simple_barcodes:
			bcol.extend([b]*len(sample_names))
		row_df['barcode'] = bcol
		row_df['counts'] = counts

		exp_means = list(row_df.groupby(variable_names + ["barcode"])['counts'].mean()) 
		exp_vars = list(row_df.groupby(variable_names + ["barcode"])['counts'].var())
		exp_scale = compute_scale(exp_means, exp_vars)
		row_df.set_index('counts', inplace=True)

		control_indicator = min(row_df[experimental_column])
		experimental_indicator = max(row_df[experimental_column])

		cm = np.nanmean(list(row_df[(row_df[experimental_column]==min(row_df[experimental_column]))].index))
		em = np.nanmean(list(row_df[(row_df[experimental_column]==max(row_df[experimental_column]))].index))
		percentNA = len([0 for x in counts if math.isnan(x)]) / float(len(counts))

		index = end_index

		try:
			model = smf.glm(experimental_formula, 
				data=row_df,
				family=sm.families.NegativeBinomial(alpha=exp_scale)).fit()
			control_model = smf.glm(control_formula, 
				data=row_df,
				family=sm.families.NegativeBinomial(alpha=exp_scale)).fit()
			fscore = abs(model.deviance - control_model.deviance) / (1.) # One degree of freedom.
			pval_f = 1 - f.cdf(fscore, 1, model.df_resid) # Degrees of freedom = df_resid.
		except:
			if log: 
				failed_df = pd.concat([failed_df, row_subset])
			continue		

		result_list.append([unique_id, row['sseqid'], simple_barcodes, model.params[len(variable_names)],
			pval_f, percentNA, cm, em, exp_scale, 
			row['library'], row['subpool'], row['guide_status'], row['donor_status'], row['bsp_status'], len(simple_barcodes)])

	result_table = pd.DataFrame(result_list, columns=
		['unique_id', 'feature', 'barcodes', 'coef',
		'pval_f', 'percentNA', 'mean_control', 'mean_exp', 'scale', 
	 	'library', 'subpool', 'guide_status', 'donor_status', 'bsp_status', 'num_barcodes'])

	rm11 = result_table[result_table['library']=="RM11_V3"]
	sk1 = result_table[result_table['library']=="SK1_V3"]
	res_rm11 = mc.multipletests(pvals=list(rm11['pval_f']), method='fdr_bh') # Perform BH FDR adjustment.
	res_sk1 = mc.multipletests(pvals=list(sk1['pval_f']), method='fdr_bh') # Perform BH FDR adjustment.
	rm11['fdr_p_sep'] =  res_rm11[1]
	sk1['fdr_p_sep'] =  res_sk1[1]
	result_table = pd.concat([rm11, sk1])

	res = mc.multipletests(pvals=list(result_table['pval_f']), method='fdr_bh') # Perform BH FDR adjustment.
	result_table['fdr_p_all'] = res[1] # FDR adjusted p-values using BH.

	if log:
		return result_table, failed_df
	else:
		return result_table