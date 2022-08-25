#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import argparse
import csv
from copy import deepcopy
from cyvcf2 import VCF
import pickle as pkl

import sklearn
from sklearn.experimental import enable_iterative_imputer
# silence 'assigning to a copy of a slice' warning
pd.options.mode.chained_assignment = None  

version = 1.1

SV_min_size = 50
small_indel_max_size = 20

# sklearn settings
seed = 42
verbose = False
impute_iterations = 10
Epsilon_Support_Vector_Regression_regularization_parameter = 50
Epsilon_Support_Vector_Regression_epsilon = 0.1

# sklearn generic processing objects
imputor = sklearn.impute.IterativeImputer(max_iter=impute_iterations, verbose=verbose, random_state=seed)
scaler = sklearn.preprocessing.StandardScaler()
regressor = sklearn.svm.SVR(C=Epsilon_Support_Vector_Regression_regularization_parameter, 
							epsilon=Epsilon_Support_Vector_Regression_epsilon,
							verbose=verbose)
#regressors = GradientBoostingRegressor(n_estimators=20, random_state=seed)

#note: Unable to obtain reference illumina reads for HG01123, HG02109, HG02486, HG02559. Other 40 reference sequences were used in the panel dataset.

statistics_inherent_to_variant = [
    'ref_graph_minor_allele_freq',
	'ref_graph_alt_alleles',
	'ref_graph_total_alleles',
	'ref_graph_missing',
    'unique_kmers',
    'number_of_alt_alleles',
	'number_of_called_alleles'
	# 'len_ref',
	# 'len_alt'
]
statistics_describing_called_variants = [
    'GQ>=200',
    # 'mean_log_GQ',
    'fraction_of_samples_that_are_heterozygous',
	'number_of_heterozygous_samples',
	'number_of_called_samples',
	'percent_calls_missing',
	'number_calls_missing'
]
statistics_requiring_ground_truth=[
	'percent_calls_correct',
	'percent_calls_incorrect',
	'number_calls_correct',
	'number_calls_incorrect',
    'percent_ref_called_ref',
	'percent_ref_called_het',
	'percent_ref_called_alt',
	'percent_ref_uncalled',
	'percent_het_called_ref',
	'percent_het_called_het',
	'percent_het_called_alt',
	'percent_het_uncalled',
	'percent_alt_called_ref',
	'percent_alt_called_het',
	'percent_alt_called_alt',
	'percent_alt_uncalled']

call_accuracy_statistics = statistics_describing_called_variants + statistics_inherent_to_variant

regression_features = statistics_requiring_ground_truth + call_accuracy_statistics

def save_model(model_object, save_location):
	with open(save_location, 'wb') as file:
		pkl.dump(model_object, file)


def load_model(save_location):
	with open(save_location, 'rb') as file:
		return pkl.load(file)


def postprocess_vcf(genotyped_vcf_file, panel_vcf_file=None, sizes_to_filter=None, trios=None, save_model_loc=None, save_models=False, load_presaved_models=False):
	# pangenie produces pangenie-{sample}_path_segments.fasta and pangenie-{sample}_genotyping.vcf
	# Include option to only process SVs
	# Otherwise separate out all ref samples from pangenie-{sample}_genotyping.vcf and merge resulting vcfs
	# together into one pangenie-unknown_genotyping.vcf which is run here
	# panel_vcf file is the reference vcf
	
	variant_df, per_sample_arrays = load_variant_df(genotyped_vcf_file, sizes_to_filter)
	panel_df, panel_arrays = load_variant_df(panel_vcf_file, sizes_to_filter)

	if panel_vcf_file:
		variant_df = calculate_call_accuracy_statistics(variant_df, per_sample_arrays, ground_truth_df=panel_df, ground_truth_arrays=panel_arrays)
	else: 
		variant_df = fill_in_dummy_statistics_requiring_ground_truth(variant_df)
	if trios:
		variant_df = calculate_mendelian_statistics(variant_df, per_sample_arrays, trios)
	
	
	variant_df = calculate_variant_statistics(variant_df, per_sample_arrays)
	variant_df = calculate_regression_filters(variant_df, use_trios=trios)
	variant_df = run_regression_scoring(variant_df, save_model_loc, skip_SNPs=True, save_models=save_models, load_presaved_models=load_presaved_models)
	return variant_df


def classify_variants(variant_df):
	# must be biallelic - otherwise classification of indel size is difficult
	snps = (variant_df.ref_size == 1) & (variant_df.alt_size == 1)#set(id for id in df.variant_id if 'SNV' in id)
	indel_size = variant_df.variant_size.abs()
	svs = ~snps & (indel_size >= SV_min_size)
	small_indels = ~snps & ~svs & (indel_size < small_indel_max_size) #set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])<=19))
	midsize_indels = ~snps & ~svs & ~small_indels #set(id for id in df.variant_id if (not 'SNV' in id) and (int(id.split('-')[-1])>=small_indel_max_size) and (int(id.split('-')[-1])<SV_min_size))
	
	insertions = (variant_df.variant_size > 0) & (variant_df.alt_size == 1)#variant_id.str.split('-')[2] == "INS"
	deletions = (variant_df.variant_size < 0) & (variant_df.ref_size == 1)
	complex_indels = (variant_df.ref_size > 1) & (variant_df.alt_size > 1)

	small_insertions = small_indels & insertions
	small_deletions = small_indels & deletions
	small_complex = small_indels & complex_indels
	assert sum(small_insertions) + sum(small_deletions) + sum(small_complex) == sum(small_indels)

	midsize_insertions = midsize_indels & insertions
	midsize_deletions = midsize_indels & deletions
	midsize_complex = midsize_indels & complex_indels
	assert sum(midsize_insertions) + sum(midsize_deletions) + sum(midsize_complex) == sum(midsize_indels)

	large_insertions = svs & insertions
	large_deletions = svs & deletions
	large_complex = svs & complex_indels
	assert sum(large_insertions) + sum(large_deletions) + sum(large_complex) == sum(svs)

	# def filter_by_size(variant_df, SNPs=None, small_indels=None, medium_indels=None, SVs=None):
	# 	if 

	variant_df['mutation_category'] = ''
	variant_df['mutation_size'] = ''
	variant_df.loc[snps, 'mutation_category'] = 'SNP'
	variant_df.loc[insertions, 'mutation_category'] = 'insertion'
	variant_df.loc[deletions, 'mutation_category'] = 'deletion'
	variant_df.loc[complex_indels, 'mutation_category'] = 'complex_indel'
	variant_df.loc[snps, 'mutation_size'] = 'SNP'
	variant_df.loc[svs, 'mutation_size'] = 'SV'
	variant_df.loc[small_indels, 'mutation_size'] = 'small_indel'
	variant_df.loc[midsize_indels, 'mutation_size'] = 'midsize_indel'

	return variant_df


def vcf_to_dataframe(path, reader_threads=1):
	'''note: assumes biallelic'''
	variants = tuple(x for x in VCF(path, gts012=True, threads=reader_threads))
	print (f"Extracting data from variants..")
	samples = VCF(path).samples
	vcf_columns = ['chrom','pos','ref','alt','ID','unique_kmers', 'number_of_called_samples', 'number_of_heterozygous_samples', 'minor_allele_frequency']
	if 'UK' in [x[0] for x in variants[0].INFO]:
		vcf_columns = ['chrom','pos','ref','alt','ID','unique_kmers', 'number_of_called_samples', 'number_of_heterozygous_samples', 'minor_allele_frequency']
		variant_df = pd.DataFrame([(var.CHROM, var.POS, var.REF, var.ALT[0], var.INFO['ID'], np.int32(var.INFO['UK']), np.int32(var.num_called), np.int32(var.num_het), np.float16(var.aaf))
						for var in variants], columns=vcf_columns)
	else:
		vcf_columns = ['chrom','pos','ref','alt','ID', 'number_of_called_samples', 'number_of_heterozygous_samples', 'minor_allele_frequency']
		variant_df = pd.DataFrame([(var.CHROM, var.POS, var.REF, var.ALT[0], var.INFO['ID'], np.int32(var.num_called), np.int32(var.num_het), np.float16(var.aaf))
						for var in variants], columns=vcf_columns)
	gts = np.fromiter((x.gt_types[0:len(samples)] for x in variants), dtype=((np.int8, len(samples))))
	quals = np.fromiter((x.gt_quals[0:len(samples)] for x in variants), dtype=((np.int16, len(samples))))
	
	gts = gts.reshape(len(variant_df), len(samples))
	quals = quals.reshape(len(variant_df), len(samples))
	variant_info = {'genotypes':gts, 'GQ':quals, 'samples':samples}
	return variant_df, variant_info


def load_variant_df(path, sizes_to_filter=None, save_as_parquet=True, threads=8):
	print (sizes_to_filter)
	if os.path.exists(path+'.parquet') and os.path.exists(path+'.gt.parquet') and os.path.exists(path+'.quals.parquet'):
		print (f"Reading {path.split('/')[-1]}.parquet into memory...")
		variant_df = pd.read_parquet(path+'.parquet').reset_index(drop=True)
		variant_info=dict()
		variant_info['genotypes']=pd.read_parquet(path+'.gt.parquet').reset_index(drop=True)
		variant_info['GQ'] = pd.read_parquet(path+'.quals.parquet').reset_index(drop=True)
		variant_info['samples'] = variant_info['genotypes'].columns
		save_as_parquet=False
	else:
		print (f"Reading {path.split('/')[-1]} into memory...")
		variant_df, variant_info = vcf_to_dataframe(path, reader_threads=threads)
		variant_info['genotypes'] = pd.DataFrame(variant_info['genotypes'], columns=variant_info['samples']).replace(3, np.nan)
		variant_info['GQ'] = pd.DataFrame(variant_info['GQ'], columns=variant_info['samples']).replace(-1, np.nan)
		variant_df = variant_df.join(variant_info['genotypes'])
		
		variant_df['ref_size'] = variant_df.ref.str.len()
		variant_df['alt_size'] = variant_df.alt.str.len()
		variant_df['variant_size'] = variant_df.ref_size - variant_df.alt_size
	
	variant_df = classify_variants(variant_df)

	if sizes_to_filter:
		for key, v in sizes_to_filter.items():
			if not v:
				variant_df = variant_df.loc[variant_df.mutation_size!=key]
		
		# variant_df = variant_df.loc[(variant_df.ref_size > small_indel_max_size)|(variant_df.alt_size > small_indel_max_size)]
		variant_info['genotypes'] = variant_info['genotypes'].iloc[variant_df.index.to_numpy()]
		variant_info['GQ'] = variant_info['GQ'].iloc[variant_df.index.to_numpy()]
	
	# print (len(variant_df.loc[variant_df.alt.str.contains(',')]))

	# print (variant_df.shape)
	# print (variant_info['genotypes'].shape)
	# print (variant_info['GQ'].shape)
	# if not SNPs:
	# 	variant_df = variant_df.loc[variant_df.variant_size.abs() < 50]
	# if not small_indels:
	# 	variant_df = variant_df.loc[(variant_df.variant_size.abs() >= 50) & (variant_df.variant_size.abs() < 20)]
	# if not medium_indels:
	# 	variant_df = variant_df.loc[variant_df.variant_size.abs() < 50]
	# if not SVs:
	# 	variant_df = variant_df.loc[variant_df.variant_size.abs() < 50]

	if save_as_parquet:
		upcast_for_parquet(variant_df).to_parquet(path + '.parquet')
		upcast_for_parquet(variant_info['GQ']).to_parquet(path+'.gt.parquet')
		upcast_for_parquet(pd.DataFrame(variant_info['GQ'], columns=variant_info['samples'])).to_parquet(path+'.quals.parquet')

	return variant_df, variant_info


def calculate_variant_statistics(variant_df, variant_info):
	gts = variant_df[variant_info['genotypes'].columns].astype(np.float16)
	# print (gts.shape)
	# print (variant_df.shape)
	gts[gts==3] = np.nan
	variant_df['number_of_alt_alleles'] = np.nansum(gts, axis=1)
	variant_df['number_of_called_alleles'] = variant_df.number_of_called_samples*2
	variant_df['fraction_of_samples_that_are_heterozygous'] = variant_df.number_of_heterozygous_samples/variant_df.number_of_called_samples
	variant_df['GQ>=200'] = (variant_info['GQ'] >= 200).sum(axis=1)
	# variant_df['mean_log_GQ'] = np.log(variant_info['GQ']).mean(axis=1)
	return variant_df
	

def fill_in_dummy_statistics_requiring_ground_truth(called_variant_df):
	called_variant_df[call_accuracy_statistics] = np.float16(np.nan)
	return called_variant_df


def calculate_call_accuracy_statistics(called_variant_df, called_variant_arrays, ground_truth_vcf=None, ground_truth_df=None, ground_truth_arrays=None):
	print (called_variant_df.shape)
	if ground_truth_vcf:
		ground_truth_df, ground_truth_arrays = load_variant_df(ground_truth_vcf)
	# print (ground_truth_df.shape)
	print (f"Calculating genotype calling statistics...")
	ground_truth_samples = list(set(ground_truth_arrays['samples']).intersection(set(called_variant_df.columns)))
	# cv_vars = called_variant_df.chrom + called_variant_df.pos.astype(str)+called_variant_df.ref+called_variant_df.alt
	# gt_vars = ground_truth_df.chrom + ground_truth_df.pos.astype(str)+ground_truth_df.ref+ground_truth_df.alt
	called_variant_df = called_variant_df.merge(ground_truth_df[['chrom','pos', 'ref','alt']], on=['chrom','pos', 'ref', 'alt'], how='inner')
	ground_truth_df = ground_truth_df.merge(called_variant_df[['chrom','pos', 'ref','alt']], on=['chrom','pos', 'ref', 'alt'], how='inner')
	
	
	# called_gt = called_variant_df[called_variant_arrays['samples']].values
	reference_gt = ground_truth_df[ground_truth_arrays['samples']].values

	nan_ref_gts = np.isnan(reference_gt).astype(bool)
	
	print ('Calculating reference stats...')
	# called_variant_df['len_ref'] = called_variant_df.ref.str.len()
	# called_variant_df['len_alt'] = called_variant_df.alt.str.len()
	called_variant_df['ref_graph_missing'] = nan_ref_gts.sum(1).astype(np.int32)
	called_variant_df['ref_graph_minor_allele_freq'] = np.nanmean(reference_gt, axis=1).astype(np.int32)
	called_variant_df['ref_graph_alt_alleles'] = np.nansum(reference_gt, axis=1).astype(np.int32)
	called_variant_df['ref_graph_total_alleles'] = ((~nan_ref_gts).sum(1)*2).astype(np.int32)
	# print (called_variant_df.shape)
	if len(ground_truth_samples) == 0:
		called_variant_df = fill_in_dummy_statistics_requiring_ground_truth(called_variant_df)
	else:
		print ('Calculating how many calls were correct...')
		common_samples_called = called_variant_df[ground_truth_samples].values
		common_samples_ref = ground_truth_df[ground_truth_samples].values
		basecall_accurate_array = (common_samples_called == common_samples_ref).astype(bool)
	
		nan_accuracy = np.isnan(basecall_accurate_array).astype(bool)
		nan_ref_gts = np.isnan(common_samples_ref).astype(bool)
		nan_called_gts = np.isnan(common_samples_called).astype(bool)
		# nan_accuracy = nan_ref_gts | nan_called_gts
		ref_zeros = (common_samples_ref==0).astype(bool)
		ref_ones = (common_samples_ref==1).astype(bool)
		ref_twos = (common_samples_ref==2).astype(bool)
		called_zeros = (common_samples_called==0).astype(bool)
		called_ones = (common_samples_called==1).astype(bool)
		called_twos = (common_samples_called==2).astype(bool)
		# with warnings.catch_warnings(): # numpy warns if performing nanmean on an empty 
		# warnings.simplefilter("ignore", category=RuntimeWarning)

		# print (called_variant_df.shape)
		called_variant_df['percent_calls_missing'] = nan_accuracy.mean(axis=1).astype(np.float16)
		called_variant_df['percent_calls_correct'] = np.nanmean(basecall_accurate_array, axis=1).astype(np.float16)
		called_variant_df['percent_calls_incorrect'] = (1 - called_variant_df.percent_calls_correct - called_variant_df.percent_calls_missing ).astype(np.float16)
		called_variant_df['number_calls_missing'] = nan_accuracy.sum(axis=1).astype(np.int32)
		called_variant_df['number_calls_correct'] = np.nansum(basecall_accurate_array, axis=1).astype(np.int32)
		called_variant_df['number_calls_incorrect'] = (basecall_accurate_array.shape[1] - called_variant_df.number_calls_correct - called_variant_df.number_calls_missing).astype(np.int32)
		called_variant_df['percent_ref_called_ref'] = (np.nansum((ref_zeros & called_zeros), axis=1)/np.nansum(ref_zeros, axis=1)).astype(np.float16)
		called_variant_df['percent_ref_called_het'] = (np.nansum((ref_zeros & called_ones), axis=1)/np.nansum(ref_zeros, axis=1)).astype(np.float16)
		called_variant_df['percent_ref_called_alt'] = (np.nansum((ref_zeros & called_twos), axis=1)/np.nansum(ref_zeros, axis=1)).astype(np.float16)
		called_variant_df['percent_ref_uncalled']   = (np.nansum((ref_zeros & nan_called_gts), axis=1)/np.nansum(ref_zeros, axis=1)).astype(np.float16)
		called_variant_df['percent_het_called_ref'] = (np.nansum((ref_ones & called_zeros), axis=1)/np.nansum(ref_ones, axis=1)).astype(np.float16)
		called_variant_df['percent_het_called_het'] = (np.nansum((ref_ones & called_ones), axis=1)/np.nansum(ref_ones, axis=1)).astype(np.float16)
		called_variant_df['percent_het_called_alt'] = (np.nansum((ref_ones & called_twos), axis=1)/np.nansum(ref_ones, axis=1)).astype(np.float16)
		called_variant_df['percent_het_uncalled']   = (np.nansum((ref_ones & nan_called_gts), axis=1)/np.nansum(ref_ones, axis=1)).astype(np.float16)
		called_variant_df['percent_alt_called_ref'] = (np.nansum((ref_twos & called_zeros), axis=1)/np.nansum(ref_twos, axis=1)).astype(np.float16)
		called_variant_df['percent_alt_called_het'] = (np.nansum((ref_twos & called_ones), axis=1)/np.nansum(ref_twos, axis=1)).astype(np.float16)
		called_variant_df['percent_alt_called_alt'] = (np.nansum((ref_twos & called_twos), axis=1)/np.nansum(ref_twos, axis=1)).astype(np.float16)
		called_variant_df['percent_alt_uncalled']   = (np.nansum((ref_twos & nan_called_gts), axis=1)/np.nansum(ref_twos, axis=1)).astype(np.float16)
	# print (called_variant_df.shape)
	return called_variant_df


def calculate_mendelian_statistics(variant_df, per_sample_arrays, trios):
	raise NotImplementedError
	return variant_df


def calculate_regression_filters(variant_df, use_trios=False):
	if use_trios:
		mendel_fail = (variant_df.pangenie_mendelian_consistency < 0.8) & (variant_df['pangenie_considered_trios']>=5)
	else:
		mendel_fail = np.zeros(len(variant_df), dtype=bool)
	ac0_fail = variant_df['minor_allele_frequency'] == 0
	variant_df['ac0_fail'] = ac0_fail
	gq_pass = variant_df['GQ>=200'] > variant_df['GQ>=200'].max()*0.9
	gq_fail = variant_df['GQ>=200'] < variant_df['GQ>=200'].max()*0.1
	self_fail = variant_df['percent_calls_correct'] < 0.9
	nonref_fail = ((variant_df[['percent_het_called_het',
								'percent_alt_called_alt']] == 0).all(axis=1) & 
				   (variant_df[['percent_het_called_ref', 
								'percent_het_called_alt',
								'percent_alt_called_ref',
								'percent_alt_called_het']] != 0).any(axis=1))
	likely_true = (~(variant_df.ac0_fail | mendel_fail | gq_pass | nonref_fail | self_fail)).astype(bool)
	likely_false = (~variant_df.ac0_fail & ((gq_fail.astype(np.int8) + mendel_fail.astype(np.int8) + nonref_fail.astype(np.int8) + self_fail.astype(np.int8)) >= 2)).astype(bool)
	
	variant_df['training_set_category'] = 0
	variant_df.loc[likely_true, 'training_set_category'] = 1
	variant_df.loc[likely_false, 'training_set_category'] = -1
	return variant_df

def run_regression_scoring(all_variants_df, save_model_loc=None, skip_SNPs=True, save_models=False, load_presaved_models=False):
	# load models
	'''high memory - seeing 100+GB when regressing against 6 samples' SNPs'''
	imputors, scalers, regressors = {}, {}, {}
	
	if load_presaved_models:
		# assert os.path.exists(save_model_loc+'.imputors.pkl'), 'No impution model at save location!'
		# assert os.path.exists(save_model_loc+'.scalers.pkl'), 'No scaler model at save location!'
		# assert os.path.exists(save_model_loc+'.regressors.pkl'), 'No regressor model at save location!'
		assert os.path.exists(save_model_loc+'.pkl'), 'No model at save location!'
		print (f"Loading models previously saved at {save_model_loc}.*.pkl")
		imputors, scalers, regressors = load_model(save_model_loc+'.pkl')
		# scalers = load_model(save_model_loc+'.scalers.pkl')
		# regressors = load_model(save_model_loc+'.regressors.pkl')
	
	# drop ac0_fails
	all_variants_df = all_variants_df.loc[~all_variants_df.ac0_fail]
	
	all_variants_df['regression_score'] = np.nan
	variant_groups_to_regress = all_variants_df.groupby(['mutation_category','mutation_size'])

	for (mutation_category, mutation_size), variant_df in variant_groups_to_regress:
		if skip_SNPs and (mutation_size in ('SNP', 'small_indel')):
			continue
		print (f'Calculating genotype confidence for {mutation_size} {mutation_category}')
		# variant_df = variant_df.sort_values('variant_id').set_index('variant_id')
		
		# impute missing values
		cat_key = mutation_category+'.'+mutation_size
		if mutation_category+'.'+mutation_size not in imputors.keys():
			imputors[cat_key] = deepcopy(imputor)
			imputors[cat_key].fit(variant_df[statistics_requiring_ground_truth+call_accuracy_statistics].values)
		else:
			# if using pre-saved model, then evaluating that model 
			# requires imputing the statistics we would not have in real samples with
			# unknown ground truth genotypes. I can either not incorporate that info into the model,
			# or I can try to impute that info from the statistics we do have,
			# which is likely the better method. (Though should probably test that!)
			variant_df.loc[:, call_accuracy_statistics] = np.nan 
		variant_df.loc[:, regression_features] = imputors[cat_key].transform(variant_df[regression_features].values)
		
		# Scale all regession points to -1,0,1
		if cat_key not in scalers.keys():
			scalers[cat_key] = deepcopy(scaler)
			scalers[cat_key].fit(variant_df[regression_features].values)
		variant_df.loc[:,regression_features] = scalers[cat_key].transform(variant_df.loc[:,regression_features].values)
		
		# Train model only on labeled points
		df_labeled = variant_df.loc[(variant_df.training_set_category != 0)]
		training_set_category = df_labeled.loc[:, 'training_set_category'].values.ravel()
		
		if cat_key not in regressors.keys():
			regressors[cat_key] = deepcopy(regressor)
			print('Training regression model')
			regressors[cat_key].fit(df_labeled.loc[:,regression_features].values, training_set_category)
		
		# Apply to all_data
		y_pred = regressors[cat_key].predict(variant_df.loc[:,regression_features].values)
		
		# Add column with variant specific scores to table
		all_variants_df.loc[(all_variants_df.mutation_category==mutation_category) &
							(all_variants_df.mutation_size==mutation_size),
							'regression_score'] = y_pred
	
	all_variants_df['confidence_level'] = 4
	all_variants_df.loc[all_variants_df.regression_score <  0.5, 'confidence_level'] = 3
	all_variants_df.loc[all_variants_df.regression_score <  0.0, 'confidence_level'] = 2
	all_variants_df.loc[all_variants_df.regression_score < -0.5, 'confidence_level'] = 1
	all_variants_df.loc[all_variants_df.training_set_category == 1, 'confidence_level'] = 4
	
	if save_models:
		print (f'Saving regression model to {save_model_loc}')
		save_model((imputors, scalers, regressors), save_model_loc+'.pkl')
		# save_model(scalers, save_model_loc+'.scalers.pkl')
		# save_model(regressors, save_model_loc+'.regressors.pkl')
	return all_variants_df


def save_passing_variant_csv(variant_df, outfile_prefix, regression_score_cutoff):
	passing_variants = variant_df.loc[(variant_df.regression_score > regression_score_cutoff) | 
									  (variant_df.training_set_category == 1)]
	passing_variants.ID.to_csv(outfile_prefix + '.csv', index=False, header=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar="\\")


def upcast_for_parquet(df):
	# parquet is a very space-efficient method of storing data,
	# but it does not handle np.float16. 
	# So, all float16s must be upcast to float32s.
	for i, dtype in enumerate(df.dtypes):
		if dtype==np.float16:
			df[df.columns[i]] = df[df.columns[i]].astype(np.float32)
	return df


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--genotype_vcf', type=str,
						help='Combined VCF of variant calls made with pangenie')
	parser.add_argument('-p', '--panel_vcf', type=str,
						help='VCF of ground truth variant calls used to make reference graph')
	parser.add_argument('-o', '--output', type=str,
						help='Save prefix of parquet results file and csv of passing filters')
	parser.add_argument('-c', '--cutoff', type=float, default=0.75,
						help='Cutoff score for what is considered to be a hit. Range -1 to 1.')
	parser.add_argument('--SNPs', action='store_true',
						help='analyze SNPs.')
	parser.add_argument('--small_indels', action='store_true',
						help='analyze small indels (< 20bp).')
	parser.add_argument('--medium_indels', action='store_true',
						help='analyze moderately sized indels (indels 20-50 bp)')
	parser.add_argument('--SVs', action='store_true',
						help='analyze SVs (indels above 50 bp)')
	parser.add_argument('-t', '--trios', type=str, default=None,
						help='.ped file of trios to analyze.')
	parser.add_argument('--model', type=str, default=None,
						help='Location to load or save regression model parameters.')
	parser.add_argument('--save_model', action='store_true',
						help='Save any generated models.')
	parser.add_argument('--load_presaved_model', action='store_true',
						help='Load model at model.')
	parser.add_argument('--version', action='store_true',
						help='Print version and quit.')

	args = parser.parse_args()
	if args.version:
		print (version)
	else:
		sizes_to_filter={'SNP':args.SNPs, 'small_indel':args.small_indels, 'medium_indel':args.medium_indels, 'SV':args.SVs}
		print (sizes_to_filter)
		variant_df = postprocess_vcf(args.genotype_vcf, panel_vcf_file=args.panel_vcf, sizes_to_filter=sizes_to_filter, trios=args.trios, save_model_loc=args.model, save_models=args.save_models, load_presaved_models=args.load_presaved_models)
		
		save_passing_variant_csv(variant_df, args.output, args.cutoff)
		
		upcast_for_parquet(variant_df).to_parquet(args.output + '.parquet')
	



# genotyped_vcf_file="/mnt/data/lalli/nf_stage/sarek_JLL/local_references/pangenie_og_test/test.samples.merged.bcf.gz"
# panel_vcf_file="/mnt/data/lalli/nf_stage/sarek_JLL/local_references/pangenie_og_test/panel.vcf"
# output="/mnt/data/lalli/nf_stage/sarek_JLL/local_references/pangenie_og_test/test_postprocess_output"
# SV_only=False
# trios=None
# save_model_loc=None
	# I think what I want is, for different cutoffs of regression_score, how do the above values vary
	# Can plot mulitple lines, each line different category of variant or different ancestry background
	# Longer term, it would be cool to plot the same thing for different methods of basecalling/different references
	# for (mutation_category, mutation_size), category_df in variant_df.groupby(['mutation_category','mutation_size']):
	# 	category_indexes = category_df.index
	# 	cat_precisions = var_precision[category_indexes]
	# 	cat_recalls = var_recall[category_indexes]
	# 	cat_gene_concordances = var_gene_concordance[category_indexes]
	# 	cat_Fscores = var_F_score[category_indexes]
	# 	modeled_confidence = category_df.regression_score
