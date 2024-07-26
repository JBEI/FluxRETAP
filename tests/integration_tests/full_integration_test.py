import cobra
from cobra import Reaction, Metabolite
import numpy as np
import pandas as pd
import os
import pytest

from core.FluxRETAP import FluxRETAP

def compareDFs(x,frozen_x):
	'''
	Auxiliary function
	Compares x a data frame with flux RETAP recommendations, to the version of x frozen in disk, frozen_x
	'''

	# Limit the comparison to only the amount found in the frozen recommendations
	L = frozen_x.shape[0]
	x = x[0:L].copy()

	## Comparison with frozen data
	# that the reactions returned are in the same order. minor stochastic changes can lead to changes in position of the overall dataframe, but the set returned should be in order. 
	assert set(frozen_x['index']) == set(x['index'])

	#assert set(frozen_x.index) == set(x.index)
	#x = x.reindex(frozen_x.index)

	# We now separate dataframe numerical and string columns so they can be separately compared (need tolerance for numerical part) 
	frozen_x_Number = frozen_x.select_dtypes(include=np.number).to_numpy()
	frozen_x_String = frozen_x.select_dtypes(include=[np.dtype,'bool'])

	x_Number = x.select_dtypes(include=np.number).to_numpy()
	x_String = x.select_dtypes(include=[np.dtype,'bool'])

	# reindex the string for comparisons. 
	x_String.index = frozen_x_String.index

	## Ensure all string columns are the same
	try:
		assert frozen_x_String.eq(x_String).all().all()
	# instances of where the minor changes reorder the returned reaction order but do not change the values
	except:
		x_String_reordered = x_String.set_index('index').loc[frozen_x_String['index']].reset_index()
		x_String_reordered.index = frozen_x_String.index
		x_String_reordered['geneName']=x_String_reordered.geneName.astype(str)

		assert frozen_x_String.eq(x_String_reordered).all().all()

	# Ensure numeric columns are approximately the same (allow for 0.1% tolerance)
	assert np.isclose(frozen_x_Number,x_Number,rtol=1e-03).all()


def test_Tian_2019():
	'''
	Integration test based on data by Tian, Tian, et al. "Redirecting metabolic flux via combinatorial multiplex 
	CRISPRi-mediated repression for isopentenol production in Escherichia coli." 
	ACS synthetic biology 8.2 (2019): 391-402.
	'''

	### Ground truth data for comparison
	# These are the genes and reactions that worked in the paper (i.e. increased production)
	worked_genes = ['gldA','asnA','thrC','maeB','acs','poxB','ackA','pta']
	worked_reactions = ['PTA2','GLYCDx','ACCOAL','PTAr','APPLDHr','LALDO2x','THRS','ME2','ACKr','POX','4HTHRS',
						'ACS','ASNS2','ALR4x']

	## Load the E. coli Genome-scale Model (GSM)
	# Using E. coli GSM iJO1366 with isoprenol production integrated
	libraryPath = "../.." 
	file_name = libraryPath+'/models/iJO1366_MVA.json' 
	model = cobra.io.load_json_model(file_name)

	# Choose metabolic subsystems in central carbon metabolism, as done in the paper
	subsystems = ['Glycolysis/Gluconeogenesis','Citric Acid Cycle','Pyruvate Metabolism','Pentose Phosphate Pathway',
				  'Threonine and Lysine Metabolism','Alternate Carbon Metabolism','Pyruvate Metabolism',
				  'Glycerophospholipid Metabolism','Methylglyoxal Metabolism','Alanine and Aspartate Metabolism',
				  'Oxidative Phosphorylation','Anaplerotic Reactions']

	### Flux RETAP steps
	# Flux RETAP input
	biomassRxn = 'BIOMASS_Ec_iJO1366_core_53p95M'
	uptakeRxn = 'EX_glc__D_e'
	UpOrDown = 'Down'
	productRxn = 'EX_isoprenol_e'

	# Run Flux RETAP
	isoprenol = FluxRETAP.getRecommendations(model=model,
								   UpOrDown=UpOrDown,
								   productRxn=productRxn,
								   carbonSourceRxn=uptakeRxn,
								   biomassRxn=biomassRxn,
								   desiredSystems=subsystems,
								   N=10, 
								   referenceCutOff=440,
								   optimalFraction=0.0,
								   fluxRangeDiff= 0.0005,
								   minGrowth = -0.01,
								   Ors = True)

	### Compare results with previous run (frozen data)
	## Process data to be comparable with frozen data
	# Add to recommendations dataframe the reactions that worked in the paper and limit to those
	isoprenol2 = isoprenol.reset_index()
	worked_inds = isoprenol2['index']==worked_reactions[0]
	for rxn in worked_reactions:
		worked_inds = worked_inds | (isoprenol2['index']==rxn)
	isoprenolFinal = isoprenol2.loc[worked_inds]

	# Let's write the new recommendations to a csv and then reimport to avoid format issues
	dummyFileName = './dummyFile.csv'
	isoprenolFinal.to_csv(dummyFileName)
	selected = pd.read_csv(dummyFileName,index_col=0,keep_default_na=False)
	os.remove(dummyFileName) # Erase dummy file

	## Get frozen data
	recommendationsFilename = libraryPath+'/tests/integration_tests/files/Tian_2019_recommendations.csv' 
	frozenSelected = pd.read_csv(recommendationsFilename,index_col=0,keep_default_na=False)	 

	# Compare dataframes
	compareDFs(selected,frozenSelected)

def test_Boghigian_2012():
	'''
	Integration test based on data by Boghigian, Brett A., et al. "Computational identification of gene 
	over-expression targets for metabolic engineering of taxadiene production." Applied microbiology and
	biotechnology 93 (2012): 2063-2073.
	'''

	### Ground truth data for comparison
	# Gene targes identified for overexpression by algorithm in paper
	paper_gene_targets_identified = ['ppk', 'sthA', 'purN', 'dxs', 'ispE', 'dxr','ispG','ispF','ispD','folD','ispH','ispA'] 

	# Gens that experimentally worked
	experimentally_worked = ['idi','ppk','sthA']

	## Create the E. coli Genome-scale Model (GSM)
	# Load E. coli GSM iJO1366 with isoprenol production integrated
	libraryPath = "../.." 
	file_name = libraryPath+'/models/iAF1260.json' 
	model = cobra.io.load_json_model(file_name)

	# Add heterologous reactions
	# geranylgeranyl diphosphate synthase
	GGPS = cobra.Reaction('GGPS')
	GGPS.name = 'GGPS'
	GGPS.subsystem = 'Cytosol'
	GGPS.lower_bound = 0
	GGPS.upper_bound = 1000
	model.add_reactions([GGPS])	 # add rxn to the model
	# Create non-native metabolites and add to reaction
	ggpp_c = cobra.Metabolite( 'ggpp_c', formula ='', name = 'GGPP', compartment='c')
	GGPS.add_metabolites({'frdp_c': -1, 'ipdp_c': -1, 'ppi_c': 1, ggpp_c: 1})

	#Taxadiene production 
	taxa = cobra.Reaction('taxa')
	taxa.name = 'taxa'
	taxa.subsystem = 'Cytosol'
	taxa.lower_bound = 0
	taxa.upper_bound = 1000
	model.add_reactions([taxa])
	# Create non-native metabolites and add to reaction
	taxadiene_c = cobra.Metabolite('taxadiene_c', formula ='', name = 'Taxadiene', compartment='c')
	taxa.add_metabolites({ 'ggpp_c': -1, 'ppi_c': 1, taxadiene_c: 1})
	# Make taxadiene a demand reaction so it can accumulate inside the cell:
	demand = model.add_boundary(model.metabolites.taxadiene_c,type="demand")

	# Set model to experimentally replicate theirs when growing on glycerol
	# (growth should now equal 0.2671/h)
	model.reactions.EX_glc__D_e.bounds=0,0
	model.reactions.EX_glyc_e.bounds = -3,1000

	amino_acids = ['EX_ala__L_e',
				   'EX_gly_e',
				   'EX_val__L_e',
				   'EX_leu__L_e',
				   'EX_ile__L_e',
				   'EX_pro__L_e',
				   'EX_his__L_e',
				   'EX_thr__L_e',
				   'EX_ser__L_e',
				   'EX_phe__L_e',
				   'EX_met__L_e',
				   'EX_glu__L_e',
				   'EX_gln__L_e',
				   'EX_asp__L_e',
				   'EX_asn__L_e',
				   'EX_arg__L_e',
				   'EX_cys__L_e',
				   'EX_trp__L_e',
				   'EX_tyr__L_e',
				   'EX_lys__L_e']
	for i in amino_acids:
		model.reactions.get_by_id(i).lower_bound = -0.1

	### Flux RETAP steps
	# Flux RETAP input
	uptakeRxn = 'EX_glyc_e'
	UpOrDown = 'Up'
	productRxn = 'taxa'
	optimalFraction = 0.0
	N=10
	fluxRangeDiff = 0
	biomassRxn = 'BIOMASS_Ec_iAF1260_core_59p81M'
	referenceCutOff = 0
	# Let's exclude a few subsystems that are not relevant:
	subsystems = [
		'~ ','~Intracellular demand','~Extracellular exchange','~Biomass and maintenance functions',
		'~Intracellular source/sink'
	]

	# Run Flux RETAP
	taxadiene = FluxRETAP.getRecommendations(model=model,
								   UpOrDown=UpOrDown,
								   productRxn=productRxn,
								   carbonSourceRxn=uptakeRxn,
								   biomassRxn=biomassRxn,
								   desiredSystems=subsystems,
								   N=N, 
								   referenceCutOff=referenceCutOff,
								   optimalFraction=optimalFraction,
								   fluxRangeDiff=fluxRangeDiff)

	# Round scores to avoid numerical noise
	taxadiene["score"] = round(taxadiene["score"])
	taxadiene["fluxDiff"] = round(taxadiene["fluxDiff"],5)
	taxadiene["FinalFlux"] = round(taxadiene["FinalFlux"],5)
	taxadiene["growth %"] = round(taxadiene["growth %"],5)

	# Find the reactions corresponding to the genes identified and the ones that 
	# experimentally worked
	worked_reactions	 =set([])
	identified_reactions =set([])

	for gene in paper_gene_targets_identified:
		for rxn in model.reactions:
			if gene in rxn.gene_name_reaction_rule:
				identified_reactions.add(rxn.id)
		
	for gene in experimentally_worked:
		for rxn in model.reactions:
			if gene in rxn.gene_name_reaction_rule:
				worked_reactions.add(rxn.id)
		
	worked_reactions = list(worked_reactions)
	identified_reactions = list(identified_reactions)

	### Compare results with previous run (frozen data)
	## Process data to be comparable with frozen data
	# Add to recommendations dataframe the reactions that worked in the paper 
	# and limit dataframe to those
	taxadiene2 = taxadiene.reset_index()
	taxadiene2 = taxadiene2.sort_values(by=['score', 'index'])
	taxadiene2 = taxadiene2.reset_index(drop=True)
	identified_inds = taxadiene2['index']==identified_reactions[0]
	for rxn in identified_reactions:
		identified_inds = identified_inds | (taxadiene2['index']==rxn)
	taxadineIdentified = taxadiene2.loc[identified_inds]
	taxadineIdentified.reset_index(drop=True)

	taxadiene2 = taxadiene.reset_index()
	worked_inds = taxadiene2['index']==worked_reactions[0]
	for rxn in worked_reactions:
		worked_inds = worked_inds | (taxadiene2['index']==rxn)
	taxadineWorked = taxadiene2.loc[worked_inds]

	# Let's write the new recommendations to a csv and then reimport to avoid format issues
	dummyFileName = './dummyFile.csv'
	taxadineIdentified.to_csv(dummyFileName)
	taxadineIdentified = pd.read_csv(dummyFileName,index_col=0,keep_default_na=False)
	os.remove(dummyFileName) # Erase dummy file
	taxadineWorked.to_csv(dummyFileName)
	taxadineWorked = pd.read_csv(dummyFileName,index_col=0,keep_default_na=False)
	os.remove(dummyFileName) # Erase dummy file

	## Get frozen data
	IdFilename = libraryPath+'/tests/integration_tests/files/Boghigian_2012_identified.csv' 
	WorkedFilename = libraryPath+'/tests/integration_tests/files/Boghigian_2012_worked.csv' 

	frozenId = pd.read_csv(IdFilename,index_col=0,keep_default_na=False)  
	frozenWorked = pd.read_csv(WorkedFilename,index_col=0,keep_default_na=False)  

	# Compare dataframes
	compareDFs(taxadineIdentified,frozenId)
	compareDFs(taxadineWorked,frozenWorked)

	
# Indigoidine Test 
def test_Banerjee_2020():
	'''
	Integration test based on data by Banerjee, Deepanwita, et al. "Genome-scale metabolic 
	rewiring improves titers rates and yields of the non-native product indigoidine at scale." 
	Nature Communications 11 (2020):5385.	
	'''

	### Ground truth data for comparison
	# Gene targets identified for overexpression by algorithm in paper
	paper_targets_identified = ['PP_2082',
			   'PP_1444',
			   'PP_0654',
			   'PP_1251',
			   'PP_5085',
			   'PP_2168',
			   'PP_5176',
			   'PP_0864',
			   'PP_4947',
			   'PP_5003',
			   'PP_0434',
			   'PP_4734',
			   'PP_2925',
			   'PP_5005',
			   'PP_0751',
			   'PP_0100']



	## Create the P. putida Genome-scale Model (GSM)
	# Load P. putida GSM iJN1463 with isoprenol production integrated
	libraryPath = "../.." 
	file_name = libraryPath+'/models/iJN1463_modified.json' 
	model = cobra.io.load_json_model(file_name)


	## add heterologous production steps to the model. 
	APNPT = Reaction('APNPT')
	APNPT.name = 'ATP:pantetheine 4\'-phosphotransferase'
	APNPT.subsystem = 'Cytosol'
	APNPT.lower_bound = -1000
	APNPT.upper_bound = 1000


	# initialize a rxn object and set the bounds
	product_rxn = Reaction('indigoidine_production')
	product_rxn.name = 'Indigoidine production'
	product_rxn.subsystem = 'Cytosol'
	product_rxn.lower_bound = 0
	product_rxn.upper_bound = 1000


	model.add_reactions([product_rxn,APNPT])

	# initialize the non-natvive metabolites 
	ptth_c = Metabolite(
		'ptth_c',
		formula='C11H22N2O4S',
		name='Pantetheine',
		compartment='c'
	)



	indigoidine_c = Metabolite(
		'indigoidine_c',
		formula ='C10H8N4O4',
		name = 'Indigoidine',
		compartment='c'
	)

	# add the stochiometrically balanced reaction to the model
	## note - the names are metabolites names as they appear in the GSM
	product_rxn.add_metabolites({
		'gln__L_c': -2,
		'atp_c': -2,
		'coa_c': -2,
		'fmn_c': -2,
		'o2_c': -2.5,
		'pap_c':2,
		'fmnh2_c': 2,
		'ppi_c': 2,
		'amp_c': 2,
		ptth_c: 2,
		'h2o_c': 1,
		'pi_c': 2,
		indigoidine_c: 1   
	})

	APNPT.add_metabolites({
		'atp_c': -1,
		ptth_c: -1,
		'adp_c': 1,
		'h_c':1,
		'pan4p_c':1
	})


	## allow the product to accumulate in the model 
	# - allows for flux to flow through the reaction in silico
	demand = model.add_boundary(model.metabolites.indigoidine_c,type="demand")


	### Flux RETAP steps
	# Flux RETAP input
	uptakeRxn = 'EX_glc__D_e'
	UpOrDown = 'Both'
	productRxn = 'indigoidine_production'
	optimalFraction = 0.5
	N=10
	fluxRangeDiff = 0
	biomassRxn = 'BIOMASS_KT2440_WT3'
	referenceCutOff = 0
	# Let's exclude a few subsystems that are not relevant:
	subsystems = [
		'~ ','~Intracellular demand','~Extracellular exchange','~Biomass and maintenance functions',
		'~Intracellular source/sink'
	]

	# Run Flux RETAP
	indigoidine = FluxRETAP.getRecommendations(model=model,
								   UpOrDown=UpOrDown,
								   productRxn=productRxn,
								   carbonSourceRxn=uptakeRxn,
								   biomassRxn=biomassRxn,
								   desiredSystems=subsystems,
								   N=N, 
								   referenceCutOff=referenceCutOff,
								   optimalFraction=optimalFraction,
								   fluxRangeDiff=fluxRangeDiff)

	# Round scores to avoid numerical noise
	indigoidine["score"] = round(indigoidine["score"])
	indigoidine["fluxDiff"] = round(indigoidine["fluxDiff"],5)
	indigoidine["FinalFlux"] = round(indigoidine["FinalFlux"],5)
	indigoidine["growth %"] = round(indigoidine["growth %"],5)

	# Find the reactions corresponding to the genes identified and the ones that 
	# experimentally worked
	worked_reactions	 =set([])
	identified_reactions =set([])

	for gene in paper_targets_identified:
		for rxn in model.reactions:
			if gene in rxn.gene_reaction_rule:
				identified_reactions.add(rxn.id)

	identified_reactions = list(identified_reactions)			


	### Compare results with previous run (frozen data)
	## Process data to be comparable with frozen data
	# Add to recommendations dataframe the reactions that worked in the paper 
	# and limit dataframe to those
	indigoidine2 = indigoidine.reset_index()
	indigoidine2 = indigoidine2.sort_values(by=['score', 'index'],ascending=False)
	indigoidine2 = indigoidine2.reset_index(drop=True)

	## Get frozen data
	IdFilename = libraryPath+'/tests/integration_tests/files/Banerjee_2020_identified.csv' 
 
	frozenSelectedInitial = pd.read_csv(IdFilename,index_col=0,keep_default_na=False)	
 
	# Create an empty DataFrame to store the results
	indigoidineIdentified = pd.DataFrame()
	
	# copy of frozenSelected to iterate over
	indigoidineIdentified = indigoidine2[indigoidine2['index'].isin(frozenSelectedInitial['index'])]
	
 
	
	#update frozenSelected	
	frozenSelected = frozenSelectedInitial[frozenSelectedInitial['index'].isin(indigoidineIdentified['index'])]
		
	# Let's write the new recommendations to a csv and then reimport to avoid format issues
	dummyFileName = './dummyFile.csv'
	indigoidineIdentified.to_csv(dummyFileName)
	selected = pd.read_csv(dummyFileName,index_col=0,keep_default_na=False)
	os.remove(dummyFileName) # Erase dummy file
 
	if 'level_0' in selected.columns:
		selected = selected.drop('level_0', axis=1)



	# Compare dataframes
	compareDFs(selected,frozenSelected)

