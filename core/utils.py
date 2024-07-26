from doctest import IGNORE_EXCEPTION_DETAIL
import cobra
from cobra.flux_analysis import flux_variability_analysis
from cobra import Metabolite, Reaction
import numpy as np
from scipy.stats import norm
import pandas as pd
import time
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

pd.options.mode.chained_assignment = None  # default='warn'


def getRanges(model, N, biomassRxn, productRxn, carbonSourceRxn, optimalFraction,fast):
    '''
    Gets flux ranges for different levels of flux forcing through final product reaction

    Args:
        model (cobra.core.model.Model):
            Genome-scale model to perform simulations
        N (int)
            Number of flux forcing steps to simulate
        biomassRxn (str):
            biomass reaction name
        productRxn (str):
            product reaction name
        carbonSourceRxn (str):
            carbon source reaction name
        optimalFraction (float):
            Fraction of biomass to use as minimal biomass constraint
        fast (bool):
        specifies if we want the fast calculation (calculates only the first and last two fractions) 
        or the complete one (calculates all fractions)

    Returns:
        ranges (pd.DataFrame Dictionary):
            FVA values from simulation
    '''

    ## First, we define the step sizes of increasing product we will simulate
    fractions = [x/N for x in range(0,N+1)]
    # Keep only the ones we need for calculations (two first and two last)
    if fast:
        fractions = [fractions[0],fractions[1],fractions[N-1],fractions[N]]  
    # Print flux fractions to be used
    print(f'Fractions of flux: {fractions}')
    

    ## Second, we find the maximum flux possible through the final product reaction
    maxProductFlux = 0
    with model:
        # First set lower bound so biomass does not reach zero
        model.objective = model.reactions.get_by_id(biomassRxn)
        maxBio = model.optimize().objective_value
        model.reactions.get_by_id(biomassRxn).lower_bound = maxBio*optimalFraction
        
        # Now set objective to product reaction and maximize its flux
        model.objective = model.reactions.get_by_id(productRxn)
        solution = model.optimize()
        maxProductFlux = solution.objective_value  # This is the maximum flux allowed through the reaction 
        print(f'Maximum flux to product: {maxProductFlux}') # print maximum product flux

    
    ## Third, force increasing levels of flux through final product pathway
    start_time = time.time()  # to time full simulation
    print('\nstarting FVA simulations . . .')
    ranges = {} # holds the results of the FVA
    
    # Use FVA to find maximum and minimum fluxes compatible flux constraint
    for num, x in enumerate(fractions):
        with model:
            # Constrains flux to fraction x of total possible 
            model.reactions.get_by_id(productRxn).upper_bound = x*maxProductFlux    
            model.reactions.get_by_id(productRxn).lower_bound = x*maxProductFlux
            model.reactions.get_by_id(biomassRxn).lower_bound =  maxBio*optimalFraction
                
            # FBA with constraints to get the 'optimal fit'
            solution = model.optimize()                                                                 

            # FVA with constraints to get the flux ranges
            FVAminmax = flux_variability_analysis(model, model.reactions, fraction_of_optimum=0.999,loopless=False)   
            optimals = [solution.fluxes[i] for i,j in FVAminmax.iterrows()]
            FVAminmax.insert(1, "Optimal", optimals) 

            # Normalize FVA ranges by input flux 
            ranges[num] = FVAminmax/abs(solution.fluxes[carbonSourceRxn])
        
        # Print progress of simulations
        if ((num+1)/len(fractions))%0.25==0:
            print('{} % complete'.format((num+1)/len(fractions)*100))
    
    # Print time the full simulation took
    print("--- %s seconds ---" % (time.time() - start_time))

    return ranges

def getFluxScoresTable(model, ranges, desiredSystems):
    '''
    Produces table with flux scores according to the overlap of distribution of initial and final fluxes
    (initial and final in terms of flux fraction forced into final production reaction)

    Args:
        model (cobra.core.model.Model):
            Genome-scale model to perform simulations
        ranges (pd.DataFrame Dictionary):
            FVA values from simulation
        returnGaussian (bool):
                specifies how much information to return. Default is to return three list and a float. 
        desiredSystems (list):
                subsystems of interest
    Returns:
        fluxScores (pd.DataFrame): 
                A dataframe with fluxes and their corresponding scores
    '''
    
    ## Get gaussian distributions that describe fluxes around the intial and final points of the flux forcing scheme
    initGaussian, finGaussian, reference = getDistributions(model, ranges)

    ## Get overlap between the initial and final gaussian distributions that can be used to define scores
    overlap = getOverlap(model, ranges, initGaussian, finGaussian, reference)

    ## Use overlap to obtain score that will be used to rank reactions
    scores = overlap.copy()
    
    # scores['score'] = scores['score'].map(lambda x: 1./x)    # Putting this outside of getOverlap adds flexibility
    # when step sizes are vary large (i.e, N=100), there can be no flux overlap. 
    scores['score'] = scores['score'].map(lambda x: 1./x if x != 0 else np.inf)


    # Return scores, and initial and final gaussian distributions
    return scores, initGaussian, finGaussian
    

def getDistributions(model, ranges):
    '''
    Produces gaussian distributions that reflect the flux ranges for each flux in the initial and final 
    flux fractions forced into the production reaction

    Args:
        model (cobra.core.model.Model):
            Genome-scale model to perform simulations
        ranges (pd.DataFrame dictionary):
            FVA values from simulation
    Returns:
        initGaussian (pd.DataFrame):
            A dataframe with means and variances for the initial gaussians for each flux
        finGaussian (pd.DataFrame): 
            A dataframe with means and variances for the final gaussians for each flux
        reference (pd.DataFrame): 
            A dataframe containing the average over fractions of the FVA values captures in ranges
    '''


    ### Define relevant scale (reference) by finding average flux through time
    Ntimes   = len(ranges)
    reference = ranges[0].copy()
    for i,rangeDF in ranges.items():
        reference += rangeDF
    reference = reference.iloc[:,1]/Ntimes

    # Calculate criteria for sorting
    # Find maximum and minimum for each flux for the first two time points
    initial      = pd.DataFrame(np.zeros((len(reference.index),2)), index =reference.index, columns = ['minimum', 'maximum']) 
    initGaussian = pd.DataFrame(np.zeros((len(reference.index),2)), index =reference.index, columns =['loc', 'scale'])    
    
    # initialize the normal gaussians distributed around the FVA reaction values
    for flux, span in initial.iterrows():
        initial.loc[flux,'minimum']    = min(ranges[0].loc[flux,'minimum'],ranges[1].loc[flux,'minimum'])
        initial.loc[flux,'maximum']    = max(ranges[0].loc[flux,'maximum'],ranges[1].loc[flux,'maximum'])

        # gaussian µ is located betwen the FVA min and max of the reaction
        initGaussian.loc[flux,'loc']   =    (initial.loc[flux,'minimum']+initial.loc[flux,'maximum'])/2.0 

        # guassian variance is the abs span between the FVA min and max of the reaction
        initGaussian.loc[flux,'scale'] = abs(initial.loc[flux,'maximum']-initial.loc[flux,'minimum'])

    # Find maximum and minimum for each flux for the final two time points
    final       = pd.DataFrame(np.zeros((len(reference.index),2)), index =reference.index, columns =['minimum', 'maximum']) 
    finGaussian = pd.DataFrame(np.zeros((len(reference.index),2)), index =reference.index, columns =['loc', 'scale'])
    
    # initialize the max product cases of normal gaussians distributed around the FVA reaction values
    for flux, span in final.iterrows():
        final.loc[flux,'minimum']     = min(ranges[Ntimes-2].loc[flux,'minimum'],ranges[Ntimes-1].loc[flux,'minimum'])
        final.loc[flux,'maximum']     = max(ranges[Ntimes-2].loc[flux,'maximum'],ranges[Ntimes-1].loc[flux,'maximum'])

        # gaussian µ is located betwen the FVA min and max of the reaction
        finGaussian.loc[flux,'loc']   =    (final.loc[flux,'minimum']+final.loc[flux,'maximum'])/2.0 

        # guassian variance is the abs span between the FVA min and max of the reaction
        finGaussian.loc[flux,'scale'] = abs(final.loc[flux,'maximum']-final.loc[flux,'minimum'])    # original score

    return initGaussian, finGaussian, reference


def getOverlap(model, ranges, initGaussian,finGaussian, reference):
    '''
    Produces the overlap among initial and final gaussian distributions

    Args:
        model (cobra.core.model.Model):
            Genome-scale model to perform simulations
        ranges (pd.DataFrame dictionary):
            FVA values from simulation
        initGaussian (pd.DataFrame): 
            A dataframe with means and variances for the initial gaussians for each flux
        finGaussian (pd.DataFrame): 
            A dataframe with means and variances for the final gaussians for each flux
        reference (pd.DataFrame): 
            A dataframe containing the average over fractions of the FVA values captures in ranges    

    Returns:
        scores (pd.DataFrame): 
            A dataframe that keeps information on overlaps (scores) for each flux as well as 
            subsystem, trends... etc
    '''

    Ntimes   = len(ranges)
    ## Fit gaussian halfway through the minimum and maximum for the initial and final points
	# Create DataFrame with specified data types for each column
    scores = pd.DataFrame({
		'score': np.zeros(len(reference.index), dtype=np.float64),	# float64
		'trend': '',  # String
		'ors?': '',	# Boolean
		'ands?': '', # Boolean
		'subsystem': '',  # String
		'fluxDiff': np.zeros(len(reference.index), dtype=np.float64)  # float64
	}, index=reference.index)
	
    for flux, span in initGaussian.iterrows():

        # reaction distribution with near zero product flux
        iniDist = norm(initGaussian.loc[flux,'loc'],initGaussian.loc[flux,'scale'])

        # reaction distribution with near theoretical max product flux
        finDist = norm(finGaussian.loc[flux,'loc'], finGaussian.loc[flux,'scale'])        

        #
        x = np.linspace(min(iniDist.ppf(0.001),finDist.ppf(0.001)),max(iniDist.ppf(0.999),finDist.ppf(0.999)), 1000)

        # Find overlap among the initial and final gaussians and use the 1/overlap to sort fluxes
        overlap = sum(iniDist.pdf(x)*finDist.pdf(x))/sum(iniDist.pdf(x))/sum(finDist.pdf(x))
        scores.loc[flux,'score'] = overlap
        #scores.loc[flux,'score'] = 1/overlap  # This will be one outside since it is a particular choice of doing scores

        # Find if flux decrease or increases as fraction increases and store in trend
        trend = 'none'
        if np.abs(initGaussian.loc[flux,'loc']) < np.abs(finGaussian.loc[flux,'loc']):
            trend = 'Up'
        else:
            trend = 'Down'
        scores.loc[flux,'trend'] = trend

        # Find if the gene rule involve any 'or's
        ors = False
        if 'or' in model.reactions.get_by_id(flux).gene_reaction_rule:
            ors = True
        scores.loc[flux,'ors?'] = ors

        # determine absolute flux difference between start and end
        #fluxDiff = ranges[Ntimes-1].loc[flux,'Optimal'] - ranges[0].loc[flux,'Optimal']  # Old way
        initial  = (ranges[0].loc[flux,'maximum'] + ranges[0].loc[flux,'minimum'])/2.0   # Average flux bounds rather than using maximum biomass one
        final    = (ranges[Ntimes-1].loc[flux,'maximum'] + ranges[Ntimes-1].loc[flux,'minimum'])/2.0
        fluxDiff = np.abs(final - initial)
        scores.loc[flux,'fluxDiff'] = fluxDiff

        # Add final flux to table (useful to decide full KOs or not)
        scores.loc[flux,'FinalFlux'] = (ranges[Ntimes-1].loc[flux,'minimum']+ranges[Ntimes-1].loc[flux,'maximum'])/2.0

        # Find if the gene rule involve any 'and's
        ands = False
        if 'and' in model.reactions.get_by_id(flux).gene_reaction_rule:
            ands = True
        scores.loc[flux,'ands?'] = ands
        
        
        # Add gene subsystem for filtering later on.
        rxn = model.reactions.get_by_id(flux)
        subsystem = rxn.subsystem

        if subsystem:
            scores.loc[flux,'subsystem'] = subsystem
        else:
            scores.loc[flux,'subsystem'] = 'N.A.'                                                                               

    ## Return sorted version of ss
    scores = scores.sort_values(by=['score'],ascending=False)

    return scores
    

def filterReactions(UpOrDown,referenceCutOff,fluxScores,fluxRangeDiff,desiredSystems,model,minGrowth,biomassRxn, Ors=False):
    '''
    Filters reactions according to user criteria

    Args:
        UpOrDown (str): 
            select up ('Up') or down ('Down') regulation targets, or both ('Both')
        referenceCutOff (list or int):
            reactions in the production pathway that can be used as a reference 
                            during the scoring evaluation of target reaction significance
                            A number can also be specified as the cutoff 
        fluxScores (pd.DataFrame): 
            A dataframe with fluxes and their corresponding scores
        fluxRangeDiff (float):
            specifies the minimal range between the flux optimal to be considered relevant.   
        desiredSystems (list): 
            subsystems of interest
        Ors (bool): 
            Select reactions that have GPR rules with 'Ors' in them

    Returns:
        selected (pd.DataFrame): 
            A dataframe with fluxes and their corresponding scores for the selected reactions 
            (i.e. a subset of fluxScores)
    '''

    #### First filter reaction through desired subsystems
    fluxScores = filterSubSystem(fluxScores,desiredSystems)

    ##### Calculate desired reference scores
    
    # Cases where a cutoff number was provided
    if isinstance(referenceCutOff,int):
        referenceScore = referenceCutOff
        print('##################')
        print('Your cuttoff score is {:.0f}'.format(referenceCutOff))
        print('The median reaction score in the dataset is {:.0f} and the average is {:.0f}'.format(fluxScores['score'].astype(float).median(),fluxScores['score'].astype(float).mean()))
        print('##################\n')
    
    # Cases where a list of reactions were provided as references to create a cutoff score
    elif referenceCutOff:
        averageList = []
        print('##################')
        print('Your reference reactions scores are')
        
        # Get the score for the reference reactions
        for ref in referenceCutOff:
            f=fluxScores.loc[ref,'score']
            print(f'Reference reaction {ref} score {f}')
            averageList.append(fluxScores.loc[ref,'score'])
        
        # Calculate the average score
        referenceScore = sum(averageList)/len(averageList)
        print('For an cuttoff score of {:.0f}'.format(referenceScore))
        print('The median reaction score in the dataset is {:.0f} and the average is {:.0f}'.format(fluxScores['score'].astype(float).median(),fluxScores['score'].astype(float).mean()))
        print('##################') 
    
    # Cases where reference score was not provided
    else:
        print('No cuttoff score chosen - setting to default of 500')
        print('Use referenceCutOff to specify reactions to use or an interger for a cutoff score if so desired')
        
        # set the reference score to a default of 500 in this instance
        referenceScore = 500

    
    #### Apply filtering criteria

    # Exclude reactions that can be catalyzed by more than one isozymes if user specified
    if Ors==False:
        
        # select reactions based without regard to whether they are correlated or anti-correlated
        if UpOrDown!='Both':
            if print==True:
                print('------------------------------------------\n')
                print('you have selected targets that are {} correlated with biomass flux'.format(UpOrDown))
            selected = fluxScores.loc[(fluxScores['In']==True) & 
                                         (fluxScores['score']>referenceScore) & 
                                         (fluxScores['trend']==UpOrDown) & 
                                         (fluxScores['ors?']==Ors) &
                                          (abs(fluxScores['fluxDiff'])>fluxRangeDiff)]
    
        # select reactions that are either correlated or anti-correlated                   
        else:
            if print==True:
                print('------------------------------------------\n')
                print('you have selected targets that are both up and down- correlated with biomass flux')
            selected = fluxScores.loc[(fluxScores['In']==True) & 
                                         (fluxScores['score']>referenceScore) & 
                                         (fluxScores['ors?']==Ors) &
                                        (abs(fluxScores['fluxDiff'])>fluxRangeDiff)]
    
    # Include reactions that can be catalyzed by more than one isozyme if user specified   
    else:
  
        # select reactions based without regard to whether they are correlated or anti-correlated
        if UpOrDown!='Both':        
            if print==True:
                print('------------------------------------------\n')
                print('you have selected targets that are {} correlated with biomass flux'.format(UpOrDown))
            selected = fluxScores.loc[(fluxScores['In']==True) & 
                                         (fluxScores['score']>referenceScore) & 
                                         (fluxScores['trend']==UpOrDown) &
                                          (abs(fluxScores['fluxDiff'])>fluxRangeDiff)]
        
        # select reactions that are either correlated or anti-correlated                   
        else:
            if print==True:
                print('------------------------------------------\n')
                print('you have selected targets that are both up and down- correlated with biomass flux')
            selected = fluxScores.loc[(fluxScores['In']==True) & 
                                         (fluxScores['score']>referenceScore) &
                                         (abs(fluxScores['fluxDiff'])>fluxRangeDiff)]

    #### Eliminate cases that do not surpass a given growth rate when knocked out
    # Find growth rates when reactions in table are knocked out
    reactions = [x for x in selected.index]
    deletion_results = cobra.flux_analysis.single_reaction_deletion(model, reactions)

    # Add new column to final table
    selected['growth %'] = np.zeros(selected.shape[0])

    # Find maximum growth rate to normalize
    model.objective = model.reactions.get_by_id(biomassRxn)
    solution = model.optimize()
    maxGrowth = solution[biomassRxn]

    # Add growth normalized to maximum growth to table
    for ind in deletion_results.index:
        rxn = next(iter(deletion_results.loc[ind].ids))  # get reaction name from deletion table
        selected.loc[rxn,'growth %'] = deletion_results.loc[ind].growth/maxGrowth  # add normalized growth rate to table

    # Keep only those that pass minimum grow rate and trend down (only cases susceptible to knock out) 
    selected2  = selected.loc[((selected['growth %'] > minGrowth) & (selected['trend']=='Down')) | (selected['trend'] == 'Up')] 


    return selected2


def filterSubSystem(fluxScores,desiredSystems):
    '''
    Filters the flux scores and returns the reactions that are in the desired subsystem
    
    Args:
        fluxScores (pd.DataFrame): 
            Dataframe containing all the scores and information for each reaction
        desiredSystems (list): 
            subsystems interested in exploring
    
    Returns:
        fluxScoresNew (pd.DataFrame):
            A new dataframe with only fluxes in the desiredSystems
    '''
    
    # generate a copy of the scors to work on 
    fluxScoresNew = fluxScores.copy()
    fluxScoresNew['In'] = False 
    
    # check if the reaction is in the desired subsystems
    for flux, span in fluxScoresNew.iterrows():
        if fluxScoresNew.loc[flux,'subsystem'] in desiredSystems:
            fluxScoresNew.loc[flux,'In'] = True
    
    return fluxScoresNew     

def getGenes(model, selected):
    '''
    Gets the genes corresponding to the selected reactions
    
    Args:
        model (cobra.core.model.Model):
            Genome-scale model to perform simulations
        selected (pd.DataFrame): 
            A dataframe with fluxes and their corresponding scores for the selected reactions 
    
    Returns:
        selected (pd.DataFrame):
            An expanded dataframe with new columsn for genes and gene names
    '''

    #### Obtain genes associated with selected reactions
    rxnList = []
    selected['bgene'] = ""
    selected['geneName'] = ""
    
    # iterate through the returned reactions
    for flux, span in selected.iterrows():
        geneName = []

        # obtain reaction GPR
        rxn = model.reactions.get_by_id(flux)
        tempGPR = selected.loc[flux,'bgene']= rxn.gene_reaction_rule
        
        # Divive gene reaction rule into the constituent genes
        systematicNames = [x.strip(' ') for x in tempGPR.replace('(','').replace(')','').replace('and','or').split('or')]
        geneName = [model.genes.get_by_id(g).name for g in systematicNames if g!='']

        # # check if the GPR has multiple associated genes
        # # return multiple genes if so
        # print(tempGPR)
        # if 'and' in tempGPR:
        #     temp = tempGPR.split('and')
        # elif 'or' in tempGPR:
        #     temp = tempGPR.split('or')
        # else:
        #     temp = [tempGPR]

        # print(temp)
        # if len(temp)>1:
        #     for g in temp:
        #         g = g.replace('(','').replace(')','').replace(' ','')
        #         geneName.append(model.genes.get_by_id(g).name)

        # # only 1 gene        
        # elif (temp[0]):
        #     tempGPR = temp[0]
        #     tempGPR = tempGPR.replace('(','').replace(')','').replace(' ','')
        #     geneName.append(model.genes.get_by_id(tempGPR).name)
        
        # # no associated gene
        # else:
        #     geneName.append('')
        selected.at[flux,'geneName']=geneName

    # sort scores descending 
    selected.sort_values(by='score',ascending=False,inplace=True)

    return selected
             