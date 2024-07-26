import cobra
from cobra.flux_analysis import flux_variability_analysis
from cobra import Metabolite, Reaction
import numpy as np
from scipy.stats import norm
import pandas as pd
import time
import matplotlib.pyplot as plt
import core
print('You are using cobra version: ',cobra.__version__)

#np.warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None  # default='warn'

# HGM: Full refactor of code by Hector
class FluxRETAP: 
    ''' 
    Attributes:
        getRecommendations   - Main FluxRETAP routine to obtain gene recommendations
        plotRanges           - Plot the FVA min and max results
        plotDist             - Plot the calculated distributions for each reaction
    '''

    def getRecommendations(model,UpOrDown,productRxn,carbonSourceRxn,biomassRxn,desiredSystems=None,optimalFraction=0.1,
        N=10, Ors=False, referenceCutOff=None, returnGaussian=False,fluxRangeDiff=0.01,minGrowth=0.0,fast=True):

        '''Main routine that produces gene targets for interference/KO or augmentation

        Args:
            model (cobra.core.model.Model):
                Genome-scale metabolic model to perform simulations
            UpOrDown (str):
                select up ('Up') or down ('Down') regulation targets, or both ('Both')
            productRxn (str):
                name of target product to maximize
            carbonSourceRxn (str):
                name of the carbon uptake reaction in the GSM
            biomassRxn (str):
                name of the biomass reaction in the GSM
            desiredSystems (list):
                subsystems of interest, undesired subsystems can be marked with a '~' in the beginning of the name
            optimalFraction (float):
                the fraction of biomass to use
            N (int)
                Number of flux increases through production to simulate
            Ors (bool): 
                Select reactions that have GPR rules with 'Ors' in them
            referenceCutOff (list or int):
                reactions in the production pathway that can be used as a reference 
                                during the scoring evaluation of target reaction significance
                                A number can also be specified as the cutoff
            returnGaussian (bool):
                specifies how much information to return. Default is to return three list and a float.
            fluxRangeDiff (float):
                specifies the minimal range between the flux optimal to be considered relevant.
            fast (bool):
                specifies if we want the fast calculation (calculates only the first and last two fractions) 
                or the complete one (calculates all fractions)

        Returns:
            selectedGenes (pd.DataFrame): dataframe with recommended genes, as well as corresponding scores and other data
            ranges (pd.DataFrame dictonary): flux ranges for each flux and production increase step
            iniDist (pd.DataFrame): Mean and variance for initial gaussian (initial in production step increases)
            finDist (pd.DataFrame): Mean and variance for final gaussian (final in production step increases) 
        '''

        #### Make a list of desired and undesired subsystems. Reactions in undesired subsystems will be eliminated
        allSubsystems = set([r.subsystem for r in model.reactions])  # Create set of all subsystems
        if desiredSystems==None:  # if no desired systems are given, all are returned
            desiredSystems = allSubsystems
        else:
            undesiredSystems = filter(lambda x: len(x) > 0 and x[0]=='~', desiredSystems)  # Undesired subsystems start with '~'
            undesiredSystems = set([x.strip('~') for x in undesiredSystems]) # Strip the initial "~" and put in set

            if undesiredSystems:   # If there are any undesired system, just focus on that, if not just keep what we have
                desiredSystems = allSubsystems.difference(undesiredSystems)


        #### Get ranges of fluxes for each production increase step
        ranges = core.getRanges(model, N, biomassRxn, productRxn, carbonSourceRxn, optimalFraction, fast)

        #### Calculate scores for each flux based on the overlap of flux ranges for initial and final production increase steps
        print("\nCalculating significance scores")

        fluxScores,iniDist,finDist = core.getFluxScoresTable(model, ranges, desiredSystems)

        #### Filter reactions according to user criteria
        selectedRxns = core.filterReactions(UpOrDown,referenceCutOff,fluxScores,fluxRangeDiff,desiredSystems,model,minGrowth,biomassRxn,Ors)

        #### Obtain genes associated with selected reactions
        selectedGenes = core.getGenes(model, selectedRxns)

        #### Return calculations
        # user specifies to return distributions for visualization
        if returnGaussian==False:
        	returnGaussian=None
        if returnGaussian:
            return selectedGenes, ranges, iniDist, finDist
        # else return reactions
        return selectedGenes


    def plotRanges(rxn, fractions, ranges):
        '''
        Displays the FVA min and max flux values for a specified reaction as the product flux is increased.
        
        Args:
            rxn (string): 
                reaction name as present in the GSM
            fractions (list):
                fraction of max product flux at which to display pionts
            ranges (pd.DataFrame):
                FVA values from simulation
        '''

        
        x = np.arange(len(fractions))
        
        # get the FVA optimal data point
        y           = np.array([ranges[i].loc[rxn,'Optimal'] for i in x])
        
        # error bars will be the FVA simulated min and max
        upper_error = abs(np.array([ranges[i].loc[rxn,'maximum'] for i in x]) - y)
        lower_error = abs(y - np.array([ranges[i].loc[rxn,'minimum'] for i in x]))
        asymmetric_error = [lower_error, upper_error]

        # plot
        fig = plt.figure()
        plt.errorbar(fractions, y, yerr=asymmetric_error, fmt='o')
        plt.title(f'Reaction: {rxn}')
        plt.ylabel('Flux')
        plt.xlabel('Fraction to product')
        plt.show()        

    def plotDists(rxn, initGaussian, finGaussian):    
        '''
        Displays the distributions of a specified reactions as production is increased
        
        Args:
            rxn (list): 
                reactions name as present in the GSM to plot               
            initGaussian (pd.DataFrame):
                initial reaction average and spread values to create gaussian distributions
            finGaussian pd.DataFrame):
                final reaction average and spread values to create gaussian distributions

        '''
        for r in rxn:
            
            for flux, span in initGaussian.iterrows(): 
                iniDist = norm(initGaussian.loc[r,'loc'],initGaussian.loc[r,'scale'])
                finDist = norm(finGaussian.loc[r,'loc'], finGaussian.loc[r,'scale'])        
                x = np.linspace(min(iniDist.ppf(0.001),finDist.ppf(0.001)),max(iniDist.ppf(0.999),finDist.ppf(0.999)), 1000)

            ini = iniDist.pdf(x)
            fin = finDist.pdf(x)
            
            # Find overlap among the initial and final gaussians and use the 1/overlap to sort fluxes
            overlap = sum(ini*fin)/sum(ini)/sum(fin)
            score = 1/overlap
            
            fig = plt.figure()

            plt.plot(x, ini, 'r-')   # Red for initial gaussian
            plt.plot(x, fin, 'b-')   # Blue for final gaussian 

            plt.title(f'Reaction: {r} Score: {score}')
            plt.ylabel('Prob')
            plt.xlabel('x')
            plt.show()        



