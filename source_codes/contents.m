% HMM+DM Toolbox
% 
% -----------------------------------------------------------  
% Please refer the following technical report:
% Le Chen, David Barber, and Jean-Marc, 
% Version 1.0  10-May-07
% By Le Chen (chenle02@gmail.com)
% -----------------------------------------------------------
% 
% 
% Data Generation: 
%    GenDir.m               -- Generat samples from single Dirichlet distribution.
%    GenMixtureDir.m        -- Generating samples from a mixture of Dirichlet distributions.
%    GenDynamicMixtureDir.m -- Generate one sequence of random samples from HMM+DM model.
% 
% Model Estimation:
%    MomentMatchingInitDM.m -- Initialize the mixture of Dirichlet by Kmeans + Moment Matching.
%    Dirichlet_loglike.m    -- Evaluate the data's log-likelihoods relative to a Dirichlet model.
%    dataLikelihood_DM.m	  -- Calculate the data likelihood for the Dirichlet mixture.
%    EstDirchlet.m	  -- Estimating single Dirichlet distribution by Newton method.
%    EstMixDirichlet.m	  -- Estimate the parameters of the mixture of Dirichlet by EM algorithm.
%    entropy_base_e.m	  -- Calculate the Entropy using log based on e, instead 2.
%    forback.m		  -- Inference procedure for general HMM + Mixture of density model.
%    EstHMMDM.m		  -- Estimating parameters of HMM+DM by EM algorithm.
%    HardClusterGamma.m     -- Calculating the hard clustering results from the smoothed posteriors.
% 
% Some Scripts:
%    Script_Test_SingleDirichlet.m -- Demo and test script for EstDirchlet.m
%    Script_Test_EstHMMDM.m -- Script to test EstHMMDM.m
% 
% 




