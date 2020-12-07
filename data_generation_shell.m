% A Shell function to generate EC response data for a specific odor panel,
% inhibition strength, and number of trials. 
% Saves the responses of EC & SAC, the adjacency matrices of the networks
% used, and the residuals of the numerical algebraic equation solver

%Load the odor panel
load('omp56_high_no_zeros.mat')
eps = 0.004; % Define the strength of inhibtion
trials = 3; % Number of network iterations to test (typically ~100 on a  high performance computing cluster)
data_generator(I,eps,trials)
