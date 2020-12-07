function [] = data_generator(I,eps,trials)
%Computes the responses to an odor panel of across realizations of
%selective and nonselective networks

%Inputs
        %I: odor panel, rows correspond to OSN, columns to odors
        %eps: the overall strength of inhibition in the network
        %trials: the number of network iterations (per selectivity value)
        %        to test
%OUTPUTS
        %Saves the EC,SAC, networks, and  numerical algebra residuals of the network

%% Parameters


% Network parameters
[N,~] = size(I); % Number of nodes in network
hp = 20;    % Number of connections made by polyglom. SACs
ho = 4;     % Number of connections made by oligoglom. SACs
l = 40;     % Number of SACs per glom.
m_vect = [20 N-1]; % Target size (selectivity)
q = 0.8;           % Prob. of a SAC
lambda = 0.8;      % Exponential weight parameter

% SAC (indexed 1), EC (indexed 2),  and common response equation params.
% Equations are a logistic function describing change in firing rate from
% baseline
b1=70; % SAC Growth Rate
b2=10; % EC Growth Rate
w1 = -0.1; % SAC activity lower bound (negative activity = suppression)
w2 = -0.05; % EC activity lower bound
v = 2.5; % Determines where (compared to the upper/lower bounds)
a = 1; %strength of exciation from EC to SAC in each glom. 


%% Preallocate data matrices
ec_data = cell(trials,length(m_vect));
sa_data = cell(trials,length(m_vect));
c_data = cell(trials,length(m_vect));
fvals = cell(trials,length(m_vect));


rng shuffle
%% Computations
for i=1:length(m_vect)
    for j=1:trials
        C = weights_net(N,hp,ho,m_vect(i),l,q,lambda);
        c_data{j,i} = sparse(C);
        [ec_data{j,i},sa_data{j,i},fvals{j,i}] = feedback_comps(I,w1,w2,b1,b2,v,eps,a,C);
        
    end
end



%% Save Data
lambda_val = num2str(lambda);
ep_val = num2str(eps);
a_val = num2str(a);

lambda_val(lambda_val=='.') = [];
ep_val(ep_val == '.') = [];
a_val(a_val == '.') = [];

filename1 = ['weights_dF_feedback_' 'ec_' 'lambda_' lambda_val 'ep_' ep_val 'a_val_' a_val];
filename2 = ['weights_dF_feedback_' 'library' 'lambda_' lambda_val 'ep_' ep_val 'a_val_' a_val];
filename3 = ['weights_dF_feedback_' 'fvals_' 'lambda_' lambda_val 'ep_' ep_val 'a_val_' a_val];
filename4 = ['weights_dF_feedback_' 'sa_' 'lambda_' lambda_val 'ep_' ep_val 'a_val_' a_val];

save(filename1,'ec_data')
save(filename2,'c_data')
save(filename3,'fvals')
save(filename4,'sa_data')



end  
        
        
