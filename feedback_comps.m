function [ec,sac,fvals] = feedback_comps(I,w1,w2,b1,b2,v,ep1,a,C)
%Loops through the odors in the odor panel to compute the responses of EC
%and SACs. 
%INPUTS
%See the descriptions of the parameters in "data_generator.m"
%OUTPUTS
%       ec: array of EC responses to odor panel
%       sac: array of SAC responses to odor panel
%       fvals: residuals of the numerical algebraic equation solver

[glomnum,odornum] = size(I);
ec = zeros(glomnum,odornum);
sac = zeros(glomnum,odornum);
fvals = zeros(glomnum*2,odornum);

parfor i=1:odornum
    [ec(:,i),sac(:,i),fvals(:,i)] = dF_dynamics(I(:,i),w1,w2,b1,b2,v,ep1,a,C);
end



end
