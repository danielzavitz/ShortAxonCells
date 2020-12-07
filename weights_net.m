function [C] = weights_net(N,hp,ho,m,l,q,lambda)
%Creates a SAC network adjacency matrix
%INPUTS
%       N: number of glom in the network
%       ho: number of processes per oligo SA cell
%       hp: number of processes per poly SA cell
%       m: number of target glomer per glom. 
%       q: probability that a SA cell is oligo
%       l: number of short axon cells per glom
%       lambda: rate parameter of the exponential distribution of SA cell
%               strengths. The mean of strength of connections is 1/lambda
%OUTPUTS
%       C: an NxN adjacency matrix of a SA cell network. 
%% Preallocate adjacency matrix
C =zeros(N);

%% Assign connections
%For each glomerulus
for i=1:N
    %Find the target set of the glomerulus
    v = 1:N;
    v(v==i) = [];
    targ_list =  v(randperm(N-1,m));
     
     for j=1:l
         sa_identity = rand();
         %If SA cells are oligo, assign connections
         if sa_identity < q
             cell_targs = randperm(length(targ_list),ho);
             C((targ_list(cell_targs)),i) = C((targ_list(cell_targs)),i) + exprnd(1/lambda,ho,1);
         end
         %If SA cells are poly, assign connections
         if sa_identity > q
             cell_targs = randperm(length(targ_list),hp);
             C((targ_list(cell_targs)),i) = C((targ_list(cell_targs)),i) + exprnd(1/lambda,hp,1);
         end
             

     end
         
end

end

