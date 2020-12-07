function [ec,sac,fvals] = dF_dynamics(I,w1,w2,b1,b2,v,eps,a,C)
%Computes the response of glomeruli  with feedback 
%INPUTS
%       I = The OSN inputs for one odor 
%       w1: lower bounds for the change in firing rate of ET cells
%       w2: lower bounds for the change in firing rate of SA cells
%       b1: The steepness of response of ET cells
%       b2: The steepness of response of SA cells
%       eps: The strength of connection from SA to ET
%       a: The strength of connection from ET to SA (typically set to 1)
%       C: The SAC network
%OUTPUTS
%       et: The responses of EC 
%       sa: The responses of SACs
%       fvals: residuals of fsolve


[n,~] = size(I);

    function F = activity(x)
        %The  first 1 to n entries in x are ET cells. The n+1 to 2n entries
        %are SA cells
        
        input1 = I-eps.*C*x(n+1:2*n,1);
        input2 = I+a.*x(1:n,1);
        
        F(1:n) = x(1:n) - dF_sigmoid(w1,b1,v,input1);
        F(n+1:2*n) = x(n+1:2*n) - dF_sigmoid(w2,b2,v,input2);

    end


options = optimoptions('fsolve','Display','off');
[cells,fvals] = fsolve(@activity,rand(2*n,1),options);

ec = cells(1:n);
sac = cells(n+1:2*n);

end
