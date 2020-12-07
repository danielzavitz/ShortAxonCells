function [output] = dF_sigmoid(w,b,v,input)
%The logistic function that describes the activity of SAC/EC 
[n,~] = size(input);

%condition such that dF_sigmoid(0) =0;
q = -1+((w-1)./w)^v;
w = w.*ones(n,1);
output = w+(ones(n,1)-w)./(ones(n,1)+q.*exp(-b.*input)).^(1/v);

end