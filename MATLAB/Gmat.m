function G = Gmat(r)
arguments
    r {mustBeNumeric, mustBePositive, mustBeInteger}
end
% Gmat Generate the generalized reduced design matrix
%   G = Gmat(r) generates the r-by-(2^r) generalized reduced design
%       matrix for r basic factors.

%G = de2bi(1:(2^r)-1,'right-msb')';
G = zeros(r,2^r);
G(1,1) = 1;
for ii = 2:(2^r-1)
    if log2(ii) == floor(log2(ii))
        G(log2(ii)+1,ii) = 1;
    else
        a = 2^(floor(log2(ii)));
        b = ii-a;
        G(:,ii) = mod(G(:,a)+G(:,b),2);
    end
end
end
