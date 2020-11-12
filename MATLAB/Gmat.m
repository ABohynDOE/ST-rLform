function G = Gmat(r)
    arguments
        r {mustBeNumeric, mustBePositive, mustBeInteger}
    end
    % Gmat Generate the generalized reduced design matrix
    %   G = Gmat(r) generates the r-by-(2^r)-1 generalized reduced design
    %       matrix for r basic factors.
    G = de2bi(1:(2^r)-1,'right-msb')';
end
