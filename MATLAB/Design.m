function D = Design(N,cols)
    arguments
        N double {mustBeNumeric, mustBePositive, mustBeInteger}
        cols (1,:) {mustBeVector, mustBeInteger, mustBeInRange(cols,1,N,"exclude-upper")}
    end
    % DESIGN Generates the design matrix of a N-run regular two-level
    % design.
    %   Design(N,c) generates a N-by-len(c) matrix where the columns in c
    %   are picked from the generalized design matrix B.
    r = log2(N);
    B = mod(Rmat(r)*Gmat(r),2);
    D = B(:,sorted(cols));
end