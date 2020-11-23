function D = Design(N,cols)
    % DESIGN Generates the design matrix of a N-run regular two-level
    % design.
    %   Design(N,c) generates a N-by-len(c) matrix where the columns in c
    %   are picked from the generalized design matrix B.
    r = log2(N);
    G = Gmat(r);
    R = cat(1,zeros(1,r),flip(G',2));
    D = mod(R*G(:,cols),2);
end