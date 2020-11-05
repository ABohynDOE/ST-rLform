function M = STmat(D)
    arguments
        D (:,:) {mustBeNumeric, mustBeInRange(D,0,1)}
    end
    % STmat Generate the search-table matrix of a design.
    %   M = STmat(D) generates the (N-1)-by-(N-1) search-table of the
    %   N-by-n design D.
    Dp = (D*2)-1;
    r = log2(size(D,1));
    Bp = (mod(Rmat(r)*Gmat(r),2)*2)-1;
    M = (transpose(Bp)*Dp)/size(D,1);
end