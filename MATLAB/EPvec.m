function E = EPvec(D)
    arguments
        D (:,:) {mustBeNumeric, mustBeInRange(D,0,1)}
    end
    % EPVEC Generate the entry-point vector of a design.
    %   EPvec(D) generates the entry-point vector (or the column-number
    %   vector) of a design D.
    Dp = (D*2)-1;
    r = log2(size(D,1));
    Bp = (mod(Rmat(r)*Gmat(r),2)*2)-1;
    M = (transpose(Bp)*Dp)/size(D,1);
    E = (1:size(D,1)-1)*M;
end