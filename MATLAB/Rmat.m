function R = Rmat(r)
    arguments
        r {mustBeNumeric, mustBePositive, mustBeInteger}
    end
    % Rmat Generate the basic factor matrix.
    %   R = Rmat(r) generates the (2^r)-by-r basic factor matrix for r
    %   basic factor.
    R = zeros(2^r,r);
    for ii = 1:r
        a = (2^r)/(2^ii);
        b = (2^r)/(2*a);
        R(:,ii) = repmat(cat(1,repmat(0,a,1),repmat(1,a,1)),b,1);
    end
end