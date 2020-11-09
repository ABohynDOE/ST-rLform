function candicols = STselect(N,cols,lastfac,order,sorted)
    arguments
        N double {mustBeNumeric, mustBePositive, mustBeInteger}
        cols (1,:) {mustBeVector, mustBeInteger, mustBeInRange(cols,1,N,"exclude-upper")}
        lastfac double {mustBeNumeric, mustBePositive, mustBeInteger, mustBeInRange(lastfac,1,N,"exclude-upper")}
        order {mustBeText, mustBeMember(order,{'rL','cL'})} = 'rL'
        sorted {mustBeNumericOrLogical} = true
    end
    % STselect Generate candidate designs using the search-table selection 
    % method.
    %   ST(N,cols,l) generates a matrix where each row contains the
    %   columns for a candidate design, based on the columns in cols. 'l'
    %   indicates the last column added in the parent designs.
    %   ST(N,cols,l,order) generates the matrix of candidate designs, using
    %   the order specified, which can be either conditional lexicographic
    %   'cL' or reversed lexicogrpahic 'rL'. Default is 'rL'.
    %   ST(N,cols,l,order,sorted) sorts the columns of the candidate
    %   designs if sorted is true, and puts the last added columns in last 
    %   position otherwise. Default is true.
    candicols = [];
    idx = 1;
    % Initiate the different orders
    if strcmp(order,'cL')
        [~,ref] = sort(sum(num2col(1:N-1)));
    elseif strcmp(order,'rL')
        ref = 1:N-1;
    end
    for ii = 1:N-1 % All possible added columns
        if ismember(ii,cols) % Reject if already used in the design
            continue
        elseif find(ref==ii) <= find(ref==lastfac) % Reject if upper in the search-table
            continue
        end
        if sorted
            outcols = sort([cols ii]); % Sort columns if needed
        else
            outcols = [cols ii];
        end
        candicols(idx,:) = outcols;
        idx = idx + 1;
    end
end