function candicols = STselect(N,cols,lastfac,order,sorted)
    arguments
        N double {mustBeNumeric, mustBePositive, mustBeInteger}
        cols (1,:) {mustBeVector, mustBeInteger, mustBeInRange(cols,1,N,"exclude-upper")}
        lastfac double {mustBeNumeric, mustBePositive, mustBeInteger, mustBeInRange(lastfac,1,N,"exclude-upper")}
        order {mustBeText, mustBeMember(order,{'rL','cL'})} = 'rL'
        sorted {mustBeNumericOrLogical} = true
    end
    candicols = [];
    idx = 1;
    if strcmp(order,'cL')
        [~,ref] = sort(sum(de2bi(1:N-1),2));
    elseif strcmp(order,'rL')
        ref = 1:N-1;
    end
    for ii = 1:N-1
        if ismember(ii,cols)
            continue
        elseif find(ref==ii) <= find(ref==lastfac)
            continue
        end
        if sorted
            outcols = sort([cols ii]);
        else
            outcols = [cols ii];
        end
        candicols(idx,:) = outcols;
        idx = idx + 1;
    end
end