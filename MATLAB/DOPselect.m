function candicols = DOPselect(N,cols)
    arguments
        N double {mustBeNumeric, mustBeInteger, mustBePositive}
        cols (1,:) {mustBeVector, mustBeNumeric, mustBeInRange(cols,1,N,"exclude-upper")}
    end
    outcols = zeros(N-1-log2(N),length(cols)+1);
    index = 1;
    for ii = 1:N-1
        if ismember(ii,cols)
            continue
        end
        wlp = zeros(length(cols)+1,N-1);
        for jj = 1:length(cols)+1
            t = [cols ii];
            t(jj) = [];
            w = WLP(t,2);
            wlp(jj,:) = [w zeros(1,(N-1-length(w)))];
        end
        [~,idx] = sortrows(flip(wlp,1));
        if idx(1) == 1
            outcols(index,:) = [cols ii];
            index = index +1;
        end
    end
    candicols = outcols(sum(outcols,2)>0,:);
end