function candicols = DOPselect(N,cols)
    arguments
        N double {mustBeNumeric, mustBeInteger, mustBePositive}
        cols (1,:) {mustBeVector, mustBeNumeric, mustBeInRange(cols,1,N,"exclude-upper")}
    end
    % DOPselect Generate candidate designs using the DOP-selection method.
    %   DOPselect(N,cols) generates a matrix where each row contains the
    %   columns for a candidate design, based on the columns in cols.
    outcols = zeros(N-1-log2(N),length(cols)+1);
    index = 1;
    for ii = 1:N-1 % All possible added factors
        if ismember(ii,cols) % Reject if already used in the design
            continue
        end
        wlp = zeros(length(cols)+1,N-1); % Initiate WLP matrix
        for jj = 1:length(cols)+1
            t = [cols ii];
            t(jj) = [];
            w = WLP(t,2);
            wlp(jj,:) = [w zeros(1,(N-1-length(w)))]; % Compute WLP for each DOP
        end
        [~,idx] = sortrows(flip(wlp,1)); % Find the DOP with smallest aberration
        if idx(1) == 1 % Keep candidate if parent designs has MA
            outcols(index,:) = [cols ii];
            index = index +1;
        end
    end
    candicols = outcols(sum(outcols,2)>0,:);
end