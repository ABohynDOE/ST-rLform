function rL = rLmin(cols,lastfac)
    arguments
        cols (1,:) {mustBeVector, mustBeNumeric, mustBePositive}
        lastfac double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = cols(end)
    end
    % RLFORM Determine if a reduced design matrix in rL-form is rL-minimal.
    %   rLmin(cols) outputs true if the reduced design matrix, in rL-form, 
    %   formed the columns in cols is rL-minimal.
    %   rLform(cols,lastfac) defines lastfac as the number of the last
    %   column added to the candidate design.
    
    S = num2col(sort(cols)); % Reduced design matrix in rL order
    rL = true;
    rlst = nchoosek(1:size(S,2),size(S,1)); % All possible set of basic factors
    p = perms(1:size(S,1)); % All possible row permutations
    Lref = S(:,sum(S,1) > 1); % Added factors of the candidate design
    for ii = 1:length(rlst) % For each set of r basic factors
%         if ~ismember(lastfac,rlst(ii,:)) % Reject if last factor not included
%             continue
%         end
        R = S(:,rlst(ii,:)); % New basic factor matrix
        if cond(R) >= 1/eps % Test conditional number for inversion
            continue; % Reject if it is singular
        end
        K = S(:,setdiff(1:1:size(S,2),rlst(ii,:))); % Matrix of non-basic factors
        L = mod(R\K,2); % Recompute new added factors
        if any(any(~ismember(L,[0 1]))) % Checks that the L matrix is binary
            continue
        end
        for jj = 1:length(p) % For all row permutations
            Lstar = L(p(jj,:),:); % Permute the rows
            [~,idx] = sort(col2num(Lstar));
            Lstarmin = Lstar(:,idx); % Sort the colums of L in rL-order
            if rLsmaller(Lstarmin,Lref) % Test if new L is rL-smaller
                rL = false;
                return
            end
        end
    end
end
