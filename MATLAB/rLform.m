function rL = rLform(cols,lastfac)
    arguments
        cols (1,:) {mustBeVector, mustBeNumeric, mustBePositive}
        lastfac double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = cols(end)
    end
    % RLFORM Determine if a reduced design matrix is in rL-form.
    %   rLform(cols,lastfac) outputs true is S is in rL-form and the l-th column of
    %   S if included in the basic factor set.
    %   rLform(S) outputs true if S is in rL-form.
    
    S = de2bi(sort(cols))'; % Reduced design matrix in rL order
    rL = true;
    rlst = nchoosek(1:size(S,2),size(S,1)); % B.F. set possibilities
    p = perms(1:size(S,1)); % Row permutations possibilities
    Lref = S(:,sum(S,1) > 1); % Added factor only
    for ii = 1:length(rlst) % For each set of r basic factors
%         if ~ismember(lastfac,rlst(ii,:)) % Reject if last factor not included
%             continue
%         end
        R = S(:,rlst(ii,:));
        if det(R) == 0
            continue; % Reject if the set forms a singular matrix
        end
        K = S(:,setdiff(1:1:size(S,2),rlst(ii,:))); % Not chosen as basic factors
        L = mod(R\K,2); % New matrix of the added factors
        if any(any(~ismember(L,[0 1]))) %Checks that the L matrix is binary (not singular)
            continue
        end
        for jj = 1:length(p) % For all row permutations
            Lstar = L(p(jj,:),:); % Permute the rows
            [~,idx] = sort(bi2de(Lstar'));
            Lstarmin = Lstar(:,idx); % Sort the colums of L in rL-order
            if rLsmaller(Lstarmin,Lref) % Test if new L is rL-smaller
                rL = false;
                return
            end
        end
    end
end
