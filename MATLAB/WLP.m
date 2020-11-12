function wlp = WLP(cols,s)
    % WLP Word-length pattern of a design.
    %   WLP(cols) outputs the word-length pattern, a vector where the ith
    %   element is the number of words of length i, given the columns of
    %   the design in cols.
    %
    %   WLP(cols,s) outputs the word-length pattern, ommiting the words of
    %   length less than s.
    
    arguments
         cols (:,:) {mustBeNumeric, mustBePositive}
        s double {mustBeInteger, mustBePositive, mustBeNumeric} = 3
    end
    % Basic values
    [nd n] = size(cols);
    r = sum(log2(cols(1,:)) == floor(log2(cols(1,:))));
    k = n-r;
    
    % Basic matrix
    G = Gmat(r);
    W = Gmat(k);
    I = eye(k);
    
    % Initiate WLP matrix
    addcols = cols(:,r+1:end);
    wlp = zeros(nd,n-2);
    len = zeros(nd,2^k-1);
    
    % Compute each word
    for ii = 1:nd
        subG = G(:,addcols(ii,:));
        wordsubG = cat(1,subG,I);
        words = mod(wordsubG*W,2);
        len(ii,:) = sum(words,1);
    end
    len(:,2^k) = n+1;
    
    % Count frequency for WLP
    for ii = 1:nd
        tab = tabulate(len(ii,:));
        wlp(ii,:) = tab(s:end-1,2)';
    end
end