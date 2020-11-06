function w = WLP(cols,s)
    % WLP Word-length pattern of a design.
    %   WLP(cols) outputs the word-length pattern, a vector where the ith
    %   element is the number of words of length i, given the columns of
    %   the design in cols.
    %
    %   WLP(cols,s) outputs the word-length pattern, ommiting the words of
    %   length less than s.
    
    arguments
        cols (1,:) {mustBeVector, mustBeNumeric, mustBePositive}
        s double {mustBeInteger, mustBePositive, mustBeNumeric} = 3
    end
    S = de2bi(cols)'; % Reduced design matrix
    L = S(:,sum(S)>1); % Keep the added factors
    Ls = (L*(-2))+1; 
    model = ff2n(size(L,2)); % Enumerate all added factor combination
    modmat = model(sum(model,2)>0,:); % Only keep 2FI and more
    X = x2fx(Ls,modmat); % Compute wordlength
    Xf = [abs((X-1)/(-2));sum(modmat,2)'];
    vec = sum(Xf,1);
    edges = 1:max(vec)+1;
    w = histcounts(vec,edges); % Summarize into WLP
    w = w(s:end);
end