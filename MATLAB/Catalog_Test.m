%% ST - rL-form combination
% Input parameters
N = 32; % Nbr of runs
verbose = true;

% Initiate the catalog
bf = 2.^(0:log2(N)-1); % basic factors
n = log2(N);
k = n-log2(N);
cols = bf;
nMat = [n];
kMat = [k];
colsCell = {cols};

% Iterate through values of k
startcols = cols;
for n = log2(N)+1:15
    if verbose
        fprintf('%i factors\n',n)
    end
    k = n-log2(N);
    candicols = [];
    % Candidate columns for first case
    if n == log2(N)+1
        for ii = 1:N-1
            if ismember(ii,startcols)
                continue
            end
            candicols = [candicols; [startcols ii]];
        end
        % Candidate columns for second case
    else
        for ii = 1:size(startcols,1)
            tempcandicols = STselect(N,startcols(ii,:),startcols(ii,end),"rL",false);
            candicols = [candicols; tempcandicols];
        end
    end
    if verbose
        fprintf('\t%i candidates\n',size(candicols,1))
    end
    outcols = [];
    for ii = 1:size(candicols,1)
        if rLform(candicols(ii,:))
            outcols = [outcols; candicols(ii,:)];
        end
    end
    for ii = 1:size(outcols,1)
        nMat = [nMat; n];
        kMat = [kMat; k];
        colsCell{end+1} = outcols(ii,:);
    end
    startcols = outcols;
    if verbose
        fprintf('\t%i unique designs\n',size(startcols,1))
    end
end
T = table(nMat(:),kMat(:),colsCell','VariableNames',{'n','k','cols'});
