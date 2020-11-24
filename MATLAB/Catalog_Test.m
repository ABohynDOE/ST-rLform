%% ST - rL-form combination
% Input parameters
N = 32; % Nbr of runs
verbose = true; % Print info about enumeration
maxfac = 15; % Max number of factors in the catalog

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
totalTime = 0;
for n = log2(N)+1:maxfac
    tic;
    if verbose
        fprintf('%i factors:',n)
    end
    k = n-log2(N);
    candicols = [];
    % Candidate columns for first added factor
    if n == log2(N)+1
        for ii = 1:N-1
            if ismember(ii,startcols)
                continue
            end
            candicols = [candicols; [startcols ii]];
        end
    % Candidate columns for the rest of the added factors
    else
        % Candidate columns selection
        for ii = 1:size(startcols,1)
            tempcandicols = STselect(N,startcols(ii,:),startcols(ii,end),"rL",false); %Search-table for candidate columns
            candicols = [candicols; tempcandicols];
        end
    end
    if verbose
        fprintf('\t%i candidates',size(candicols,1))
    end
    outcols = [];
    % Isomorphism reduction of the candidate designs set
    for ii = 1:size(candicols,1)
        if rLmin(candicols(ii,:)) % Only keep rL-minimal designs
            outcols = [outcols; candicols(ii,:)];
        end
    end
    for ii = 1:size(outcols,1) % Add designs to the catalog
        nMat = [nMat; n];
        kMat = [kMat; k];
        colsCell{end+1} = outcols(ii,:);
    end
    startcols = outcols;
    t = toc;
    totalTime = totalTime +t;
    if verbose
        fprintf('\t%i representatives\t(%.2f sec.)\n',size(startcols,1),t)
    end
end
fprintf('Catalog for %i runs generated in %.2f seconds\n',N,totalTime);
T = table(nMat(:),kMat(:),colsCell','VariableNames',{'n','k','cols'});


%% Summary table of the number of designs
sumMat = zeros((maxfac-log2(N)),2);
for ii = log2(N):maxfac
    rows = T.n == ii;
    sumMat(ii-log2(N)+1,1) = ii;
    sumMat(ii-log2(N)+1,2) = numel(T(rows,1));
end
% Dipslay the table
if verbose
    sumMat
end
