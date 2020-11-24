function out = rLsmaller(A,B)
    % rLsmaller Determine if a matrix is rL-smaller than the other.
    %   rLsmaller(A,B) returns true if A is rL-smaller than B.
    %   
    %   Let c be the first column in which A and B differ.
    %   Let r be the last element in which A(:,c) and B(:,c) differ
    %   A is rLsmaller than B, is A(r,c) < B(r,c).
    
    boolMat = A~=B;
    if size(A,2) == 1
        colInd = 1;
    else
        colInd = find(any(boolMat,1),1,'first');
    end
    rowInd = find(boolMat(:,colInd),1,'last');
    out = A(rowInd,colInd) < B(rowInd,colInd);
end