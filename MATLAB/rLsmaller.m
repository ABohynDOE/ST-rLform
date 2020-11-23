function out = rLsmaller(A,B)
    % rLsmaller Determine if a matrix is rL-smaller than the other.
    %   rLsmaller(A,B) returns true if A is rL-smaller than B.
    %   
    %   Let c be the first column in which A and B differ.
    %   Let r be the last element in which A(:,c) and B(:,c) differ
    %   A is rLsmaller than B, is A(r,c) < B(r,c).
    
    if size(A,2) == 1
        ind = find(A~=B,1,'last'); 
        out = A(ind) < B(ind);
        return
    end
    colInd = find(all(A~=B),1,'first');
    rowInd = find(A(:,colInd)~=B(:,colInd),1,'last');
    out = A(rowInd,colInd) < B(rowInd,colInd);
end
