function out = rLsmaller(L1,L2)
    % rLsmaller Test if a matrix is rL-smaller than another.
    %   rLsmaller(A,B) returns True if A is rL-smaller than B, and False
    %   otherwise.
    
    a = bi2de(L1');
    b = bi2de(L2');
    mL1 = a(a~=b); % Find first non-identical column in A
    if isempty(mL1)
        out = false; % If A and B are similar, returns False
        return
    end
    mL2 = b(a~=b); % Find first non-identical column in B
    out = mL1(1) < mL2(1); % Check which if is column in A is rL-smaller than column in B
end