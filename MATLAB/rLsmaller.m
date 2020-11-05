function out = rLsmaller(L1,L2)
    a = bi2de(L1');
    b = bi2de(L2');
    mL1 = a(a~=b);
    if isempty(mL1)
        out = false;
        return
    end
    mL2 = b(a~=b);
    out = mL1(1) < mL2(1);
end