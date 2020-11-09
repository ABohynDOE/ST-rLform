function numvec = num2col(v)
    rmax = floor(log2(max(v)))+1;
    G = Gmat(rmax);
    numvec = G(:,v);
end