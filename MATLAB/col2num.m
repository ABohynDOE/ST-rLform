function colvec = col2num(L)
    r = size(L,1);
    rvec = 2.^(0:r-1);
    colvec = rvec*L;
end