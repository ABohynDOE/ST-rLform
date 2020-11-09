function colvec = col2num(L)
    r = size(L,1);
    n = size(L,2);
    rvec = 2.^(0:r-1);
    colvec = zeros(1,n);
    for ii = 1:n
        colvec(ii) = rvec*L(:,ii);
    end
end