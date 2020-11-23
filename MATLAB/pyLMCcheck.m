function iso = pyLMCcheck(N,cols)
    D = Design(N,cols);
    npMat = py.numpy.array(D).astype('int');
    oaMat = py.oapackage.array_link(npMat);
    LMCval = py.oapackage.LMCcheck(oaMat);
    iso = LMCval == py.oapackage.LMC_MORE();
end