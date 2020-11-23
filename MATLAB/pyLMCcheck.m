function iso = pyLMCcheck(Dmat)
    npMat = py.numpy.array(Dmat).astype('int');
    oaMat = py.oapackage.array_link(npMat);
    LMCval = py.oapackage.LMCcheck(oaMat);
    iso = LMCval == py.oapackage.LMC_MORE();
end