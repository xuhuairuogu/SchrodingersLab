function createDefaultPref()
    pref.Nx = 256;
    pref.Nt = 300;
    pref.Lmode = 'periodic';
    pref.pwr = 2;
    pref.plotType = 'density'; 
    pref.fourierPlotType = 'density';
    pref.fourierLines = 7;
    pref.L = 20;
    pref.aLMult = 1;
    pref.dec = 5; %#ok<STRNU>
    
    save('pref.mat', 'pref');  
    
end