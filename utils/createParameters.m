function [parameters]=createParameters(initialParameters)
% function [mesh,particles]=createInitial(initialImport,NTX,MTXS,NCX,MCXS)

parameters.model.fM = initialParameters.value(initialParameters.name == 'fM');
parameters.model.fT = initialParameters.value(initialParameters.name == 'fT');
parameters.model.dP = initialParameters.value(initialParameters.name == 'dP');
parameters.model.g = initialParameters.value(initialParameters.name == 'g');
parameters.model.W0 = initialParameters.value(initialParameters.name == 'W0');
parameters.model.P_ini = initialParameters.value(initialParameters.name == 'P_ini');
parameters.model.k0 = initialParameters.value(initialParameters.name == 'k0');
parameters.model.n = initialParameters.value(initialParameters.name == 'n');
parameters.model.isen = initialParameters.value(initialParameters.name == 'isen');
parameters.model.inter = initialParameters.value(initialParameters.name == 'inter');
parameters.model.method = initialParameters.value(initialParameters.name == 'method');
parameters.model.diff = initialParameters.value(initialParameters.name == 'diff');
parameters.model.part = initialParameters.value(initialParameters.name == 'part');
parameters.model.zou = initialParameters.value(initialParameters.name == 'zou');
parameters.model.parallel = initialParameters.value(initialParameters.name == 'parallel');
parameters.model.rho_m = initialParameters.value(initialParameters.name == 'rho_m');
parameters.model.rho_l = initialParameters.value(initialParameters.name == 'rho_l');
parameters.model.alpha_m = initialParameters.value(initialParameters.name == 'alpha_m');
parameters.model.Cp_m = initialParameters.value(initialParameters.name == 'Cp_m');
parameters.model.mu_l = initialParameters.value(initialParameters.name == 'mu_l');
parameters.model.beta = initialParameters.value(initialParameters.name == 'beta');
parameters.source.oxides.SiO2 = initialParameters.value(initialParameters.name == 'SiO2');
parameters.source.oxides.Al2O3 = initialParameters.value(initialParameters.name == 'Al2O3');
parameters.source.oxides.Fe2O3 = initialParameters.value(initialParameters.name == 'Fe2O3');
parameters.source.oxides.FeO = initialParameters.value(initialParameters.name == 'FeO');
parameters.source.oxides.MnO = initialParameters.value(initialParameters.name == 'MnO');
parameters.source.oxides.MgO = initialParameters.value(initialParameters.name == 'MgO');
parameters.source.oxides.CaO = initialParameters.value(initialParameters.name == 'CaO');
parameters.source.oxides.Na2O = initialParameters.value(initialParameters.name == 'Na2O');
parameters.source.oxides.Cr2O3 = initialParameters.value(initialParameters.name == 'Cr2O3');
parameters.source.oxides.TiO2 = initialParameters.value(initialParameters.name == 'TiO2');
parameters.source.traces.La = initialParameters.value(initialParameters.name == 'La');
parameters.source.traces.Ce = initialParameters.value(initialParameters.name == 'Ce');
parameters.source.traces.Pe = initialParameters.value(initialParameters.name == 'Pe');
parameters.source.traces.Pr = initialParameters.value(initialParameters.name == 'Pr');
parameters.source.traces.Nd = initialParameters.value(initialParameters.name == 'Nd');
parameters.source.traces.Sm = initialParameters.value(initialParameters.name == 'Sm');
parameters.source.traces.Eu = initialParameters.value(initialParameters.name == 'Eu');
parameters.source.traces.Gd = initialParameters.value(initialParameters.name == 'Gd');
parameters.source.traces.Tb = initialParameters.value(initialParameters.name == 'Tb');
parameters.source.traces.Dy = initialParameters.value(initialParameters.name == 'Dy');
parameters.source.traces.Ho = initialParameters.value(initialParameters.name == 'Ho');
parameters.source.traces.Er = initialParameters.value(initialParameters.name == 'Er');
parameters.source.traces.Yb = initialParameters.value(initialParameters.name == 'Yb');
parameters.source.traces.Lu = initialParameters.value(initialParameters.name == 'Lu');



