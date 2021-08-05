function [newss,melt_1,PTF,flag_error,ierror] = melting_master(params,data_REE,INVCOV,grid,stencil)
%

ierror=0;

%%% Setup melting parameters
parametersImport    = importParameters('model_parameters.txt');
[parameters]        = createParameters(parametersImport);
parameters.mcmc     = params;

%% Input from inverse solver
Tp      = params(1);             % [C]
ztop    = params(2)*1e3;        % [m]

%%% Source composition
if length(params) <= 17  % if only Al and Na are inverted
    SiO2_ox     = parameters.source.oxides.SiO2;
    Al2O3_ox    = params(3);
    Fe2O3_ox    = parameters.source.oxides.Fe2O3;
    FeO_ox      = parameters.source.oxides.FeO;
    MnO_ox      = parameters.source.oxides.MnO;
    MgO_ox      = parameters.source.oxides.MgO;
    CaO_ox      = parameters.source.oxides.CaO;
    Na2O_ox     = params(4);
    Cr2O3_ox    = parameters.source.oxides.Cr2O3;
    TiO2_ox     = parameters.source.oxides.TiO2;
    
    FeOt_ox     = FeO_ox + 0.9*Fe2O3_ox;
    
    [SiO2,Al2O3,FeOt,MnO,MgO,CaO,Na2O,Cr2O3,TiO2] = major_creator(SiO2_ox,Al2O3_ox,FeOt_ox,MnO_ox,MgO_ox,CaO_ox,Na2O_ox,Cr2O3_ox,TiO2_ox);
    
    %%% Save parameters
    parameters.source.oxides.SiO2 = SiO2;
    parameters.source.oxides.Al2O3 = Al2O3;
    parameters.source.oxides.FeO = FeOt;
    parameters.source.oxides.MnO = MnO;
    parameters.source.oxides.MgO = MgO;
    parameters.source.oxides.CaO = CaO;
    parameters.source.oxides.Na2O = Na2O;
    parameters.source.oxides.Cr2O3 = Cr2O3;
    parameters.source.oxides.TiO2 = TiO2;
    parameters.source.oxides.FeOt = FeOt;
    
else  % if all the majors are inverted
    Al2O3 = params(3);
    MnO = parameters.source.oxides.MnO;
    MgO = params(19);
    CaO = params(20);
    Na2O = params(4);
    Cr2O3 = params(21);
    TiO2  = parameters.source.oxides.TiO2;
    FeOt = params(18);
    SiO2 = 100 - Al2O3 - MnO - MgO - CaO - Na2O - Cr2O3 - TiO2 - FeOt;
    
    %%% Save parameters
    parameters.source.oxides.Al2O3 = params(3);
    parameters.source.oxides.FeO = params(18);
    parameters.source.oxides.MgO = params(19);
    parameters.source.oxides.CaO = params(20);
    parameters.source.oxides.Na2O = params(4);
    parameters.source.oxides.Cr2O3 = params(21);
    parameters.source.oxides.FeOt = params(18);
    parameters.source.oxides.SiO2 = 100 ...
        - parameters.source.oxides.Al2O3...
        - parameters.source.oxides.MnO...
        - parameters.source.oxides.MgO...
        - parameters.source.oxides.CaO...
        - parameters.source.oxides.Na2O...
        - parameters.source.oxides.Cr2O3...
        - parameters.source.oxides.TiO2...
        - parameters.source.oxides.FeOt;
end
    

parameters.source.oxides.grid = grid;
parameters.source.oxides.stencil = stencil;


%%%%%%% to invert for source composition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 C_source = [SiO2 Al2O3 FeOt MnO MgO CaO Na2O Cr2O3 TiO2];
 C_source_REE = [params(5) params(6) params(7) params(8) params(9) params(10) params(11) params(12) params(13) params(14) params(15) params(16) params(17)];
 C_source_REE = C_source_REE./100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run melting model

[C_solid, C_fluid, C_extracted, C_ave, C_solid_major, C_fluid_major, C_extracted_major, C_ave_major, PTF, flag_error] = perplex_lookup_isentrope(Tp,ztop,parameters,C_source_REE);

% NOTE: If either C_ave or C_ave_major are = 1001 -> there is no melting 
%       If either C_ave or C_ave_major are = 2002 -> isentrope fail
%       If either C_ave or C_ave_major are = 3003 -> traces fail

%% Add data used here
% Trace element composition in melt stored in local vector x
ll=size(C_ave); 
x=C_ave(ll(1),1:ll(2)); 

% major element composition in melt stored in local vector xx
l=size(C_ave_major); 
xx=C_ave_major(l(1),1:l(2)); 

% assemble the right vector for misfit calculations. Must be consistent
% with what was defined in mcmc_wrapper as data

if length(data_REE) == 22
    XX = [x xx]; % traces and majors
elseif length(data_REE) == 13
    XX = x;    % only trace elements
elseif length(data_REE) == 9
    XX = xx;   % only major elements
end

%% CALCULATE REE MISFIT

if x == 1001  % no melt produced
    newss = Inf;
    disp('first option pause - melting master')
elseif x == 2002  % 
    newss = Inf;
    ierror=1;
    disp('second option pause - melting master')
elseif x == 3003  % 
    newss = Inf;
    ierror=1;
    disp('third option pause - melting master')
elseif any(flag_error == 1)
    newss = Inf;
    ierror=1;
    disp('fourth option pause - melting master')    
else

[newss] = misfit(data_REE,INVCOV,XX);

end


% convert to another variable...
melt_1=XX;


end
