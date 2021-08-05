% -------------------------------------------------------------------------
%  mcmc_wrapper.m
%  
%  This is the main driver for running MCMC inversions of trace and/or 
%  major elements from mantle-derived melts for the thermochemical state
%  of the source (mantle).
%
%  After defining some initial parameters that control the inversion,
%  the script calls the function 'amrun', which uses the Adaptive
%  Metropolis algorithm and runs the inversion.
% 
%  The main melting model parameters are defined in 'model_parameters.txt'
%  within 'setup' folder. Other trace element model parameters are defined
%  in 'trace_elements_properties'.
%
%  There are two ways to compute the isentropic melting path (isen = 1 and
%  2 in model_parameters.txt). 
%
%  isen == 1 ; uses precomputed isentropic paths over a grid.
%  isen == 2 ; computes and saves isentropic paths on the fly (over a grid)
%
%  Default parameters reproduce the results in Section 5, where both major
%  and trace element compositional data has been used to infer mantle
%  potential temperature (Tp), final depth of melting (ztop), and source
%  composition (Al2O3, Na2O and REE's). Follow the comments below if you
%  want to include more oxides in the inversion (i.e. using option isen =
%  2, and expanding 'Sig_sol', 'par0' and 'bounds'). 
%
%  Contant the authors with any question or comment.
% 
%  Copyright (C) 2021  B. Oliveira, J.C. Afonso & M. Klocking
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  Reference: Oliveira, B., Afonso, J.C., Klocing, M. (2021), Melting 
%             conditions and mantle source composition from probabilistic 
%             joint inversion of major and rare earth element concentrations,
%             Geochim. Cosmochim. Acta, doi:
% -------------------------------------------------------------------------

%% file paths & names
clear all;
clc;
matlabrc;

addpath([pwd,filesep,'amcode/']); 
addpath([pwd,filesep,'utils/']);
addpath([pwd,filesep,'setup/']);
addpath([pwd,filesep,'perplex-model/'])
addpath([pwd,filesep,'perplex-model/tables'])
% Add the path were the folder 'Variables_Table_extra_aux' is:
addpath([pwd,filesep,'perplex-model/Variables_Table_extra_aux'])

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool;
else
    poolsize = p.NumWorkers;
end

%% COMPILE MEX FILES
working_folder = pwd;
specific_folder_1 = 'perplex-model/mex_solid';
specific_folder_2= 'perplex-model/mex_melt';
perplex_folder_1=sprintf('%s/%s',working_folder,specific_folder_1);
perplex_folder_2=sprintf('%s/%s',working_folder,specific_folder_2);

%% Compile Only Solid
cd(perplex_folder_1)

mex meemum_fun_solid_gateway.F meemum_fun.f

cd(working_folder)

%% Compile Solid and Melt
cd(perplex_folder_2)

mex meemum_fun_melt_gateway.F meemum_fun.f

cd(working_folder)

%rng shuffle
%% Input parameters for MCMC inversion

nsimu       = 220000;            % number of mcmc steps  
npar        = 4+13;              % number of model parameters 
adascale    = 2.4/sqrt(npar);    % scale factor for adaptation
save_chain  = 5000;              % how often do we save the chain?     

% Set up the initial proposal covariance matrix for sampling during the
% pre-adaptive stage. After adaptation, this covariance matrix will be
% modified. Below we start with a simple diagonal matrix. The ordering of 
% the elements in the matrix is as follows:
% Tp, z_top, Al2O3, Na2O, La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Yb,
% Lu.

Sig_sol = diag([180 3 0.002 0.0006  0.8 1 0.6 1 ...
              0.6 0.25 0.6 0.2 0.8 0.2 0.6 0.6 0.06]);
          
%
% Note that this matrix has to change depending on whether we are inverting
% source composition of majors + traces, traces only, majors only, reduced 
% majors (i.e. majors as a function of Al and Na content) or full majors 
% (i.e. all 7 major elements are part of the parameter vector). The 
% ordering of the elements in the matrix for the majors +traces is as 
% follows:
% Tp, z_top, Al2O3, Na2O, La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Yb,
% Lu, FeOt, MgO, CaO, Cr2O3
% 
% TiO2 and MnO are part of the major element composition and assumed to be
% constant (as defined in model_parameters.txt). SiO2 is computed as:
% SiO2 = 100 - sum(other oxides);

         
%% Chemical data/input
% The ordering of the data vector is
% La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Yb, Lu + majors in the order
% SiO2 Al2O3 FeO MnO MgO CaO Na2O Cr2O3 TiO2
% Assign zero to the values that you do not want to account for

mu = [131.5 97.6 81.8 67.1 41.9 35.2 29.6 24.4 19.9 17.2 15.0 13.3 13.8 ...;  
     0.4824  0.1447  0.0954  0   0.1346  0.1104   0.0325     0       0]; 

% Assign the variance (uncertainties) for data. These should also include
% the modelling uncertainties. In the following we use a simple estimate
% as a percentage of the inputs.

% vector of sigmas for traces
int1=[mu(1)*(13/100) mu(2)*(12/100) mu(3)*(11/100) mu(4)*(9/100) mu(5)*(7.5/100) mu(6)*(7.5/100) mu(7)*(6.3/100) mu(8)*(7/100) mu(9)*(6.26/100) mu(10)*(7./100) mu(11)*(8./100) mu(12)*(7.6/100) mu(13)*(14./100)];
int1 = int1.^2; % variance
% vector of sigmas for majors
%      SiO2           Al2O3          FeO              MnO              MgO             CaO             Na2O             Cr2O3            TiO2
int2=[mu(14)*(3/100) mu(15)*(5.5/100) mu(16)*(8.6/100) mu(17)*(7.5/100) mu(18)*(15./100) mu(19)*(7.5/100) mu(20)*(11./100) mu(21)*(6/100) mu(22)*(9./100)];
int2=int2.^2;   % variance

var=[int1 int2]; % use data from both traces and majors - change as necessary
var(17)=1; var(21)=1; var(22)=1; % necessary for data points with entries = 0

Sig = eye(length(mu)).*var';     %... covariance matrix as diagonal matrix
Lam = inv(Sig);

%% Hard bounds
% The inversion will not explore parameters outside these bounds
Tmin      = 1250;    Tmax      = 1520;      % parameter 1
zmin      = 25;      zmax      = 110;       % parameter 2
Al2O3min  = 3.1;     Al2O3max  = 4.2;       % parameter 3
Na2Omin   = 0.1;     Na2Omax   = 0.6;       % parameter 4
%---------- Traces --------------------------------------------------------
Lamin = 19;      Lamax = 68.3;              % parameter 5
Cemin = 55;      Cemax = 175.2;             % parameter 6
Prmin = 10.7;    Prmax = 26.5;              % parameter 7
Ndmin = 58.1;    Ndmax = 134.1;             % parameter 8
Smmin = 23.9;    Smmax = 43.4;              % parameter 9
Eumin = 9.6;     Eumax = 16.6;              % parameter 10
Gdmin = 35.8;    Gdmax = 58.5;              % parameter 11
Tbmin = 7;       Tbmax = 10.7;              % parameter 12 
Dymin = 50.5;    Dymax = 72.4;              % parameter 13
Homin = 11.5;    Homax = 15.9;              % parameter 14
Ermin = 34.8;    Ermax = 46.8;              % parameter 15 
Ybmin = 36.5;    Ybmax = 47.7;              % parameter 16
Lumin = 5.0;     Lumax = 7.5;               % parameter 17
%---------- Majors --------------------------------------------------------
FeOmin   = 6;    FeOmax   = 9;              % parameter 18
MgOmin   = 36;   MgOmax   = 48;             % parameter 19
CaOmin   = 0.25; CaOmax   = 3.75;           % parameter 20
Cr2O3min = 0.2;  Cr2O3max = 2.2;            % parameter 21
% SiO2 as a function of all the others
SiO2min  = 100 - Al2O3max - Na2Omax - FeOmax - MgOmax - CaOmax - Cr2O3max;
SiO2max  = 100 - Al2O3min - Na2Omin - FeOmin - MgOmin - CaOmin - Cr2O3min;

%% format input parameters for amrun
% create input arguments for the amrun function
clear model params options

% Define the sum-of-squares function (i.e. misfit function)
model.ssfun    = @melting_master;

% arguments for ssfun are in data
data  = struct('mu',mu,'Sig',Sig,'Lam',Lam);

% Underlying discretization - DO NOT TOUCH without consulting the authors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dTp     = 1;            % [C]
dZtop   = 0.1;          % [km]
dAl2O3  = 0.01;         % [%wt]
dNa2O   = 0.01;         % [%wt]
dFeO    = 0.01;         % [%wt]
dMgO    = 0.01;         % [%wt]
dCaO    = 0.01;         % [%wt]
dCr2O3  = 0.01;         % [%wt]
dSiO2   = 0.01;         % [%wt]

grid_Tp     =  (1100:dTp:1700);
grid_Ztop   =  (0:dZtop:200);
grid_Al2O3  =  (0:dAl2O3:100);
grid_Na2O   =  (0:dNa2O:100);
grid_FeO    =  (0:dFeO:100);
grid_MgO    =  (0:dMgO:100);
grid_CaO    =  (0:dCaO:100);
grid_Cr2O3  =  (0:dCr2O3:100);
grid_SiO2   =  (0:dSiO2:100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stencil for grid - Needs to be multiple of 0.01
dTp_stencil     = 5;   % [C]
dZtop_stencil   = 1;   % [km]
dAl2O3_stencil  = 0.2;
dNa2O_stencil   = 0.05;
dFeO_stencil    = 0.01;
dMgO_stencil    = 0.01;
dCaO_stencil    = 0.01;
dCr2O3_stencil  = 0.01;
dSiO2_stencil   = 0.01;

% Check if stencil is valid
if fix(dTp_stencil/dTp) ~= dTp_stencil/dTp || ...
        fix(dZtop_stencil/dZtop) ~= dZtop_stencil/dZtop || ...
        fix(dAl2O3_stencil/dAl2O3) ~= dAl2O3_stencil/dAl2O3 || ...
        fix(dNa2O_stencil/dNa2O) ~= dNa2O_stencil/dNa2O || ...
        fix(dFeO_stencil/dFeO) ~= dFeO_stencil/dFeO || ...
        fix(dMgO_stencil/dMgO) ~= dMgO_stencil/dMgO || ...
        fix(dCaO_stencil/dCaO) ~= dCaO_stencil/dCaO || ...
        fix(dCr2O3_stencil/dCr2O3) ~= dCr2O3_stencil/dCr2O3 || ...
        fix(dSiO2_stencil/dSiO2) ~= dSiO2_stencil/dSiO2
    error('Warning: Check stencil as it is not multiple of the underlying grid')
end

grid_Tp_stencil     =  [1100:dTp_stencil:1700];
grid_Ztop_stencil   =  [0:dZtop_stencil:200];
grid_Al2O3_stencil  =  [0:dAl2O3_stencil:100];
grid_Na2O_stencil   =  [0:dNa2O_stencil:100];
grid_FeO_stencil    =  [0:dFeO_stencil:100];
grid_MgO_stencil    =  [0:dMgO_stencil:100];
grid_CaO_stencil    =  [0:dCaO_stencil:100];
grid_Cr2O3_stencil  =  [0:dCr2O3_stencil:100];
grid_SiO2_stencil   =  [0:dSiO2_stencil:100];


% random starting location (initial point in the mcmc inversion)
% extend with other oxide information if inverting for all majors
% Order of parameters and bounds same as in Sig_sol.
params.par0    = [1390,55,3.5,0.31,53.75,125.1,18.6,100.1,35.65,13.1,47.15, ...
                 8.85,61.45,13.7,40.8,42.1,6.45] ;     % initia l Tp, ztop, Al2O3, Na2O, and REE [C, km, %wt, %wt, ppm]
params.bounds  = [Tmin zmin Al2O3min Na2Omin Lamin Cemin Prmin Ndmin ...
                  Smmin Eumin Gdmin Tbmin Dymin Homin Ermin Ybmin Lumin;
                  Tmax zmax Al2O3max Na2Omax Lamax Cemax Prmax Ndmax Smmax ...
                  Eumax Gdmax Tbmax Dymax Homax Ermax Ybmax Lumax];
              
if length(params.par0) ~= size(params.bounds,2) || ...
        length(params.par0) ~= size(Sig_sol,2) 
    error('Error: Dimensions of par0, bounds and Sig_sol are inconsistent')
end
              
params.grid.Tp       = grid_Tp;
params.grid.Ztop     = grid_Ztop;
params.grid.Al2O3    = grid_Al2O3;
params.grid.Na2O     = grid_Na2O;
params.grid.FeO      = grid_FeO;
params.grid.MgO      = grid_MgO;
params.grid.CaO      = grid_CaO;
params.grid.Cr2O3    = grid_Cr2O3;
params.grid.SiO2     = grid_SiO2;

params.stencil.Tp    = grid_Tp_stencil;
params.stencil.Ztop  = grid_Ztop_stencil;
params.stencil.Al2O3 = grid_Al2O3_stencil;
params.stencil.Na2O  = grid_Na2O_stencil;
params.stencil.FeO   = grid_FeO_stencil;
params.stencil.MgO   = grid_MgO_stencil;
params.stencil.CaO   = grid_CaO_stencil;
params.stencil.Cr2O3 = grid_Cr2O3_stencil;
params.stencil.SiO2  = grid_SiO2_stencil;              

%% Options for the AM algorithm
options.nsimu       = nsimu;
options.adaptint    = 15000;                    % adapt every "adaptint" simulations                     
options.qcov        = Sig_sol.*2.4^2./npar;     %...initial proposal covariance matrix
options.int_fac     = 2;                        % multiplier for first non-adaptive stage (*)
options.burn        = 30000;                    % burn-in 
options.adascale    = adascale;        
options.save_chain  = save_chain;               % how often do we save the chain?
options.printint    = 200;
options.plot        = 0;                        % if 0 no real-time plots; if = 1 -> real time plots (for testing)

% (*) it means that the first pre-adaptive stage will be (int_fac*adaptint)
% long. After that, adaptivity will be erformed every adaptint mcmc steps

%% run inversion 
% control the random number generator
rng shuffle;

[results,chain,misfit,misfit2,melt_erupt,melt_path,s2chain,boundError,error] = amrun(model,data,params,options);

%% plot
% plot results
% plotmatrix(chain);
plots_results_RGR;

%% save output

filename = strcat(pwd,'/outputs/results.mat');
save(filename,'-struct', 'results');
filename = strcat(pwd,'/outputs/options.mat');
save(filename, 'options');
filename = strcat(pwd,'/outputs/params.mat');
save(filename, 'params');
filename = strcat(pwd,'/outputs/data.mat');
save(filename, 'data');
filename = strcat(pwd,'/outputs/chain.mat');
save(filename, 'chain');
filename = strcat(pwd,'/outputs/misfit.mat');
save(filename, 'misfit');
filename = strcat(pwd,'/outputs/misfit2.mat');
save(filename, 'misfit2');
filename = strcat(pwd,'/outputs/melt_erupt.mat');
save(filename, 'melt_erupt');
filename = strcat(pwd,'/outputs/melt_path.mat');
save(filename, 'melt_path');
filename = strcat(pwd,'/outputs/s2chain.mat');
save(filename, 's2chain');
filename = strcat(pwd,'/outputs/boundError.mat');
save(filename, 'boundError');
filename = strcat(pwd,'/outputs/error.mat');
save(filename, 'error');

