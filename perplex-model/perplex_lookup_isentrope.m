function [C_solid, C_fluid, C_extracted, C_ave, C_solid_major, C_fluid_major, C_extracted_major, C_ave_major, PTF, flag_error] = perplex_lookup_isentrope(Tp,ztop,parameters,C_source_REE)
% PERPLEX_LOOKUP_ISENTROPE
%
% [C_solid, C_fluid, C_extracted, C_ave, C_solid_major, C_fluid_major, C_extracted_major, C_ave_major, PTF] = perplex_lookup_isentrope(Tp,ztop,oxides,parameters,C_source_REE)
%
% INPUT:
%
% Tp
% ztop
% oxides
% parameters
% C_source_REE
%
% OUTPUT:
%
% C_solid:          Trace element composition within the solid (normalized)
% C_fluid:          Trace element composition within the liquid (normalized)
% C_extracted:      Trace element composition within the extracted liquid (normalized)
% C_ave:            Trace element composition within the total liquid in average (normalized)
%       Columns:    1    2    3    4    5    6    7    8    9    10   11   12   13
%                   La - Ce - Pr - Nd - Sm - Eu - Gd - Tb - Dy - Ho - Er - Yb - Lu
%
% C_solid_major:    Major element composition within the solid
% C_fluid_major:    Major element composition within the liquid
% C_extracted_major Major element composition within the extracted liquid
% C_ave_major:      Major element composition within the total liquid in average (normalized)
%       Columns:    1      2       3     4     5     6     7      8       9
%                   SiO2 - Al2O3 - FeO - MnO - MgO - CaO - Na2O - Cr2O3 - TiO2

% PTF:              Pressure (bars) - Temperature (Celsius) - F (mass percentage) path


%% file paths & names
working_folder = pwd;
specific_folder = 'perplex-model';
perplex_folder=sprintf('%s/%s',working_folder,specific_folder);

batch = 0;
flag_error = zeros(1,8);

tic
%% input variables
Tp       = Tp + 273.15;                                      % potential temperature [K]
isen     = parameters.model.isen;                            % which calculation
parallel = parameters.model.parallel;                       % how many parallel nodes

% if ztop<15e3; error('ztop must be greater than 15km'); end

%% type of computed isentrope
if isen == 1   
    
    Al2O3   = parameters.source.oxides.Al2O3;
    Na2O    = parameters.source.oxides.Na2O;
    
    Na2O_vector  = [0.01 0.05:0.05:1];
    [~,index_Na2O] = mink(abs(Na2O_vector-Na2O),2);
    val_Na2O = Na2O_vector(index_Na2O);
    Al2O3_vector = [0.2:0.2:4.6];
    [~,index_Al2O3] = mink(abs(Al2O3_vector-Al2O3),2);
    val_Al2O3 = Al2O3_vector(index_Al2O3);
    Tp_vector = [1200:5:1600];
    [~,index_Tp] = mink(abs(Tp_vector-(Tp-273.15)),2);
    val_Tp = Tp_vector(index_Tp);
    
    
    
    index_grid = [index_Na2O' index_Al2O3' index_Tp'];
    val_grid = [val_Na2O' val_Al2O3' val_Tp'];
    val_true = [Na2O Al2O3 Tp-273.15];
    
elseif isen == 2
    
    %%% source composition ... defined in melting_master
    SiO2    = parameters.source.oxides.SiO2;
    Al2O3   = parameters.source.oxides.Al2O3;
    Fe2O3   = parameters.source.oxides.Fe2O3;
    FeO     = parameters.source.oxides.FeOt;
    MnO     = parameters.source.oxides.MnO;
    MgO     = parameters.source.oxides.MgO;
    CaO     = parameters.source.oxides.CaO;
    Na2O    = parameters.source.oxides.Na2O;
    Cr2O3   = parameters.source.oxides.Cr2O3;
    TiO2    = parameters.source.oxides.TiO2;
    
    Na2O_grid     = parameters.source.oxides.grid.Na2O;
    Al2O3_grid    = parameters.source.oxides.grid.Al2O3;
    FeO_grid      = parameters.source.oxides.grid.FeO;
    MgO_grid      = parameters.source.oxides.grid.MgO;
    CaO_grid      = parameters.source.oxides.grid.CaO;
    Cr2O3_grid    = parameters.source.oxides.grid.Cr2O3;
    SiO2_grid     = parameters.source.oxides.grid.SiO2;
    Tp_grid       = parameters.source.oxides.grid.Tp;
    
    Na2O_stencil     = parameters.source.oxides.stencil.Na2O;
    Al2O3_stencil    = parameters.source.oxides.stencil.Al2O3;
    FeO_stencil      = parameters.source.oxides.stencil.FeO;
    MgO_stencil      = parameters.source.oxides.stencil.MgO;
    CaO_stencil      = parameters.source.oxides.stencil.CaO;
    Cr2O3_stencil    = parameters.source.oxides.stencil.Cr2O3;
    SiO2_stencil     = parameters.source.oxides.stencil.SiO2;
    Tp_stencil       = parameters.source.oxides.stencil.Tp;
    
    [~,index_Na2O_stencil]  = mink(abs(Na2O_stencil-Na2O),2);
    val_Na2O        = Na2O_stencil(index_Na2O_stencil);
    [~,index_Na2O_grid(1)]  = mink(abs(Na2O_grid-val_Na2O(1)),1);
    [~,index_Na2O_grid(2)]  = mink(abs(Na2O_grid-val_Na2O(2)),1);
    
    [~,index_Al2O3_stencil] = mink(abs(Al2O3_stencil-Al2O3),2);
    val_Al2O3       = Al2O3_stencil(index_Al2O3_stencil);
    [~,index_Al2O3_grid(1)] = mink(abs(Al2O3_grid-val_Al2O3(1)),1);
    [~,index_Al2O3_grid(2)] = mink(abs(Al2O3_grid-val_Al2O3(2)),1);
    
    [~,index_FeO_stencil]   = mink(abs(FeO_stencil-FeO),1);
    val_FeO         = FeO_stencil(index_FeO_stencil);
    [~,index_FeO_grid]   = mink(abs(FeO_grid-val_FeO(1)),1);
    
    [~,index_MgO_stencil]   = mink(abs(MgO_stencil-MgO),1);
    val_MgO         = MgO_stencil(index_MgO_stencil);
    [~,index_MgO_grid]   = mink(abs(MgO_grid-val_MgO(1)),1);
    
    [~,index_CaO_stencil]   = mink(abs(CaO_stencil-CaO),1);
    val_CaO         = CaO_stencil(index_CaO_stencil);
    [~,index_CaO_grid]   = mink(abs(CaO_grid-val_CaO(1)),1);
    
    [~,index_Cr2O3_stencil] = mink(abs(Cr2O3_stencil-Cr2O3),1);
    val_Cr2O3       = Cr2O3_stencil(index_Cr2O3_stencil);
    [~,index_Cr2O3_grid] = mink(abs(Cr2O3_grid-val_Cr2O3(1)),1);
    
    [~,index_SiO2_stencil]  = mink(abs(SiO2_stencil-SiO2),1);
    val_SiO2        = SiO2_stencil(index_SiO2_stencil);
    [~,index_SiO2_grid]  = mink(abs(SiO2_grid-val_SiO2(1)),1);
    
    [~,index_Tp_stencil]    = mink(abs(Tp_stencil-(Tp-273.15)),2);
    val_Tp          = Tp_stencil(index_Tp_stencil);
    [~,index_Tp_grid(1)]    = mink(abs(Tp_grid-val_Tp(1)),1);
    [~,index_Tp_grid(2)]    = mink(abs(Tp_grid-val_Tp(2)),1);
    
    
    index_grid = {index_Na2O_grid index_Al2O3_grid index_FeO_grid index_MgO_grid index_CaO_grid index_Cr2O3_grid index_SiO2_grid index_Tp_grid};
    val_grid = {val_Na2O val_Al2O3 val_FeO val_MgO val_CaO val_Cr2O3 val_SiO2 val_Tp};
    val_true = {Na2O Al2O3 FeO MgO CaO Cr2O3 SiO2 Tp};
    
else
    error('Contact the authors regarding other options to compute isentropes!')
    
end


%% Calculate melt compositions
tic
if isen == 1 || isen == 2
    [C_solid, C_fluid, C_extracted, C_ave, C_solid_major, C_fluid_major, C_extracted_major, C_ave_major, C_Ol, C_Cpx, C_Opx, C_Grt, C_Sp, C_Pl, w_j_cell, phi_j_cell, PTF_cell, M_out_smooth, F_out_smooth, PTF, grid_z, phi_j, v_s, v_f, trace_fail, flag_file, time_mat, flag_error] = ...
        run_nearest_diseqm(index_grid, val_grid, val_true,ztop,parameters,C_source_REE,isen,batch,flag_error,parallel,working_folder,perplex_folder);
end

% Check if trace element computation failed
if trace_fail==1
    disp('Error in computing traces - exit')
    flag_error(1) = 1;
    flag_error = single(flag_error>0);
    C_solid = [3003];       C_fluid = [3003];           C_extracted = [3003]; C_ave = [3003]; C_solid_major = [3003];
    C_fluid_major = [3003]; C_extracted_major = [3003]; C_ave_major = [3003]; PTF = [3003];
    cd(working_folder);
    return
end

%% Times 
if isen == 1
    verbouse_X = ['Isen = 1; Cube (',num2str(6-sum(flag_file)),'/6); Time computing solidus+isentropes: ',num2str(time_mat(1)),'s; Time smoothing isentrope: ',num2str(time_mat(2)),'s; Time computing trace elements: ',num2str(time_mat(3)),'s'];
    disp(verbouse_X)
elseif isen == 2
    verbouse_X = ['Isen = 2; Cube (',num2str(6-sum(flag_file)),'/6); Time computing solidus+isentropes: ',num2str(time_mat(1)),'s; Time smoothing isentrope: ',num2str(time_mat(2)),'s; Time computing trace elements: ',num2str(time_mat(3)),'s'];
    disp(verbouse_X)
elseif isen == 7
    verbouse_X = ['Isen = 7; Time computing solidus+isentrope: ',num2str(timeIsentrope_toc),'s; Time smoothing isentrope: ',num2str(timeSmooth_toc),'s; Time computing trace elements: ',num2str(timeTrace_toc),'s'];
    disp(verbouse_X)
end

end
