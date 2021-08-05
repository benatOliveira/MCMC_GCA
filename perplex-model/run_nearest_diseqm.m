function [C_solid_diseqm, C_fluid_diseqm, C_extracted_diseqm, C_ave, C_solid_major_diseqm, C_fluid_major_diseqm, C_extracted_major_diseqm, C_ave_major, C_Ol_diseqm, C_Cpx_diseqm, C_Opx_diseqm, C_Grt_diseqm, C_Sp_diseqm, C_Pl_diseqm, w_j_diseqm, phi_j_diseqm, PTF_cell, M_out, F_out, PTF, grid_z, phi_j, v_s, v_f, trace_fail, flag_file, time_mat, flag_error] = ...
    run_nearest_diseqm(index_grid, val_grid, val_true,ztop,parameters_new,C_source_REE,isen,batch,flag_error,parallel,working_folder,perplex_folder)



path_filename = pwd;
g       = parameters_new.model.g;                               % gravitational acceleration [m/s2]
rho_m   = parameters_new.model.rho_m;                           % mantle reference density [kg/m3]
    
% Predefine flag
flag_file = zeros(1,6);

%% type of isentrope
if isen == 1
    
    % Retrieve indexes & values
    index_Na2O = index_grid(:,1);
    index_Al2O3 = index_grid(:,2);
    index_Tp = index_grid(:,3);
    val_Na2O = val_grid(:,1);
    val_Al2O3 = val_grid(:,2);
    val_Tp = val_grid(:,3);
    
    index_Na2O_min  = min(index_Na2O);
    index_Na2O_max  = max(index_Na2O);
    index_Al2O3_min = min(index_Al2O3);
    index_Al2O3_max = max(index_Al2O3);
    index_Tp_min    = min(index_Tp);
    index_Tp_max    = max(index_Tp);
    val_Na2O_min  = min(val_Na2O);
    val_Na2O_max  = max(val_Na2O);
    val_Al2O3_min = min(val_Al2O3);
    val_Al2O3_max = max(val_Al2O3);
    val_Tp_min    = min(val_Tp);
    val_Tp_max    = max(val_Tp);
    
    % True values
    Na2O = val_true(:,1);
    Al2O3 = val_true(:,2);
    Tp = val_true(:,3);
    
    % Compute weights for linear interpolation
    DeltaNa2O = (Na2O-val_Na2O_min)/(val_Na2O_max-val_Na2O_min);
    DeltaAl2O3 = (Al2O3-val_Al2O3_min)/(val_Al2O3_max-val_Al2O3_min);
    DeltaTp = (Tp-val_Tp_min)/(val_Tp_max-val_Tp_min);
    
elseif isen == 2
    
    % Retrieve grid
    Na2O_grid     = parameters_new.source.oxides.grid.Na2O;
    Al2O3_grid    = parameters_new.source.oxides.grid.Al2O3;
    FeO_grid      = parameters_new.source.oxides.grid.FeO;
    MgO_grid      = parameters_new.source.oxides.grid.MgO;
    CaO_grid      = parameters_new.source.oxides.grid.CaO;
    Cr2O3_grid    = parameters_new.source.oxides.grid.Cr2O3;
    SiO2_grid     = parameters_new.source.oxides.grid.SiO2;
    Tp_grid       = parameters_new.source.oxides.grid.Tp;
    
    % Retrieve indexes & values
    index_Na2O = index_grid{1};
    index_Al2O3 = index_grid{2};
    index_FeO = index_grid{3};
    index_MgO = index_grid{4};
    index_CaO = index_grid{5};
    index_Cr2O3 = index_grid{6};
    index_SiO2 = index_grid{7};
    index_Tp = index_grid{8};
    val_Na2O = val_grid{1};
    val_Al2O3 = val_grid{2};
    val_FeO = val_grid{3};
    val_MgO = val_grid{4};
    val_CaO = val_grid{5};
    val_Cr2O3 = val_grid{6};
    val_SiO2 = val_grid{7};
    val_Tp = val_grid{8};
    
    % indexes in underlying grid
    index_Na2O_min  = min(index_Na2O);
    index_Na2O_max  = max(index_Na2O);
    index_Al2O3_min = min(index_Al2O3);
    index_Al2O3_max = max(index_Al2O3);
    index_FeO_min   = min(index_FeO);
    index_FeO_max   = max(index_FeO);
    index_MgO_min   = min(index_MgO);
    index_MgO_max   = max(index_MgO);
    index_CaO_min   = min(index_CaO);
    index_CaO_max   = max(index_CaO);
    index_Cr2O3_min = min(index_Cr2O3);
    index_Cr2O3_max = max(index_Cr2O3);
    index_SiO2_min  = min(index_SiO2);
    index_SiO2_max  = max(index_SiO2);
    index_Tp_min    = min(index_Tp);
    index_Tp_max    = max(index_Tp);
    val_Na2O_min    = min(val_Na2O);
    val_Na2O_max    = max(val_Na2O);
    val_Al2O3_min   = min(val_Al2O3);
    val_Al2O3_max   = max(val_Al2O3);
    val_FeO_min   = min(val_FeO);
    val_FeO_max   = max(val_FeO);
    val_MgO_min   = min(val_MgO);
    val_MgO_max   = max(val_MgO);
    val_CaO_min   = min(val_CaO);
    val_CaO_max   = max(val_CaO);
    val_Cr2O3_min   = min(val_Cr2O3);
    val_Cr2O3_max   = max(val_Cr2O3);
    val_SiO2_min   = min(val_SiO2);
    val_SiO2_max   = max(val_SiO2);
    val_Tp_min      = min(val_Tp);
    val_Tp_max      = max(val_Tp);
    
    % True values
    Na2O = val_true{1};
    Al2O3 = val_true{2};
    FeO = val_true{3};
    MgO = val_true{4};
    CaO = val_true{5};
    Cr2O3 = val_true{6};
    SiO2 = val_true{7};
    Tp = val_true{8};
    
    Fe2O3   = parameters_new.source.oxides.Fe2O3;
    MnO     = parameters_new.source.oxides.MnO;
    TiO2    = parameters_new.source.oxides.TiO2;
    
    % Compute weights for linear interpolation
    DeltaNa2O = (Na2O-val_Na2O_min)/(val_Na2O_max-val_Na2O_min);
    DeltaAl2O3 = (Al2O3-val_Al2O3_min)/(val_Al2O3_max-val_Al2O3_min);
    DeltaFeO = (FeO-val_FeO_min)/(val_FeO_max-val_FeO_min);
    DeltaMgO = (MgO-val_MgO_min)/(val_MgO_max-val_MgO_min);
    DeltaCaO = (CaO-val_CaO_min)/(val_CaO_max-val_CaO_min);
    DeltaCr2O3 = (Cr2O3-val_Cr2O3_min)/(val_Cr2O3_max-val_Cr2O3_min);
    DeltaSiO2 = (SiO2-val_SiO2_min)/(val_SiO2_max-val_SiO2_min);
    DeltaTp = ((Tp-273.15)-val_Tp_min)/(val_Tp_max-val_Tp_min);
    
    % Define 8 vertex of cube
    vertex_1 = [index_Na2O_min index_Al2O3_min index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_min];
    if DeltaNa2O>=DeltaAl2O3
        vertex_2 = [index_Na2O_max index_Al2O3_min index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_min];
    else
        vertex_2 = [index_Na2O_min index_Al2O3_max index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_min];
    end
    vertex_3 = [index_Na2O_max index_Al2O3_max index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_min];
    vertex_4 = [index_Na2O_min index_Al2O3_min index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_max];
    if DeltaNa2O>=DeltaAl2O3
        vertex_5 = [index_Na2O_max index_Al2O3_min index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_max];
    else
        vertex_5 = [index_Na2O_min index_Al2O3_max index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_max];
    end
    vertex_6 = [index_Na2O_max index_Al2O3_max index_FeO_min index_MgO_min index_CaO_min index_Cr2O3_min index_SiO2_min index_Tp_max];
    
    index_vertex = [vertex_1;
        vertex_2;
        vertex_3;
        vertex_4;
        vertex_5;
        vertex_6];%
    
    
    % Check if the isentrope file exists
    for index_filename = 1:size(index_vertex,1)
        %             path_filename = 'C:\Users\Benat\Dropbox\MCMC_EPSL\Code\version\MCMC_melting_27072020\perplex-model\Variables_Table_extra_aux';
        filename = strcat(path_filename,'/perplex-model/tables/input_trace_smooth_isentrope_Na2O_',num2str(index_vertex(index_filename,1)),'_Al2O3_',num2str(index_vertex(index_filename,2)),'_FeO_',num2str(index_vertex(index_filename,3)),'_MgO_',num2str(index_vertex(index_filename,4)),'_CaO_',num2str(index_vertex(index_filename,5)),'_Cr2O3_',num2str(index_vertex(index_filename,6)),'_SiO2_',num2str(index_vertex(index_filename,7)),'_Tp_',num2str(index_vertex(index_filename,8)),'.mat');
        
        if isfile(filename)
            S = load(filename);
            grid_z = S.grid_z;
            if ztop>=grid_z(end)
                flag_file(index_filename) = false;% File exists.
            else
                flag_file(index_filename) = true;% File exists but not long enough
            end
            clear grid_z
        else
            flag_file(index_filename) = true;% File does not exist.
        end
    end
    
    M_out_cell = cell(1,6);
    F_out_cell = cell(1,6);
    grid_z_cell = cell(1,6);
    PTF_cell = cell(1,6);
    
    % run 6 isentropes (only 6 are required for intterpolation)
    cd(perplex_folder);
    for index_isen = 1:6
        if flag_file(index_isen)  % doesnt exist
            [M_out, F_out, Error_out, grid_z, flag_error] = ...
                perplex_isentrope_opti(SiO2_grid(index_vertex(index_isen,7)),Al2O3_grid(index_vertex(index_isen,2)),Fe2O3,FeO_grid(index_vertex(index_isen,3)),MnO,MgO_grid(index_vertex(index_isen,4)),CaO_grid(index_vertex(index_isen,5)),Na2O_grid(index_vertex(index_isen,1)),Cr2O3_grid(index_vertex(index_isen,6)),TiO2,Tp_grid(index_vertex(index_isen,8)),ztop,parameters_new,batch,flag_error,parallel);
            
            if any(flag_error == 1) % define output if there is any error in isentrope computation
                C_solid_diseqm = [2002]; C_fluid_diseqm = [2002]; C_extracted_diseqm = [2002]; C_ave = [2002]; C_solid_major_diseqm = [2002]; C_fluid_major_diseqm = [2002]; C_extracted_major_diseqm = [2002]; C_ave_major = [2002];
                C_Ol_diseqm = [2002]; C_Cpx_diseqm = [2002]; C_Opx_diseqm = [2002]; C_Grt_diseqm = [2002]; C_Sp_diseqm = [2002]; C_Pl_diseqm = [2002];
                w_j_diseqm = [2002]; phi_j_diseqm = [2002]; phi_j = [2002]; v_s = [2002]; v_f = [2002]; PTF_cell = [2002];
                M_out = [2002]; F_out = [2002]; PTF = [2002]; grid_z = [2002];
                
                cd(working_folder);
                return
            end
            
            % pressure to grid
            grid_z = obtainGrid(M_out,rho_m,g);
            save_list = true;
        else    % exists
            filename = strcat(path_filename,'/perplex-model/tables/input_trace_smooth_isentrope_Na2O_',num2str(index_vertex(index_isen,1)),'_Al2O3_',num2str(index_vertex(index_isen,2)),'_FeO_',num2str(index_vertex(index_isen,3)),'_MgO_',num2str(index_vertex(index_isen,4)),'_CaO_',num2str(index_vertex(index_isen,5)),'_Cr2O3_',num2str(index_vertex(index_isen,6)),'_SiO2_',num2str(index_vertex(index_isen,7)),'_Tp_',num2str(index_vertex(index_isen,8)),'.mat');
            S = load(filename);
            M_out = S.M_out;
            F_out = S.F_out;
            Error_out = S.Error_out;
            grid_z = S.grid_z;
            flag_error = S.flag_error;
            dT = 0;
            n  = size(M_out,1);
            
            if ztop>=grid_z(end)
                save_list = false;
            else
                [M_out, F_out, Error_out, flag_error] = runIsentrope_Connolly_parallel(M_out,F_out,Error_out,P_top,dP,dT,perplex_folder,working_folder,batch,n,f0,parallel,flag_error);
                
                if any(flag_error == 1) % define output if there is any error in isentrope computation
                    C_solid_diseqm = [2002]; C_fluid_diseqm = [2002]; C_extracted_diseqm = [2002]; C_ave = [2002]; C_solid_major_diseqm = [2002]; C_fluid_major_diseqm = [2002]; C_extracted_major_diseqm = [2002]; C_ave_major = [2002];
                    C_Ol_diseqm = [2002]; C_Cpx_diseqm = [2002]; C_Opx_diseqm = [2002]; C_Grt_diseqm = [2002]; C_Sp_diseqm = [2002]; C_Pl_diseqm = [2002];
                    w_j_diseqm = [2002]; phi_j_diseqm = [2002]; phi_j = [2002]; v_s = [2002]; v_f = [2002]; PTF_cell = [2002];
                    M_out = [2002]; F_out = [2002]; PTF = [2002]; grid_z = [2002];
                    cd(working_folder);
                    return
                end
                % pressure to grid
                grid_z = obtainGrid(M_out,rho_m,g);
                save_list = true;
            end
        end
        
        if save_list  % save insentrope 
            filename = strcat(path_filename,'/perplex-model/tables/input_trace_smooth_isentrope_Na2O_',num2str(index_vertex(index_isen,1)),'_Al2O3_',num2str(index_vertex(index_isen,2)),'_FeO_',num2str(index_vertex(index_isen,3)),'_MgO_',num2str(index_vertex(index_isen,4)),'_CaO_',num2str(index_vertex(index_isen,5)),'_Cr2O3_',num2str(index_vertex(index_isen,6)),'_SiO2_',num2str(index_vertex(index_isen,7)),'_Tp_',num2str(index_vertex(index_isen,8)),'.mat');
            parsave(filename, single(M_out),single(F_out),single(Error_out),single(grid_z),single(flag_error));
        end
        
        [M_out_cell{index_isen},F_out_cell{index_isen},grid_z_cell{index_isen},PTF_cell{index_isen}] = correct_out(M_out,F_out,ztop,batch,parameters_new);
        
        
    end
end

cd(working_folder);

time_isentrope = toc;

%% Smooth and refine M_out, F_out, grid_z, PTF for each vertex

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isen == 1
    load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_min,index_Tp_min))
elseif isen == 2
    M_out_smooth = M_out_cell{1};
    F_out_smooth = F_out_cell{1};
    grid_z = grid_z_cell{1};
    PTF_smooth = PTF_cell{1};
end
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
% F out of smoothed data
F = 0*F_out_smooth(1,:);
F_out_smooth = 0*F_out_smooth;
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    F = F0 + (1-F0(1)).*(X-X0); 
    F_out_smooth(index,:) = F;
end

% Refine isentrope
[M_out_000,F_out_000,grid_z_000,PTF_000] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
%
parameters_000  = parameters_new;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DeltaNa2O >= DeltaAl2O3
    if isen == 1
        load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_min,index_Tp_min))
    elseif isen == 2
        M_out_smooth = M_out_cell{2};
        F_out_smooth = F_out_cell{2};
        grid_z = grid_z_cell{2};
        PTF_smooth = PTF_cell{2};
    end
    % Inerpolate melt compositions in nodes where there is no melt
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
    C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
    C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
    C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
    C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
    C_solid(isnan(C_solid)) = 0;
    C_melt(isnan(C_melt)) = 0;
    M_out_smooth(:,39:47) = C_solid*100;
    M_out_smooth(:,48:56) = C_melt*100;
    % F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
        F0 = F;
        X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
        X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
        F = F0 + (1-F0(1)).*(X-X0);
        F_out_smooth(index,:) = F;
    end
    
    % Refine isentrope
    [M_out_100,F_out_100,grid_z_100,PTF_100] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
    %
    parameters_100  = parameters_new;
else
    if isen == 1
        load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_max,index_Tp_min))
    elseif isen == 2
        M_out_smooth = M_out_cell{2};
        F_out_smooth = F_out_cell{2};
        grid_z = grid_z_cell{2};
        PTF_smooth = PTF_cell{2};
    end
    % Inerpolate melt compositions in nodes where there is no melt
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
    C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
    C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
    C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
    C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
    C_solid(isnan(C_solid)) = 0;
    C_melt(isnan(C_melt)) = 0;
    M_out_smooth(:,39:47) = C_solid*100;
    M_out_smooth(:,48:56) = C_melt*100;
    % F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
        F0 = F;
        X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
        X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
        F = F0 + (1-F0(1)).*(X-X0);
        F_out_smooth(index,:) = F;
    end
    
    % Refine isentrope
    [M_out_010,F_out_010,grid_z_010,PTF_010] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
    %
    parameters_010  = parameters_new;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third Vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isen == 1
    load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_max,index_Tp_min))
elseif isen == 2
    M_out_smooth = M_out_cell{3};
    F_out_smooth = F_out_cell{3};
    grid_z = grid_z_cell{3};
    PTF_smooth = PTF_cell{3};
end
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
% F out of smoothed data
F = 0*F_out_smooth(1,:);
F_out_smooth = 0*F_out_smooth;
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    F = F0 + (1-F0(1)).*(X-X0);
    F_out_smooth(index,:) = F;
end

% Refine isentrope
[M_out_110,F_out_110,grid_z_110,PTF_110] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
%
parameters_110  = parameters_new;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth Vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isen == 1
    load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_min,index_Tp_max))
elseif isen == 2
    M_out_smooth = M_out_cell{4};
    F_out_smooth = F_out_cell{4};
    grid_z = grid_z_cell{4};
    PTF_smooth = PTF_cell{4};
end
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
% F out of smoothed data
F = 0*F_out_smooth(1,:);
F_out_smooth = 0*F_out_smooth;
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    F = F0 + (1-F0(1)).*(X-X0); 
    F_out_smooth(index,:) = F;
end

% Refine isentrope
[M_out_001,F_out_001,grid_z_001,PTF_001] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
%
parameters_001  = parameters_new;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fifth Vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DeltaNa2O >= DeltaAl2O3
    if isen == 1
        load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_min,index_Tp_max))
    elseif isen == 2
        M_out_smooth = M_out_cell{5};
        F_out_smooth = F_out_cell{5};
        grid_z = grid_z_cell{5};
        PTF_smooth = PTF_cell{5};
    end
    % Inerpolate melt compositions in nodes where there is no melt
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
    C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
    C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
    C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
    C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
    C_solid(isnan(C_solid)) = 0;
    C_melt(isnan(C_melt)) = 0;
    M_out_smooth(:,39:47) = C_solid*100;
    M_out_smooth(:,48:56) = C_melt*100;
    % F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
        F0 = F;
        X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
        X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
        F = F0 + (1-F0(1)).*(X-X0);
        F_out_smooth(index,:) = F;
    end
    
    % Refine isentrope
    [M_out_101,F_out_101,grid_z_101,PTF_101] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
    %
    parameters_101  = parameters_new;
else
    if isen == 1
        load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_min,index_Al2O3_max,index_Tp_max))
    elseif isen == 2
        M_out_smooth = M_out_cell{5};
        F_out_smooth = F_out_cell{5};
        grid_z = grid_z_cell{5};
        PTF_smooth = PTF_cell{5};
    end
    % Inerpolate melt compositions in nodes where there is no melt
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
    C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
    C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
    C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
    C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
    C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
    C_solid(isnan(C_solid)) = 0;
    C_melt(isnan(C_melt)) = 0;
    M_out_smooth(:,39:47) = C_solid*100;
    M_out_smooth(:,48:56) = C_melt*100;
    % F out of smoothed data
    F = 0*F_out_smooth(1,:);
    F_out_smooth = 0*F_out_smooth;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
        F0 = F;
        X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
        X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
        F = F0 + (1-F0(1)).*(X-X0);
        F_out_smooth(index,:) = F;
    end
    
    % Refine isentrope
    [M_out_011,F_out_011,grid_z_011,PTF_011] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
    %
    parameters_011  = parameters_new;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sixth Vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isen == 1
    load(sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O_max,index_Al2O3_max,index_Tp_max))
elseif isen == 2
    M_out_smooth = M_out_cell{6};
    F_out_smooth = F_out_cell{6};
    grid_z = grid_z_cell{6};
    PTF_smooth = PTF_cell{6};
end
% Inerpolate melt compositions in nodes where there is no melt
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1); ind_o = ind(end);
C_solid_aux   = M_out_smooth(:,39:47)./repmat(sum(M_out_smooth(:,39:47),2),1,9); C_solid = C_solid_aux;
C_melt_aux    = M_out_smooth(:,48:56)./repmat(sum(M_out_smooth(:,48:56),2),1,9); C_melt  = C_melt_aux;
C_solid(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_solid_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_solid(1:ind_a,:) = C_solid_aux(1:ind_a,:);
C_melt(isnan(C_melt_aux(:,1)),:) = interp1(grid_z(~isnan(C_melt_aux(:,1))),C_melt_aux(~isnan(C_melt_aux(:,1)),:),grid_z(isnan(C_melt_aux(:,1))),'linear','extrap');
C_melt(1:ind_a,:) = C_melt_aux(1:ind_a,:);
C_solid(isnan(C_solid)) = 0;
C_melt(isnan(C_melt)) = 0;
M_out_smooth(:,39:47) = C_solid*100;
M_out_smooth(:,48:56) = C_melt*100;
% F out of smoothed data
F = 0*F_out_smooth(1,:);
F_out_smooth = 0*F_out_smooth;
ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
for index = ind_a:size(M_out_smooth,1)
    F0 = F;
    X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
    X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
    F = F0 + (1-F0(1)).*(X-X0); 
    F_out_smooth(index,:) = F;
end

% Refine isentrope
[M_out_111,F_out_111,grid_z_111,PTF_111] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters_new);
%
parameters_111  = parameters_new;

%% prepare the cube

if DeltaNa2O >= DeltaAl2O3
    M_out_cell = {M_out_000; M_out_100; M_out_110; M_out_001; M_out_101; M_out_111};
    F_out_cell = {F_out_000; F_out_100; F_out_110; F_out_001; F_out_101; F_out_111};
    grid_z_cell       = {grid_z_000; grid_z_100; grid_z_110; grid_z_001; grid_z_101; grid_z_111};
    parameters_cell   = {parameters_000; parameters_100; parameters_110; parameters_001; parameters_101; parameters_111};
    PTF_cell          = {PTF_000; PTF_100; PTF_110; PTF_001; PTF_101; PTF_111};
    
else
    M_out_cell = {M_out_000; M_out_010; M_out_110; M_out_001; M_out_011; M_out_111};
    F_out_cell = {F_out_000; F_out_010; F_out_110; F_out_001; F_out_011; F_out_111};
    grid_z_cell       = {grid_z_000; grid_z_010; grid_z_110; grid_z_001; grid_z_011; grid_z_111};
    parameters_cell   = {parameters_000; parameters_010; parameters_110; parameters_001; parameters_011; parameters_111};
    PTF_cell          = {PTF_000; PTF_010; PTF_110; PTF_001; PTF_011; PTF_111};
    
end

time_smooth = toc;

%% run trace element computation for each vertex of cube
trace_fail = 0;
tic
% parfor index_vertex = 1:6
for index_vertex = 1:6
    if sum(F_out_cell{index_vertex}(:,1))>0 % there is melt in the column
        [C_solid_vertex, C_fluid_vertex, C_extracted_vertex, C_ave_vertex, ...
            C_solid_major_vertex, C_fluid_major_vertex, C_extracted_major_vertex, C_ave_major_vertex, ...
            C_Ol_vertex, C_Cpx_vertex, C_Opx_vertex, C_Grt_vertex, C_Sp_vertex, C_Pl_vertex, ...
            w_j_vertex, phi_j_vertex, v_s_vertex, v_f_vertex, trace_fail_vertex] = diseqm_trace_parallel_fractional(M_out_cell{index_vertex},F_out_cell{index_vertex},grid_z_cell{index_vertex},parameters_cell{index_vertex},C_source_REE);
        
        flag_melt{index_vertex} = 1;
    else
        [C_solid_vertex, C_fluid_vertex, C_extracted_vertex, C_ave_vertex, ...
            C_solid_major_vertex, C_fluid_major_vertex, C_extracted_major_vertex, C_ave_major_vertex, ...
            C_Ol_vertex, C_Cpx_vertex, C_Opx_vertex, C_Grt_vertex, C_Sp_vertex, C_Pl_vertex, ...
            w_j_vertex, phi_j_vertex, v_s_vertex, v_f_vertex, trace_fail_vertex] = deal(1001);
        
        flag_melt{index_vertex} = 0;
    end
    
    % Assign output
    C_solid_diseqm{index_vertex} = C_solid_vertex;
    C_fluid_diseqm{index_vertex} = C_fluid_vertex;
    C_extracted_diseqm{index_vertex} = C_extracted_vertex;
    C_ave_diseqm{index_vertex} = C_ave_vertex;
    C_solid_major_diseqm{index_vertex} = C_solid_major_vertex;
    C_fluid_major_diseqm{index_vertex} = C_fluid_major_vertex;
    C_extracted_major_diseqm{index_vertex} = C_extracted_major_vertex;
    C_ave_major_diseqm{index_vertex} = C_ave_major_vertex;
    C_Ol_diseqm{index_vertex} = C_Ol_vertex;
    C_Cpx_diseqm{index_vertex} = C_Cpx_vertex;
    C_Opx_diseqm{index_vertex} = C_Opx_vertex;
    C_Grt_diseqm{index_vertex} = C_Grt_vertex;
    C_Sp_diseqm{index_vertex} = C_Sp_vertex;
    C_Pl_diseqm{index_vertex} = C_Pl_vertex;
    w_j_diseqm{index_vertex} = w_j_vertex;
    phi_j_diseqm{index_vertex} = phi_j_vertex;
    v_s_diseqm{index_vertex} = v_s_vertex;
    v_f_diseqm{index_vertex} = v_f_vertex;
    trace_fail_diseqm{index_vertex} = trace_fail_vertex;
    
end

% chek if trace element computation failed
if any(cell2mat(trace_fail_diseqm)==1)
    trace_fail = 1;
end

%% interpolate
if sum(cell2mat(flag_melt))==6
    if DeltaNa2O >= DeltaAl2O3
        %     M_out_cell = {M_out_000; M_out_100; M_out_110; M_out_001; M_out_101; M_out_111};
        %                           1          2          3          4          5          6 ;  number of the vertex of the cube
        C_ave               = C_ave_diseqm{1}       + (C_ave_diseqm{2}-C_ave_diseqm{1})*DeltaNa2O + (C_ave_diseqm{3}-C_ave_diseqm{2})*DeltaAl2O3                            + (C_ave_diseqm{4}-C_ave_diseqm{1})*DeltaTp + ...
            (C_ave_diseqm{5}-C_ave_diseqm{4}-C_ave_diseqm{2}+C_ave_diseqm{1})*DeltaNa2O*DeltaTp                                     + (C_ave_diseqm{6}-C_ave_diseqm{5}-C_ave_diseqm{3}+C_ave_diseqm{2})*DeltaAl2O3*DeltaTp;
        C_ave_major         = C_ave_major_diseqm{1} + (C_ave_major_diseqm{2}-C_ave_major_diseqm{1})*DeltaNa2O + (C_ave_major_diseqm{3}-C_ave_major_diseqm{2})*DeltaAl2O3    + (C_ave_major_diseqm{4}-C_ave_major_diseqm{1})*DeltaTp + ...
            (C_ave_major_diseqm{5}-C_ave_major_diseqm{4}-C_ave_major_diseqm{2}+C_ave_major_diseqm{1})*DeltaNa2O*DeltaTp             + (C_ave_major_diseqm{6}-C_ave_major_diseqm{5}-C_ave_major_diseqm{3}+C_ave_major_diseqm{2})*DeltaAl2O3*DeltaTp;
    else
        %     M_out_cell = {M_out_000; M_out_010; M_out_110; M_out_001; M_out_011; M_out_111};
        %                           1          2          3          4          5          6 ;  number of the vertex of the cube
        C_ave               = C_ave_diseqm{1}       + (C_ave_diseqm{3}-C_ave_diseqm{2})*DeltaNa2O + (C_ave_diseqm{2}-C_ave_diseqm{1})*DeltaAl2O3                            + (C_ave_diseqm{4}-C_ave_diseqm{1})*DeltaTp + ...
            (C_ave_diseqm{6}-C_ave_diseqm{5}-C_ave_diseqm{3}+C_ave_diseqm{2})*DeltaNa2O*DeltaTp                                     + (C_ave_diseqm{5}-C_ave_diseqm{4}-C_ave_diseqm{2}+C_ave_diseqm{1})*DeltaAl2O3*DeltaTp;
        C_ave_major         = C_ave_major_diseqm{1} + (C_ave_major_diseqm{3}-C_ave_major_diseqm{2})*DeltaNa2O + (C_ave_major_diseqm{2}-C_ave_major_diseqm{1})*DeltaAl2O3    + (C_ave_major_diseqm{4}-C_ave_major_diseqm{1})*DeltaTp + ...
            (C_ave_major_diseqm{6}-C_ave_major_diseqm{5}-C_ave_major_diseqm{3}+C_ave_major_diseqm{2})*DeltaNa2O*DeltaTp             + (C_ave_major_diseqm{5}-C_ave_major_diseqm{4}-C_ave_major_diseqm{2}+C_ave_major_diseqm{1})*DeltaAl2O3*DeltaTp;
    end
else
    % Define output when there is no melt
    C_ave               = 1001;
    C_ave_major         = 1001;
end

time_trace = toc;

%% times
time_mat = [time_isentrope time_smooth time_trace];


%% extract interpolate isentrope (optional)

if sum(cell2mat(flag_melt))==6
    length_000      = size(F_out_000,1);
    length_110      = size(F_out_110,1);
    length_001      = size(F_out_001,1);
    length_111      = size(F_out_111,1);
    
    if DeltaNa2O >= DeltaAl2O3
        length_100      = size(F_out_100,1);
        length_101      = size(F_out_101,1);
        max_length = max([length_000 length_100 length_110 length_001 length_101 length_111]) ;
        min_length = min([length_000 length_100 length_110 length_001 length_101 length_111]) ;
        if max_length~=length_000; F_out_000 = [repmat(F_out_000(1,:),max_length-length_000,1); F_out_000]; end
        if max_length~=length_100; F_out_100 = [repmat(F_out_100(1,:),max_length-length_100,1); F_out_100]; end
        if max_length~=length_110; F_out_110 = [repmat(F_out_110(1,:),max_length-length_110,1); F_out_110]; end
        if max_length~=length_001; F_out_001 = [repmat(F_out_001(1,:),max_length-length_001,1); F_out_001]; end
        if max_length~=length_101; F_out_101 = [repmat(F_out_101(1,:),max_length-length_101,1); F_out_101]; end
        if max_length~=length_111; F_out_111 = [repmat(F_out_111(1,:),max_length-length_111,1); F_out_111]; end
        if max_length~=length_000; M_out_000 = [repmat(M_out_000(1,:),max_length-length_000,1); M_out_000]; end
        if max_length~=length_100; M_out_100 = [repmat(M_out_100(1,:),max_length-length_100,1); M_out_100]; end
        if max_length~=length_110; M_out_110 = [repmat(M_out_110(1,:),max_length-length_110,1); M_out_110]; end
        if max_length~=length_001; M_out_001 = [repmat(M_out_001(1,:),max_length-length_001,1); M_out_001]; end
        if max_length~=length_101; M_out_101 = [repmat(M_out_101(1,:),max_length-length_101,1); M_out_101]; end
        if max_length~=length_111; M_out_111 = [repmat(M_out_111(1,:),max_length-length_111,1); M_out_111]; end
        if max_length~=length_000; phi_j_000 = [repmat(phi_j_diseqm{1}(1,:),max_length-length_000,1); phi_j_diseqm{1}]; else; phi_j_000 = phi_j_diseqm{1}; end
        if max_length~=length_100; phi_j_100 = [repmat(phi_j_diseqm{2}(1,:),max_length-length_100,1); phi_j_diseqm{2}]; else; phi_j_100 = phi_j_diseqm{2}; end
        if max_length~=length_110; phi_j_110 = [repmat(phi_j_diseqm{3}(1,:),max_length-length_110,1); phi_j_diseqm{3}]; else; phi_j_110 = phi_j_diseqm{3}; end
        if max_length~=length_001; phi_j_001 = [repmat(phi_j_diseqm{4}(1,:),max_length-length_001,1); phi_j_diseqm{4}]; else; phi_j_001 = phi_j_diseqm{4}; end
        if max_length~=length_101; phi_j_101 = [repmat(phi_j_diseqm{5}(1,:),max_length-length_101,1); phi_j_diseqm{5}]; else; phi_j_101 = phi_j_diseqm{5}; end
        if max_length~=length_111; phi_j_111 = [repmat(phi_j_diseqm{6}(1,:),max_length-length_111,1); phi_j_diseqm{6}]; else; phi_j_111 = phi_j_diseqm{6}; end
        if max_length~=length_000; v_s_000 = [repmat(v_s_diseqm{1}(1,:),max_length-length_000,1); v_s_diseqm{1}]; else; v_s_000 = v_s_diseqm{1}; end
        if max_length~=length_100; v_s_100 = [repmat(v_s_diseqm{2}(1,:),max_length-length_100,1); v_s_diseqm{2}]; else; v_s_100 = v_s_diseqm{2}; end
        if max_length~=length_110; v_s_110 = [repmat(v_s_diseqm{3}(1,:),max_length-length_110,1); v_s_diseqm{3}]; else; v_s_110 = v_s_diseqm{3}; end
        if max_length~=length_001; v_s_001 = [repmat(v_s_diseqm{4}(1,:),max_length-length_001,1); v_s_diseqm{4}]; else; v_s_001 = v_s_diseqm{4}; end
        if max_length~=length_101; v_s_101 = [repmat(v_s_diseqm{5}(1,:),max_length-length_101,1); v_s_diseqm{5}]; else; v_s_101 = v_s_diseqm{5}; end
        if max_length~=length_111; v_s_111 = [repmat(v_s_diseqm{6}(1,:),max_length-length_111,1); v_s_diseqm{6}]; else; v_s_111 = v_s_diseqm{6}; end
        if max_length~=length_000; v_f_000 = [repmat(v_f_diseqm{1}(1,:),max_length-length_000,1); v_f_diseqm{1}]; else; v_f_000 = v_f_diseqm{1}; end
        if max_length~=length_100; v_f_100 = [repmat(v_f_diseqm{2}(1,:),max_length-length_100,1); v_f_diseqm{2}]; else; v_f_100 = v_f_diseqm{2}; end
        if max_length~=length_110; v_f_110 = [repmat(v_f_diseqm{3}(1,:),max_length-length_110,1); v_f_diseqm{3}]; else; v_f_110 = v_f_diseqm{3}; end
        if max_length~=length_001; v_f_001 = [repmat(v_f_diseqm{4}(1,:),max_length-length_001,1); v_f_diseqm{4}]; else; v_f_001 = v_f_diseqm{4}; end
        if max_length~=length_101; v_f_101 = [repmat(v_f_diseqm{5}(1,:),max_length-length_101,1); v_f_diseqm{5}]; else; v_f_101 = v_f_diseqm{5}; end
        if max_length~=length_111; v_f_111 = [repmat(v_f_diseqm{6}(1,:),max_length-length_111,1); v_f_diseqm{6}]; else; v_f_111 = v_f_diseqm{6}; end
    else
        length_010      = size(F_out_010,1);
        length_011      = size(F_out_011,1);
        max_length = max([length_000 length_010 length_110 length_001 length_011 length_111]) ;
        min_length = min([length_000 length_010 length_110 length_001 length_011 length_111]) ;
        if max_length~=length_000; F_out_000 = [repmat(F_out_000(1,:),max_length-length_000,1); F_out_000]; end
        if max_length~=length_010; F_out_010 = [repmat(F_out_010(1,:),max_length-length_010,1); F_out_010]; end
        if max_length~=length_110; F_out_110 = [repmat(F_out_110(1,:),max_length-length_110,1); F_out_110]; end
        if max_length~=length_001; F_out_001 = [repmat(F_out_001(1,:),max_length-length_001,1); F_out_001]; end
        if max_length~=length_011; F_out_011 = [repmat(F_out_011(1,:),max_length-length_011,1); F_out_011]; end
        if max_length~=length_111; F_out_111 = [repmat(F_out_111(1,:),max_length-length_111,1); F_out_111]; end
        if max_length~=length_000; M_out_000 = [repmat(M_out_000(1,:),max_length-length_000,1); M_out_000]; end
        if max_length~=length_010; M_out_010 = [repmat(M_out_010(1,:),max_length-length_010,1); M_out_010]; end
        if max_length~=length_110; M_out_110 = [repmat(M_out_110(1,:),max_length-length_110,1); M_out_110]; end
        if max_length~=length_001; M_out_001 = [repmat(M_out_001(1,:),max_length-length_001,1); M_out_001]; end
        if max_length~=length_011; M_out_011 = [repmat(M_out_011(1,:),max_length-length_011,1); M_out_011]; end
        if max_length~=length_111; M_out_111 = [repmat(M_out_111(1,:),max_length-length_111,1); M_out_111]; end
        if max_length~=length_000; phi_j_000 = [repmat(phi_j_diseqm{1}(1,:),max_length-length_000,1); phi_j_diseqm{1}]; else; phi_j_000 = phi_j_diseqm{1}; end
        if max_length~=length_010; phi_j_010 = [repmat(phi_j_diseqm{2}(1,:),max_length-length_010,1); phi_j_diseqm{2}]; else; phi_j_010 = phi_j_diseqm{2}; end
        if max_length~=length_110; phi_j_110 = [repmat(phi_j_diseqm{3}(1,:),max_length-length_110,1); phi_j_diseqm{3}]; else; phi_j_110 = phi_j_diseqm{3}; end
        if max_length~=length_001; phi_j_001 = [repmat(phi_j_diseqm{4}(1,:),max_length-length_001,1); phi_j_diseqm{4}]; else; phi_j_001 = phi_j_diseqm{4}; end
        if max_length~=length_011; phi_j_011 = [repmat(phi_j_diseqm{5}(1,:),max_length-length_011,1); phi_j_diseqm{5}]; else; phi_j_011 = phi_j_diseqm{5}; end
        if max_length~=length_111; phi_j_111 = [repmat(phi_j_diseqm{6}(1,:),max_length-length_111,1); phi_j_diseqm{6}]; else; phi_j_111 = phi_j_diseqm{6}; end
        if max_length~=length_000; v_s_000 = [repmat(v_s_diseqm{1}(1,:),max_length-length_000,1); v_s_diseqm{1}]; else; v_s_000 = v_s_diseqm{1}; end
        if max_length~=length_010; v_s_010 = [repmat(v_s_diseqm{2}(1,:),max_length-length_010,1); v_s_diseqm{2}]; else; v_s_010 = v_s_diseqm{2}; end
        if max_length~=length_110; v_s_110 = [repmat(v_s_diseqm{3}(1,:),max_length-length_110,1); v_s_diseqm{3}]; else; v_s_110 = v_s_diseqm{3}; end
        if max_length~=length_001; v_s_001 = [repmat(v_s_diseqm{4}(1,:),max_length-length_001,1); v_s_diseqm{4}]; else; v_s_001 = v_s_diseqm{4}; end
        if max_length~=length_011; v_s_011 = [repmat(v_s_diseqm{5}(1,:),max_length-length_011,1); v_s_diseqm{5}]; else; v_s_011 = v_s_diseqm{5}; end
        if max_length~=length_111; v_s_111 = [repmat(v_s_diseqm{6}(1,:),max_length-length_111,1); v_s_diseqm{6}]; else; v_s_111 = v_s_diseqm{6}; end
        if max_length~=length_000; v_f_000 = [repmat(v_f_diseqm{1}(1,:),max_length-length_000,1); v_f_diseqm{1}]; else; v_f_000 = v_f_diseqm{1}; end
        if max_length~=length_010; v_f_010 = [repmat(v_f_diseqm{2}(1,:),max_length-length_010,1); v_f_diseqm{2}]; else; v_f_010 = v_f_diseqm{2}; end
        if max_length~=length_110; v_f_110 = [repmat(v_f_diseqm{3}(1,:),max_length-length_110,1); v_f_diseqm{3}]; else; v_f_110 = v_f_diseqm{3}; end
        if max_length~=length_001; v_f_001 = [repmat(v_f_diseqm{4}(1,:),max_length-length_001,1); v_f_diseqm{4}]; else; v_f_001 = v_f_diseqm{4}; end
        if max_length~=length_011; v_f_011 = [repmat(v_f_diseqm{5}(1,:),max_length-length_011,1); v_f_diseqm{5}]; else; v_f_011 = v_f_diseqm{5}; end
        if max_length~=length_111; v_f_111 = [repmat(v_f_diseqm{6}(1,:),max_length-length_111,1); v_f_diseqm{6}]; else; v_f_111 = v_f_diseqm{6}; end
    end
    
    if DeltaNa2O >= DeltaAl2O3
        F_out = F_out_000 + (F_out_100-F_out_000)*DeltaNa2O + (F_out_110-F_out_100)*DeltaAl2O3 + (F_out_001-F_out_000)*DeltaTp + ...
            (F_out_101-F_out_001-F_out_100+F_out_000)*DeltaNa2O*DeltaTp + (F_out_111-F_out_101-F_out_110+F_out_100)*DeltaAl2O3*DeltaTp;
        
        M_out = M_out_000 + (M_out_100-M_out_000)*DeltaNa2O + (M_out_110-M_out_100)*DeltaAl2O3 + (M_out_001-M_out_000)*DeltaTp + ...
            (M_out_101-M_out_001-M_out_100+M_out_000)*DeltaNa2O*DeltaTp + (M_out_111-M_out_101-M_out_110+M_out_100)*DeltaAl2O3*DeltaTp;
        
        phi_j = phi_j_000 + (phi_j_100-phi_j_000)*DeltaNa2O + (phi_j_110-phi_j_100)*DeltaAl2O3 + (phi_j_001-phi_j_000)*DeltaTp + ...
            (phi_j_101-phi_j_001-phi_j_100+phi_j_000)*DeltaNa2O*DeltaTp + (phi_j_111-phi_j_101-phi_j_110+phi_j_100)*DeltaAl2O3*DeltaTp;
        
        v_s = v_s_000 + (v_s_100-v_s_000)*DeltaNa2O + (v_s_110-v_s_100)*DeltaAl2O3 + (v_s_001-v_s_000)*DeltaTp + ...
            (v_s_101-v_s_001-v_s_100+v_s_000)*DeltaNa2O*DeltaTp + (v_s_111-v_s_101-v_s_110+v_s_100)*DeltaAl2O3*DeltaTp;
        
        v_f = v_f_000 + (v_f_100-v_f_000)*DeltaNa2O + (v_f_110-v_f_100)*DeltaAl2O3 + (v_f_001-v_f_000)*DeltaTp + ...
            (v_f_101-v_f_001-v_f_100+v_f_000)*DeltaNa2O*DeltaTp + (v_f_111-v_f_101-v_f_110+v_f_100)*DeltaAl2O3*DeltaTp;
    else
        F_out = F_out_000 + (F_out_110-F_out_010)*DeltaNa2O + (F_out_010-F_out_000)*DeltaAl2O3 + (F_out_001-F_out_000)*DeltaTp + ...
            (F_out_111-F_out_011-F_out_110+F_out_010)*DeltaNa2O*DeltaTp + (F_out_011-F_out_001-F_out_010+F_out_000)*DeltaAl2O3*DeltaTp;
        
        M_out = M_out_000 + (M_out_110-M_out_010)*DeltaNa2O + (M_out_010-M_out_000)*DeltaAl2O3 + (M_out_001-M_out_000)*DeltaTp + ...
            (M_out_111-M_out_011-M_out_110+M_out_010)*DeltaNa2O*DeltaTp + (M_out_011-M_out_001-M_out_010+M_out_000)*DeltaAl2O3*DeltaTp;
        
        phi_j = phi_j_000 + (phi_j_110-phi_j_010)*DeltaNa2O + (phi_j_010-phi_j_000)*DeltaAl2O3 + (phi_j_001-phi_j_000)*DeltaTp + ...
            (phi_j_111-phi_j_011-phi_j_110+phi_j_010)*DeltaNa2O*DeltaTp + (phi_j_011-phi_j_001-phi_j_010+phi_j_000)*DeltaAl2O3*DeltaTp;
        
        v_s = v_s_000 + (v_s_110-v_s_010)*DeltaNa2O + (v_s_010-v_s_000)*DeltaAl2O3 + (v_s_001-v_s_000)*DeltaTp + ...
            (v_s_111-v_s_011-v_s_110+v_s_010)*DeltaNa2O*DeltaTp + (v_s_011-v_s_001-v_s_010+v_s_000)*DeltaAl2O3*DeltaTp;
        
        v_f = v_f_000 + (v_f_110-v_f_010)*DeltaNa2O + (v_f_010-v_f_000)*DeltaAl2O3 + (v_f_001-v_f_000)*DeltaTp + ...
            (v_f_111-v_f_011-v_f_110+v_f_010)*DeltaNa2O*DeltaTp + (v_f_011-v_f_001-v_f_010+v_f_000)*DeltaAl2O3*DeltaTp;
    end
    
    % Create grid_z
    P_grid = [M_out(:,1); 0]*1e5;
    rho_s  = M_out(end:-1:1,10);
    rho_s(1) = parameters_new.model.rho_m;
    
    g = 9.81;
    DeltaP = P_grid(end-1:-1:1)-P_grid(end:-1:2);
    dh = DeltaP./(g*rho_s);
    
    grid_z = flipud(cumsum(dh));
    
    
    [ind_zTop, ~] = find(ztop>grid_z);
    if length(ind_zTop)==length(grid_z)
        M_out = [M_out(1,:); M_out];
        F_out = [F_out(1,:); F_out];
        grid_z = [ztop; grid_z];
        [ind_zTop, ~] = find(ztop>grid_z);
    end
    if ~isempty(ind_zTop)
        DeltaZ = grid_z(ind_zTop(1)-1)-grid_z(ind_zTop(1));
        weightA = abs(grid_z(ind_zTop(1)-1)-ztop)/DeltaZ;
        weightB = abs(grid_z(ind_zTop(1))-ztop)/DeltaZ;
        % Interpolate for depth
        M_out(ind_zTop(1),:) = M_out(ind_zTop(1)-1,:)*weightA + M_out(ind_zTop(1),:)*weightB;
        F_out(ind_zTop(1),:) = F_out(ind_zTop(1)-1,:)*weightA + F_out(ind_zTop(1),:)*weightB;
        grid_z(ind_zTop(1),:) = ztop;
        
        % Remove shallower nodes
        if length(ind_zTop)>1
            M_out(ind_zTop(2):end,:) = [];
            F_out(ind_zTop(2):end,:) = [];
            grid_z(ind_zTop(2):end) = [];
        end
    end
    
    % Remove first nodes (deep nodes) with only solid
    ind    = find(M_out(:,12)>0);
    if isempty(ind)
        ind = length(M_out(:,12));
        %     keyboard
    end
    ind_a  = ind(1);
    M_out(1:ind_a-2,:)  = [];
    F_out(1:ind_a-2,:)  = [];
    grid_z(1:ind_a-2,:) = [];
    
    
    PTF    = [M_out(:,1),M_out(:,2)-273.15,M_out(:,11)];
    
else
    % Define output when there is no melt
    M_out               = 1001;
    F_out               = 1001;
    PTF                 = 1001;
    grid_z              = 1001;
    phi_j               = 1001;
    v_s                 = 1001;
    v_f                 = 1001;
    
end


