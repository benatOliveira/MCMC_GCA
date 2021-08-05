% function [C_solid_diseqm, C_fluid_diseqm, C_extracted_diseqm, C_ave, C_solid_major_diseqm, C_fluid_major_diseqm, C_extracted_major_diseqm, C_ave_major, C_Ol_diseqm, C_Cpx_diseqm, C_Opx_diseqm, C_Grt_diseqm, C_Sp_diseqm, C_Pl_diseqm, w_j_diseqm, phi_j_diseqm, PTF_cell, M_out, F_out, PTF, grid_z, phi_j, v_s, v_f] = run_nearest_diseqm_adaptative(index_table,val_table,Tp,ztop,parameters,C_source_REE)
function [M_out,F_out,Error_out,grid_z,flag_error] = perplex_isentrope_opti(SiO2,Al2O3,Fe2O3,FeO,MnO,MgO,CaO,Na2O,Cr2O3,TiO2,Tp,ztop,parameters,batch,flag_error,parallel)

working_folder = pwd;
specific_folder = 'perplex-model';
perplex_folder=working_folder;

working_folder = pwd;
specific_folder_1 = 'mex_solid';
specific_folder_2 = 'mex_melt';
perplex_folder_1=sprintf('%s/%s',working_folder,specific_folder_1);
perplex_folder_2=sprintf('%s/%s',working_folder,specific_folder_2);


M_out = zeros(1,74);
F_out =  zeros(1,7);
Error_out = zeros(1,2);
grid_z = ztop;

%% input variables
Tp      = Tp + 273.15;                                      % change to K
f0      = parameters.model.fM;                               % residual porosity [%]
dP      = parameters.model.dP;
g       = parameters.model.g;                               % gravitational acceleration [m/s2]
rho_m   = parameters.model.rho_m;                           % mantle reference density [kg/m3]
alpha   = parameters.model.alpha_m;                         % thermal expansion coefficient [1/K]
Cp      = parameters.model.Cp_m;                            % specific heat [J/(kg K)]
dTdP = (Tp * alpha * 1e5) / (rho_m * Cp);                   % adiabatic gradient [K/bar] = [1e5 K/Pa]
P_top = rho_m * g * ztop / 1e5;                             % final pressure of melting [bar]
P_top = rho_m * g * 15000 / 1e5;                            % final pressure of melting [bar]

% FeOt = FeO + 0.9*Fe2O3;
FeOt = FeO;

%% Compute the first shallow point to obtain S_0
P_shallow = 2e3;  % bars
T_init = Tp+dTdP*P_shallow;                        % initial temperature [K]


%%% prepare input for Perplex;
input_source = [SiO2 Al2O3 FeOt MnO MgO CaO Na2O Cr2O3 TiO2];
% input_source = [46.7454    3.4948    7.5432         0   38.0709    3.3273    0.3135    0.5050         0];
input_complete = [P_shallow T_init input_source];
% input_complete = 1.0e+03 *[2.0000    1.6660    0.0467    0.0035    0.0075         0    0.0381    0.0033    0.0003    0.0005         0]

v = input_complete(1:2);
cblk = input_complete(3:end);
% v = [P_shallow T_init]; 
% cblk = input_source;

cd(perplex_folder_1)
[a,b,c,d,e]=meemum_fun_solid_gateway(v,cblk);
read_mexFile;

aux_cpx = 0;
while isnan(Cp_System)
    if aux_cpx == 5
        disp('Error computing S_0 in shallowest point: T has been unsuccessfully incremented to avoid NaN - exit')
        flag_error(1) = 1;
        C_solid = [2002];       C_fluid = [2002];           C_extracted = [2002]; C_ave = [2002]; C_solid_major = [2002];
        C_fluid_major = [2002]; C_extracted_major = [2002]; C_ave_major = [2002]; PTF = [2002];
        cd(working_folder);
        return
    end
    input_complete(2) = input_complete(2)+1;    
    
    v = input_complete(1:2);
    cblk = input_complete(3:end);
    
    [a,b,c,d,e]=meemum_fun_solid_gateway(v,cblk);
    read_mexFile;
    
    aux_cpx = aux_cpx + 1;
end
S_System_0 = S_Solid;
% cd(working_folder);

%% Obtain Temperature at given pressure (at depth) with S_0

% solidus after Hirschman 2000
%         T_sol = a*P^2 + b*P + c;
a=-5.104;
b=132.899;
c=1120.661 + 273.15;  % K

% solve system of equations for P, neglecting quadratic term
P_init = (c-Tp)/(dTdP-b/1e4);   % bar
% Round
P_init = abs(dP)*ceil(P_init/abs(dP));

adiab_geo = alpha_Solid/(rho_Solid*Cp_Solid*1e-5);

exit_aux_P = 0;
d_dT = 5;
dT_vec_aux_0 = 0:d_dT:(parallel-1)*d_dT;
dT_vec_aux = (dT_vec_aux_0-dT_vec_aux_0(round(parallel/2)));

%         dT_vec_aux = [50 25 0 -25 -50 -60];

cd(perplex_folder_2);

while exit_aux_P == 0
    
    % Add some buffer
    P_init = P_init + 5000;
    
    exit_aux_T = 0;
    while exit_aux_T == 0 % && exit_aux_P == 0
        dT = adiab_geo*Tp *P_init;        % adiabatic gradient        [K]; all pressure in [bar]
        
        %        dT_vec = dT_closest + [4 3 2 1 0 -1 -2 -3 -4 -5];
        dT_vec = dT + dT_vec_aux;
        T_init_cell = num2cell(Tp + dT_vec);
        
        %         tic
        parfor index_par = 1:parallel
%                                 for index_par = 1:vertex
            
            %%% Gibbs free energy minimisation (with updated T, P, C) %%%
            
            input_complete = double([P_init T_init_cell{index_par} input_source]);
            
            [P_par{index_par}, T_par{index_par}, S_System_par{index_par}, S_Solid_par{index_par}, S_Melt_par{index_par}, V_Solid_par{index_par}, V_Melt_par{index_par}, Cp_System_par{index_par}, Cp_Solid_par{index_par}, alpha_Solid_par{index_par}, rho_Solid_par{index_par},...
                X_Melt_par{index_par}, X_Ol_par{index_par}, X_Cpx_par{index_par}, X_Opx_par{index_par}, X_Gt_par{index_par}, X_Sp_par{index_par}, X_Pl_par{index_par}, Vol_Melt_par{index_par}, Vol_Ol_par{index_par}, Vol_Cpx_par{index_par}, Vol_Opx_par{index_par}, Vol_Gt_par{index_par}, Vol_Sp_par{index_par}, Vol_Pl_par{index_par},...
                Rho_Melt_par{index_par}, Rho_Ol_par{index_par}, Rho_Cpx_par{index_par}, Rho_Opx_par{index_par}, Rho_Gt_par{index_par}, Rho_Sp_par{index_par}, Rho_Pl_par{index_par},...
                H_Solid_par{index_par}, H_Melt_par{index_par}, H_Gt_par{index_par}, H_Sp_par{index_par},... % 1/dFdT_P dTdP_F];
                C_Solid_par{index_par}, C_Melt_par{index_par}, C_Cpx_par{index_par}, C_Gt_par{index_par}, errorPerplex{index_par}] = runPerplex(perplex_folder_2,index_par,input_complete);
            
        end
               
        S_System_vec = cell2mat(S_System_par);
        if sum(isnan(S_System_vec))>0
            disp('Error obtaining P-T at depth: S_System_diff NaN -exit')
            flag_error(3) = 1;
            C_solid = [2002];       C_fluid = [2002];           C_extracted = [2002]; C_ave = [2002]; C_solid_major = [2002];
            C_fluid_major = [2002]; C_extracted_major = [2002]; C_ave_major = [2002]; PTF = [2002];
            cd(working_folder);
            return
        end
        S_System_diff = S_System_0-S_System_vec;
        
        S_System_diff_positive = S_System_diff; S_System_diff_positive(S_System_diff<0) = 1e6;
        S_System_diff_negative = S_System_diff; S_System_diff_negative(S_System_diff>0) = -1e6;
        [~,ind_positive] = min(S_System_diff_positive);
        [~,ind_negative] = max(S_System_diff_negative);
        if abs(S_System_diff(ind_positive))<abs(S_System_diff(ind_negative))
            index_Na2O = [ind_positive ind_negative];
        else
            index_Na2O = [ind_negative ind_positive];
        end
       
        val_Na2O        = S_System_0-S_System_vec(index_Na2O);
        val_Na2O_min    = min(val_Na2O);
        val_Na2O_max    = max(val_Na2O);
        
        DeltaNa2O = (0-val_Na2O_min)/(val_Na2O_max-val_Na2O_min);
        
        X_Melt = X_Melt_par{index_Na2O(1)}*DeltaNa2O + X_Melt_par{index_Na2O(2)}*(1-DeltaNa2O);
        
        if all(S_System_0-S_System_vec>0)
            dT_vec_aux = max(dT_vec_aux) + dT_vec_aux_0;
        elseif  all(S_System_0-S_System_vec<0)
            dT_vec_aux = min(dT_vec_aux) - fliplr(dT_vec_aux_0);
        else
            exit_aux_T = 1;
            if X_Melt == 0
                exit_aux_P = 1;
            end
        end
    end
    
    
end

% correct alpha
ind = find(cellfun(@(x)x >8e-5, alpha_Solid_par));
if ~isempty(ind)
        disp('Error obtaining alpha: alpha_solid replaced with mean value')

ind_mean = cellfun(@(x) x <8e-5, alpha_Solid_par);
alpha_Solid_aux = [alpha_Solid_par{ind_mean}];
alpha_Solid_par{ind}=mean(alpha_Solid_aux);
end

[P,T_init,S_System, S_Solid, S_Melt, V_Solid, V_Melt, Cp_Solid, alpha_Solid, rho_Solid,...
    X_Melt, X_Ol, X_Cpx, X_Opx, X_Gt, X_Sp, X_Pl, Vol_Melt, Vol_Ol, Vol_Cpx, Vol_Opx, Vol_Gt, Vol_Sp, Vol_Pl,...
    Rho_Melt, Rho_Ol, Rho_Cpx, Rho_Opx, Rho_Gt, Rho_Sp, Rho_Pl,...
    H_Solid, H_Melt, H_Gt, H_Sp,... % 1/dFdT_P dTdP_F];
    C_Solid, C_Melt, C_Cpx, C_Gt] = ...
    interpolatePerplex(P_par, T_par, S_System_par, S_Solid_par, S_Melt_par, V_Solid_par, V_Melt_par, Cp_Solid_par, alpha_Solid_par, rho_Solid_par,...
    X_Melt_par, X_Ol_par, X_Cpx_par, X_Opx_par, X_Gt_par, X_Sp_par, X_Pl_par, Vol_Melt_par, Vol_Ol_par, Vol_Cpx_par, Vol_Opx_par, Vol_Gt_par, Vol_Sp_par, Vol_Pl_par,...
    Rho_Melt_par, Rho_Ol_par, Rho_Cpx_par, Rho_Opx_par, Rho_Gt_par, Rho_Sp_par, Rho_Pl_par,...
    H_Solid_par, H_Melt_par, H_Gt_par, H_Sp_par,... % 1/dFdT_P dTdP_F];
    C_Solid_par, C_Melt_par, C_Cpx_par, C_Gt_par,DeltaNa2O,index_Na2O);


%% save starting data to output array
% % % % %         P_init = dataMat(log_System,2);

dTdF_P = 250;                                   % [K]                   %% placeholder only, cancels out
dTdP_F = 130e-4;                                % [K/bar] = 130 K/GPa   %% placeholder only, cancels out

F = 0*[X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl];

% 1-10:  props1 = [P T S_System S_Solid S_Melt V_Solid V_Melt Cp_Solid alpha_Solid rho_Solid];
% 11-25: props2 = [F(1) X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl Vol_Melt Vol_Ol Vol_Cpx Vol_Opx Vol_Gt Vol_Sp Vol_Pl];
% 26-32: props3 = [Rho_Melt Rho_Ol Rho_Cpx Rho_Opx Rho_Gt Rho_Sp Rho_Pl];
% 33-38: props4 = [H_Solid H_Melt H_Gt H_Sp 1/dFdT_P dTdP_F];
% 39-56: props5 = [C_Solid C_Melt];
% 57-74: props6 = [C_Cpx C_Gt];
% output = [props1 props2 props3 props4 props5];

props1 = [P_init T_init S_System S_Solid S_Melt V_Solid V_Melt Cp_Solid alpha_Solid rho_Solid];
props2 = [F(1) X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl Vol_Melt Vol_Ol Vol_Cpx Vol_Opx Vol_Gt Vol_Sp Vol_Pl];
props3 = [Rho_Melt Rho_Ol Rho_Cpx Rho_Opx Rho_Gt Rho_Sp Rho_Pl];
props4 = [H_Solid H_Melt H_Gt H_Sp dTdF_P dTdP_F];
props5 = [C_Solid C_Melt];
if length(C_Cpx)==1; C_Cpx = 0*C_Solid; end
if length(C_Gt)==1;  C_Gt = 0*C_Solid;  end
props6 = [C_Cpx C_Gt];
if length(C_Melt)==1
    props5 = [C_Solid 0*C_Solid];
end
output = [props1 props2 props3 props4 props5 props6];

intervals = abs(round((P_init-P_top)/dP));
M_out = ones(intervals,length(output));     M_out(M_out==1) = NaN;
F_out = ones(intervals,length(F));          F_out(F_out==1) = NaN;
Error_out = ones(intervals,2);      Error_out(Error_out==1) = NaN;
M_out(1,:) = output;                                 % output array
F_out(1,:) = F;
Error_out(1,:) = [0 0];

%% Find solidus

F0 = F;
X0  = [X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl]/sum([X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl]);
n = 2;

dP_vec = dP:dP:parallel*dP;
T_last = T_init;
P_last = P_init;
hh = 0;

while X_Melt <= 0 && P_last > P_top
    %         T_init = T_init_0 + dT;
    hh = hh + 1;
    if alpha_Solid>8e-5 || alpha_Solid<5e-6
        disp('Error obtaining alpha: alpha_solid replace with standard value')
        alpha_Solid = alpha;
    end
    dT_vec = (alpha_Solid*T_last/(rho_Solid*Cp_Solid*1e-5)) * dP_vec;
    
    T_init_cell = num2cell(T_last + dT_vec);
    P_init_cell = num2cell(P_last + dP_vec);
    
    %         dT_iter = dT_iter + dT;
    
%    cd(perplex_folder);
%    filename_aux = pwd;
    %         tic
    parfor index_par = 1:parallel
%                             for index_par = 1:vertex
        
        input_complete = double([P_init_cell{index_par} T_init_cell{index_par} input_source]);
        
        [P_par{index_par}, T_par{index_par}, S_System_par{index_par}, S_Solid_par{index_par}, S_Melt_par{index_par}, V_Solid_par{index_par}, V_Melt_par{index_par}, Cp_System_par{index_par}, Cp_Solid_par{index_par}, alpha_Solid_par{index_par}, rho_Solid_par{index_par},...
            X_Melt_par{index_par}, X_Ol_par{index_par}, X_Cpx_par{index_par}, X_Opx_par{index_par}, X_Gt_par{index_par}, X_Sp_par{index_par}, X_Pl_par{index_par}, Vol_Melt_par{index_par}, Vol_Ol_par{index_par}, Vol_Cpx_par{index_par}, Vol_Opx_par{index_par}, Vol_Gt_par{index_par}, Vol_Sp_par{index_par}, Vol_Pl_par{index_par},...
            Rho_Melt_par{index_par}, Rho_Ol_par{index_par}, Rho_Cpx_par{index_par}, Rho_Opx_par{index_par}, Rho_Gt_par{index_par}, Rho_Sp_par{index_par}, Rho_Pl_par{index_par},...
            H_Solid_par{index_par}, H_Melt_par{index_par}, H_Gt_par{index_par}, H_Sp_par{index_par},... % 1/dFdT_P dTdP_F];
            C_Solid_par{index_par}, C_Melt_par{index_par}, C_Cpx_par{index_par}, C_Gt_par{index_par}, errorPerplex{index_par}] = runPerplex(perplex_folder_2,index_par,input_complete);
        
    end
    ind_melt = find(cell2mat(X_Melt_par),1);
    if isempty(ind_melt)
        ind_melt = parallel;
    end
    
    %         if ind_melt==1
    %             keyboard
    %         end
        
    P_vec = cell2mat(P_par);  P_vec = P_vec(1:ind_melt)';
    T_vec = cell2mat(T_par);  T_vec = T_vec(1:ind_melt)';
    
    
        if any(isnan(T_vec))
                disp('Error obtaining solidus: Temperature isnan - exit')
%                 T_init_cell
%                 P_init_cell
%                 input_source
                flag_error(4) = 1;
                C_solid = [2002];       C_fluid = [2002];           C_extracted = [2002]; C_ave = [2002]; C_solid_major = [2002];
                C_fluid_major = [2002]; C_extracted_major = [2002]; C_ave_major = [2002]; PTF = [2002];
                cd(working_folder);
                return
        end
        
    
    S_System_vec = cell2mat(S_System_par);  S_System_vec = S_System_vec(1:ind_melt)';
    S_Solid_vec = cell2mat(S_Solid_par);  S_Solid_vec = S_Solid_vec(1:ind_melt)';
    S_Melt_vec = cell2mat(S_Melt_par);  S_Melt_vec = S_Melt_vec(1:ind_melt)';
    V_Solid_vec = cell2mat(V_Solid_par);  V_Solid_vec = V_Solid_vec(1:ind_melt)';
    V_Melt_vec = cell2mat(V_Melt_par);  V_Melt_vec = V_Melt_vec(1:ind_melt)';
    Cp_Solid_vec = cell2mat(Cp_Solid_par);  Cp_Solid_vec = Cp_Solid_vec(1:ind_melt)';
    alpha_Solid_vec = cell2mat(alpha_Solid_par);  alpha_Solid_vec = alpha_Solid_vec(1:ind_melt)';
    rho_Solid_vec = cell2mat(rho_Solid_par);  rho_Solid_vec = rho_Solid_vec(1:ind_melt)';
    
    X_Melt_vec = cell2mat(X_Melt_par);  X_Melt_vec = X_Melt_vec(1:ind_melt)';
    X_Ol_vec = cell2mat(X_Ol_par);  X_Ol_vec = X_Ol_vec(1:ind_melt)';
    X_Cpx_vec = cell2mat(X_Cpx_par);  X_Cpx_vec = X_Cpx_vec(1:ind_melt)';
    X_Opx_vec = cell2mat(X_Opx_par);  X_Opx_vec = X_Opx_vec(1:ind_melt)';
    X_Gt_vec = cell2mat(X_Gt_par);  X_Gt_vec = X_Gt_vec(1:ind_melt)';
    X_Sp_vec = cell2mat(X_Sp_par);  X_Sp_vec = X_Sp_vec(1:ind_melt)';
    X_Pl_vec = cell2mat(X_Pl_par);  X_Pl_vec = X_Pl_vec(1:ind_melt)';
    Vol_Melt_vec = cell2mat(Vol_Melt_par);  Vol_Melt_vec = Vol_Melt_vec(1:ind_melt)';
    Vol_Ol_vec = cell2mat(Vol_Ol_par);  Vol_Ol_vec = Vol_Ol_vec(1:ind_melt)';
    Vol_Cpx_vec = cell2mat(Vol_Cpx_par);  Vol_Cpx_vec = Vol_Cpx_vec(1:ind_melt)';
    Vol_Opx_vec = cell2mat(Vol_Opx_par);  Vol_Opx_vec = Vol_Opx_vec(1:ind_melt)';
    Vol_Gt_vec = cell2mat(Vol_Gt_par);  Vol_Gt_vec = Vol_Gt_vec(1:ind_melt)';
    Vol_Sp_vec = cell2mat(Vol_Sp_par);  Vol_Sp_vec = Vol_Sp_vec(1:ind_melt)';
    Vol_Pl_vec = cell2mat(Vol_Pl_par);  Vol_Pl_vec = Vol_Pl_vec(1:ind_melt)';
    
    Rho_Melt_vec = cell2mat(Rho_Melt_par);  Rho_Melt_vec = Rho_Melt_vec(1:ind_melt)';
    Rho_Ol_vec = cell2mat(Rho_Ol_par);  Rho_Ol_vec = Rho_Ol_vec(1:ind_melt)';
    Rho_Cpx_vec = cell2mat(Rho_Cpx_par);  Rho_Cpx_vec = Rho_Cpx_vec(1:ind_melt)';
    Rho_Opx_vec = cell2mat(Rho_Opx_par);  Rho_Opx_vec = Rho_Opx_vec(1:ind_melt)';
    Rho_Gt_vec = cell2mat(Rho_Gt_par);  Rho_Gt_vec = Rho_Gt_vec(1:ind_melt)';
    Rho_Sp_vec = cell2mat(Rho_Sp_par);  Rho_Sp_vec = Rho_Sp_vec(1:ind_melt)';
    Rho_Pl_vec = cell2mat(Rho_Pl_par);  Rho_Pl_vec = Rho_Pl_vec(1:ind_melt)';
    
    H_Solid_vec = cell2mat(H_Solid_par);  H_Solid_vec = H_Solid_vec(1:ind_melt)';
    H_Melt_vec = cell2mat(H_Melt_par);  H_Melt_vec = H_Melt_vec(1:ind_melt)';
    H_Gt_vec = cell2mat(H_Gt_par);  H_Gt_vec = H_Gt_vec(1:ind_melt)';
    H_Sp_vec = cell2mat(H_Sp_par);  H_Sp_vec = H_Sp_vec(1:ind_melt)';
    
    C_Solid_vec = cell2mat(C_Solid_par');  C_Solid_vec = C_Solid_vec(1:ind_melt,:);
    C_Melt_vec = cell2mat(C_Melt_par');  C_Melt_vec = C_Melt_vec(1:ind_melt,:);
    C_Cpx_vec = cell2mat(C_Cpx_par');  C_Cpx_vec = C_Cpx_vec(1:ind_melt,:);
    C_Gt_vec = cell2mat(C_Gt_par');  C_Gt_vec = C_Gt_vec(1:ind_melt,:);
    
    dFdT_P = 0*X_Melt_vec;                             % [1/K] ; should be aorund 0.004
    dTdP_F = 0*X_Melt_vec;           % [K/bar] = 130 K/GPa ; should be around 0.0130  multicomponent clapeyron slope
    
    
    X_vec  = [X_Melt_vec X_Ol_vec X_Cpx_vec X_Opx_vec X_Gt_vec X_Sp_vec X_Pl_vec]./sum([X_Melt_vec X_Ol_vec X_Cpx_vec X_Opx_vec X_Gt_vec X_Sp_vec X_Pl_vec],2);
    F_vec = 0*X_vec;
    
    %% append data to output array
    % 1-10:  props1 = [P T S_System S_Solid S_Melt V_Solid V_Melt Cp_Solid alpha_Solid rho_Solid];
    % 11-25: props2 = [F(1) X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl Vol_Melt Vol_Ol Vol_Cpx Vol_Opx Vol_Gt Vol_Sp Vol_Pl];
    % 26-32: props3 = [Rho_Melt Rho_Ol Rho_Cpx Rho_Opx Rho_Gt Rho_Sp Rho_Pl];
    % 33-38: props4 = [H_Solid H_Melt H_Gt H_Sp 1/dFdT_P dTdP_F];
    % 39-56: props5 = [C_Solid C_Melt];
    % 57-74: props6 = [C_Cpx C_Gt];
    % output = [props1 props2 props3 props4 props5 props6];
    
    props1 = [P_vec T_vec S_System_vec S_Solid_vec S_Melt_vec V_Solid_vec V_Melt_vec Cp_Solid_vec alpha_Solid_vec rho_Solid_vec];
    props2 = [F_vec(:,1) X_Melt_vec X_Ol_vec X_Cpx_vec X_Opx_vec X_Gt_vec X_Sp_vec X_Pl_vec Vol_Melt_vec Vol_Ol_vec Vol_Cpx_vec Vol_Opx_vec Vol_Gt_vec Vol_Sp_vec Vol_Pl_vec];
    %         props2 = [X_Melt X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl Vol_Melt Vol_Ol Vol_Cpx Vol_Opx Vol_Gt Vol_Sp Vol_Pl];
    props3 = [Rho_Melt_vec Rho_Ol_vec Rho_Cpx_vec Rho_Opx_vec Rho_Gt_vec Rho_Sp_vec Rho_Pl_vec];
    props4 = [H_Solid_vec H_Melt_vec H_Gt_vec H_Sp_vec 1./dFdT_P dTdP_F];
    props5 = [C_Solid_vec C_Melt_vec];
    if length(C_Cpx_vec)==1; C_Cpx_vec = 0*C_Solid_vec; end
    if length(C_Gt_vec)==1;  C_Gt_vec = 0*C_Solid_vec;  end
    props6 = [C_Cpx_vec C_Gt_vec];
    if length(C_Melt_vec)==1
        props5 = [C_Solid_vec 0*C_Solid_vec];
    end
    output = [props1 props2 props3 props4 props5 props6];
    
    %         M_out(2:ind_melt+1,:) = output;
    %         F_out(2:ind_melt+1,:) = F_vec;
    M_out(2+parallel*(hh-1):parallel*(hh-1)+ind_melt+1,:) = output;
    F_out(2+parallel*(hh-1):parallel*(hh-1)+ind_melt+1,:) = F_vec;
    %         M_out(2:ind_melt+1,:) = [M_out; output];
    %         F_out(2:ind_melt+1,:) = [F_out; F_vec];
    
    % next iteraction
    
    P_last = P_vec(end);
    T_last = T_vec(end);
    alpha_Solid = alpha_Solid_vec(end);
    rho_Solid = rho_Solid_vec(end);
    Cp_Solid = Cp_Solid_vec(end);
    X_Melt = X_Melt_vec(end);
    n = parallel*(hh-1)+ind_melt+1;
    
end

%% run 2nd set of Perple_X calculation along infinitesimal adiabatic steps (calculating T at given P and S)
%%% iterate over entire depth range to recover isentropic PT path
%%% if melting occurs, change S, volume, bulk composition etc to reflect solid residue only

cd(perplex_folder);

[M_out, F_out, Error_out, flag_error] = runIsentrope_Connolly_parallel(M_out,F_out,Error_out,P_top,dP,dT,perplex_folder_2,working_folder,batch,n,f0,parallel,flag_error);

cd(working_folder);

if any(flag_error == 1)
    C_solid = [2002];       C_fluid = [2002];           C_extracted = [2002]; C_ave = [2002]; C_solid_major = [2002];
    C_fluid_major = [2002]; C_extracted_major = [2002]; C_ave_major = [2002]; PTF = [2002];
    return
end


