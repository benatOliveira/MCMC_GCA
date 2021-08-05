function [M_out, F_out, Error_out, flag_error] = runIsentrope_Connolly_parallel(M_out,F_out,Error_out,P_top,dP,dT,perplex_folder_2,working_folder,batch,n,f0,vertex,flag_error)

if n==1
    M_out = [M_out;M_out;M_out];  M_out(2,1) = M_out(2,1)-dP;  M_out(1,1) = M_out(1,1)-2*dP;
    F_out = [F_out;F_out;F_out];
    Error_out = [Error_out;Error_out;Error_out];
    disp('Warning in Connolly: n=1 corrected to n=3')
    n = 3;
elseif n==2
    M_out = [M_out(1,:);M_out];  M_out(1,1) = M_out(1,1)-dP;
    F_out = [F_out(1,:);F_out];
    Error_out = [Error_out(1,:);Error_out];
    disp('Warning in Connolly: n=2 corrected to n=3')
    n = 3;
end
n = n - 2;

P_init = M_out(n,1);   P = P_init;
T_init = M_out(n,2);   T = T_init;
F = F_out(n,:);
X  = M_out(n,12:18)/sum(M_out(n,12:18));
S_Solid = M_out(n,4);
C_Solid = M_out(n,39:47);
S_System = M_out(n,4);
C_System = M_out(n,39:47);
input_source = M_out(n,39:47);
dFdT_P =  1/M_out(n,37);
dTdP_F =  M_out(n,38);
vertex = 11;

n = n + 1;

iter_S_outer=1;
cd(perplex_folder_2);

while P_init > P_top
    %% next Perple_X calculation
    T_init_0 = T;
    P_init_0 = P_init;
    if batch ~= 1
        S_System_0 = S_Solid;
        input_source = C_Solid;
    else
        S_System_0 = S_System;
        input_source = input_source;
    end
    
    index_Na2O = [1 1];
    exit_index = 0;
    iter_S = 1;
    while index_Na2O(1)==index_Na2O(2) && exit_index == 0
        if iter_S_outer==1
            %             dT_vec = -[1 2 3 4 5 6 7 8 9 10];
            if vertex == 4
                dT_vec = -[0 3 6 10];
            elseif vertex == 6
                dT_vec = -[0 2 4 6 8 10];
            else
                    dT_vec = -[0 1 2 3 4 5 6 7 8 9 10];
                
            end
        else
            if dT_closest<-5
                %                 dT_vec = dT_closest + [4 3 2 1 0 -1 -2 -3 -4 -5];
                if vertex == 4
                    dT_vec = dT_closest + [4 1 -2 -5];
                elseif vertex == 6
                    dT_vec = dT_closest + [4 2 0 -2 -4 -5];
                else
                    dT_vec = dT_closest + [5 4 3 2 1 0 -1 -2 -3 -4 -5];
                end
            else
                %                 dT_vec = -[1 2 3 4 5 6 7 8 9 10];
                if vertex == 4
                    dT_vec = -[0 3 6 10];
                elseif vertex == 6
                dT_vec = -[0 2 4 6 8 10];
                else
                    dT_vec = -[0 1 2 3 4 5 6 7 8 9 10];
                        
                end
            end
        end
        
        T_init_cell = num2cell(T_init_0 + dT_vec);
        P_init = P_init_0 + dP;
        
        parfor index_par = 1:vertex
            %                                     for index_par = 1:vertex
            
            %%% Gibbs free energy minimisation (with updated T, P, C) %%%
            
            input_complete = double([P_init T_init_cell{index_par} input_source]);
            
            [P_par{index_par}, T_par{index_par}, S_System_par{index_par}, S_Solid_par{index_par}, S_Melt_par{index_par}, V_Solid_par{index_par}, V_Melt_par{index_par}, Cp_System_par{index_par}, Cp_Solid_par{index_par}, alpha_Solid_par{index_par}, rho_Solid_par{index_par},...
                X_Melt_par{index_par}, X_Ol_par{index_par}, X_Cpx_par{index_par}, X_Opx_par{index_par}, X_Gt_par{index_par}, X_Sp_par{index_par}, X_Pl_par{index_par}, Vol_Melt_par{index_par}, Vol_Ol_par{index_par}, Vol_Cpx_par{index_par}, Vol_Opx_par{index_par}, Vol_Gt_par{index_par}, Vol_Sp_par{index_par}, Vol_Pl_par{index_par},...
                Rho_Melt_par{index_par}, Rho_Ol_par{index_par}, Rho_Cpx_par{index_par}, Rho_Opx_par{index_par}, Rho_Gt_par{index_par}, Rho_Sp_par{index_par}, Rho_Pl_par{index_par},...
                H_Solid_par{index_par}, H_Melt_par{index_par}, H_Gt_par{index_par}, H_Sp_par{index_par},... % 1/dFdT_P dTdP_F];
                C_Solid_par{index_par}, C_Melt_par{index_par}, C_Cpx_par{index_par}, C_Gt_par{index_par}, errorPerplex{index_par}] = runPerplex(perplex_folder_2,index_par,input_complete);
            
        end
        
        S_System_vec = cell2mat(S_System_par);
        S_System_diff = S_System_0-S_System_vec;
        
        S_System_diff_positive = S_System_diff; S_System_diff_positive(S_System_diff<0) = 1e6;
        S_System_diff_negative = S_System_diff; S_System_diff_negative(S_System_diff>0) = -1e6;
        [~,ind_positive] = min(S_System_diff_positive);
        [~,ind_negative] = max(S_System_diff_negative);
        if abs(S_System_diff(ind_positive))<abs(S_System_diff(ind_negative))
            index_Na2O = [ind_positive ind_negative];
            if S_System_diff_negative(index_Na2O(2)) == -1e6
                index_Na2O(2) = index_Na2O(1);
            end
            
        else
            index_Na2O = [ind_negative ind_positive];
            if S_System_diff_positive(index_Na2O(2)) == 1e6
                index_Na2O(2) = index_Na2O(1);
            end
        end
        
        dT_closest = dT_vec(index_Na2O(1));
        if ( index_Na2O(1)==1 && index_Na2O(2)==1 ) && dT_closest == 0
            exit_index = 1;
        end
        iter_S=iter_S+1;
        if iter_S > 10
            disp('Error in computing the isentrope: Maximum iterations reached in runIsentrope - exit')
            %                 T_init_0
            %                 P_init_0
            %                 input_source
            flag_error(5) = 1;
            return
        end
        
    end
    
    val_Na2O        = S_System_0-S_System_vec(index_Na2O);
    val_Na2O_min    = min(val_Na2O);
    val_Na2O_max    = max(val_Na2O);
    
    DeltaNa2O = (0-val_Na2O_min)/(val_Na2O_max-val_Na2O_min);
    if exit_index == 1
        DeltaNa2O = 0.5;
    end
    
    [P,T,S_System, S_Solid, S_Melt, V_Solid, V_Melt, Cp_Solid, alpha_Solid, rho_Solid,...
        X_Melt, X_Ol, X_Cpx, X_Opx, X_Gt, X_Sp, X_Pl, Vol_Melt, Vol_Ol, Vol_Cpx, Vol_Opx, Vol_Gt, Vol_Sp, Vol_Pl,...
        Rho_Melt, Rho_Ol, Rho_Cpx, Rho_Opx, Rho_Gt, Rho_Sp, Rho_Pl,...
        H_Solid, H_Melt, H_Gt, H_Sp,... % 1/dFdT_P dTdP_F];
        C_Solid, C_Melt, C_Cpx, C_Gt] = ...
        interpolatePerplex(P_par, T_par, S_System_par, S_Solid_par, S_Melt_par, V_Solid_par, V_Melt_par, Cp_Solid_par, alpha_Solid_par, rho_Solid_par,...
        X_Melt_par, X_Ol_par, X_Cpx_par, X_Opx_par, X_Gt_par, X_Sp_par, X_Pl_par, Vol_Melt_par, Vol_Ol_par, Vol_Cpx_par, Vol_Opx_par, Vol_Gt_par, Vol_Sp_par, Vol_Pl_par,...
        Rho_Melt_par, Rho_Ol_par, Rho_Cpx_par, Rho_Opx_par, Rho_Gt_par, Rho_Sp_par, Rho_Pl_par,...
        H_Solid_par, H_Melt_par, H_Gt_par, H_Sp_par,... % 1/dFdT_P dTdP_F];
        C_Solid_par, C_Melt_par, C_Cpx_par, C_Gt_par,DeltaNa2O,index_Na2O);
    
%     if X_Melt == 0
%         keyboard
%     end
    if isnan(T)
        disp('Error in computing the isentrope: Temperature isnan in runIsentrope - exit')
%         T_init_0
%         P_init_0
%         input_source
        flag_error(6) = 1;
        return
    end
    if X_Melt<0
        disp('Error in computing the isentrope: melt negative in runIsentrope - exit')
%         T_init_0
%         P_init_0
%         input_source
        flag_error(7) = 1;
        return
    end
    
    if isnan(Cp_Solid)
        disp('Error in computing the isentrope: Cp_system NaN - exit')
        flag_error(8) = 1;
        return
    end
    
    dP = dP;
    F0 = F;
    X0 = [0 X(2:end)]/sum(X(2:end));
    if batch == 1
        X0 = X(1:end)/sum(X(1:end));
    end
    X  = [X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl]/sum([X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl]);
    F = F0 + (1-F0(1)).*(X-X0);
    if batch == 1
        F = X;
        f0 = X(1);
    end
    
    %% append data to output array
    
    % 1-10:  props1 = [P T S_System S_Solid S_Melt V_Solid V_Melt Cp_Solid alpha_Solid rho_Solid];
    % 11-25: props2 = [F(1) X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl Vol_Melt Vol_Ol Vol_Cpx Vol_Opx Vol_Gt Vol_Sp Vol_Pl];
    % 26-32: props3 = [Rho_Melt Rho_Ol Rho_Cpx Rho_Opx Rho_Gt Rho_Sp Rho_Pl];
    % 33-38: props4 = [H_Solid H_Melt H_Gt H_Sp Vp_Solid Vs_Solid];
    % 39-56: props5 = [C_Solid C_Melt];
    % 57-74: props6 = [C_Cpx C_Gt];
    % output = [props1 props2 props3 props4 props5 props6];
    
    % props1 = [P T S_System S_Solid S_Melt V_Solid V_Melt Cp_Solid alpha_Solid rho_Solid];
    props1 = [P T S_System S_Solid S_Melt V_Solid V_Melt Cp_Solid alpha_Solid rho_Solid];
    props2 = [F(1) X_Melt X_Ol X_Cpx X_Opx X_Gt X_Sp X_Pl Vol_Melt Vol_Ol Vol_Cpx Vol_Opx Vol_Gt Vol_Sp Vol_Pl];
    props3 = [Rho_Melt Rho_Ol Rho_Cpx Rho_Opx Rho_Gt Rho_Sp Rho_Pl];
    props4 = [H_Solid H_Melt H_Gt H_Sp 1/dFdT_P dTdP_F];
    props5 = [C_Solid C_Melt];
    if length(C_Cpx)==1; C_Cpx = 0*C_Solid; end
    if length(C_Gt)==1;  C_Gt = 0*C_Solid;  end
    props6 = [C_Cpx C_Gt];
    if length(C_Melt)==1
        props5 = [C_Solid 0*C_Solid];
    end
    output = [props1 props2 props3 props4 props5 props6];
    
    M_out(n,:) = output;
    F_out(n,:) = F;
    % % % % Error_out(n,:) = [iter_S-1 error_accum];
    n = n +1;
    iter_S_outer=iter_S_outer+1;
    
    
    
end
