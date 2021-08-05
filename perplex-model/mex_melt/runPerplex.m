function [P, T, S_System, S_Solid, S_Melt, V_Solid,  V_Melt,  Cp_System,  Cp_Solid,  alpha_Solid,  rho_Solid,  ...
    X_Melt,  X_Ol,  X_Cpx,  X_Opx,  X_Gt,  X_Sp,  X_Pl,  Vol_Melt,  Vol_Ol,  Vol_Cpx,  Vol_Opx,  Vol_Gt,  Vol_Sp,  Vol_Pl,  ...
    Rho_Melt,  Rho_Ol,  Rho_Cpx,  Rho_Opx,  Rho_Gt,  Rho_Sp,  Rho_Pl,  ...
    H_Solid,  H_Melt,  H_Gt,  H_Sp,  ... % 1/dFdT_P dTdP_F];
    C_Solid,  C_Melt,  C_Cpx,  C_Gt, e] = ...
    runPerplex(perplex_folder_2,vertex,input_complete)

% %             filename = strcat(filename_aux,'\tpcdata_',num2str(vertex),'.dat');
% %             save(filename, 'input_complete', '-ASCII');
% %             cmd = strcat('"project_',num2str(vertex),'.exe" tpcdata_',num2str(vertex),'.dat');
% %             system(cmd);
% %             read_calcdata;
            
            v = input_complete(1:2);
            cblk = input_complete(3:end);
    
% load('test_data_doubleOPX.mat')
% v = aa(1:2);
% cblk = aa(3:end);
%             cd(perplex_folder_2)
            [a,b,c,d,e]=meemum_fun_melt_gateway(v,cblk);
            read_mexFile;
            
    if length(C_Solid)==1; C_Solid = 0*cblk; end
    if length(C_Cpx)==1; C_Cpx = 0*C_Solid; end
    if length(C_Gt)==1;  C_Gt = 0*C_Solid;  end
    
    if length(C_Melt)==1
        C_Melt = 0*C_Solid;
    end
    
  