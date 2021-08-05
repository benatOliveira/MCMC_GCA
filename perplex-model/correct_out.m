function [M_out_smooth,F_out_smooth,grid_z_smooth,PTF_smooth] = correct_out(M_out,F_out,ztop,batch,parameters)

%% F always increases? Correct it

index_noMelt = M_out(:,12)==0 & M_out(:,11)>0;
aux_vec = [];
correct_vec = [];
rho_m   = parameters.model.rho_m;                           % mantle reference density [kg/m3]
g   = parameters.model.g;                           % mantle reference density [kg/m3]
for index_nz = 2:length(index_noMelt)
    if index_noMelt(index_nz)
        aux_vec = [aux_vec index_nz];
        correct_vec = [correct_vec index_nz];
        if index_nz~=length(index_noMelt)
            if ~index_noMelt(index_nz+1)
                % Correct melt wt
                M_out(aux_vec,12) = M_out(index_nz+1,12)/(length(aux_vec)+1);       % distribute melt
                M_out(index_nz+1,12) = M_out(index_nz+1,12)/(length(aux_vec)+1);    % correct melt content in node index_nz
                M_out(index_nz+1,13:18) = M_out(index_nz+1,13:18) - repmat(M_out(index_nz+1,12),1,6).*M_out(index_nz+1,13:18)./repmat(sum(M_out(index_nz+1,13:18),2),1,6);      % correct TP in node index_nz
                M_out(index_nz+1,13:18) = repmat(100-M_out(index_nz+1,12),1,6).*M_out(index_nz+1,13:18)./repmat(sum(M_out(index_nz+1,13:18),2),1,6);      % normalize TP keeping melt cte in node index_nz
                M_out_solid0 = M_out(aux_vec(1)-1,13:18)./repmat(sum(M_out(aux_vec(1)-1,13:18),2),1,6)*100;    % save departing TP composition (normalized to solid)
                M_out_solid1 = M_out(index_nz+1,13:18)./repmat(sum(M_out(index_nz+1,13:18),2),1,6)*100;        % save arriving TP composition  (normalized to solid)
                M_out(aux_vec,13:18) = repmat(M_out_solid0,length(aux_vec),1) +  repmat([1:length(aux_vec)]'/(length(aux_vec)+1),1,6).*repmat(M_out_solid1-M_out_solid0,length(aux_vec),1);    % correct with linear interpolation
                M_out(aux_vec,13:18) = M_out(aux_vec,13:18) - repmat(M_out(aux_vec,12),1,6).*M_out(aux_vec,13:18)./repmat(sum(M_out(aux_vec,13:18),2),1,6);     % add melt into total composition
                M_out(aux_vec,13:18) = repmat(100-M_out(aux_vec,12),1,6).*M_out(aux_vec,13:18)./repmat(sum(M_out(aux_vec,13:18),2),1,6);      % normalize TP keeping melt cte in interpolated nodes
                % Correct melt vol   - Same steps as before
                M_out(aux_vec,19) = M_out(index_nz+1,19)/(length(aux_vec)+1);
                M_out(index_nz+1,19) = M_out(index_nz+1,19)/(length(aux_vec)+1);
                M_out(index_nz+1,20:25) = M_out(index_nz+1,20:25) - repmat(M_out(index_nz+1,19),1,6).*M_out(index_nz+1,20:25)./repmat(sum(M_out(index_nz+1,20:25),2),1,6);  % normalize
                M_out(index_nz+1,20:25) = repmat(100-M_out(index_nz+1,19),1,6).*M_out(index_nz+1,20:25)./repmat(sum(M_out(index_nz+1,20:25),2),1,6);
                M_out_solid0 = M_out(aux_vec(1)-1,20:25)./repmat(sum(M_out(aux_vec(1)-1,20:25),2),1,6)*100;
                M_out_solid1 = M_out(index_nz+1,20:25)./repmat(sum(M_out(index_nz+1,20:25),2),1,6)*100;
                M_out(aux_vec,20:25) = repmat(M_out_solid0,length(aux_vec),1) +  repmat([1:length(aux_vec)]'/(length(aux_vec)+1),1,6).*repmat(M_out_solid1-M_out_solid0,length(aux_vec),1);
                M_out(aux_vec,20:25) = M_out(aux_vec,20:25) - repmat(M_out(aux_vec,19),1,6).*M_out(aux_vec,20:25)./repmat(sum(M_out(aux_vec,20:25),2),1,6);
                M_out(aux_vec,20:25) = repmat(100-M_out(aux_vec,19),1,6).*M_out(aux_vec,20:25)./repmat(sum(M_out(aux_vec,20:25),2),1,6);
                % Correct F
                M_out(aux_vec,11) = M_out(aux_vec,11)+ [1:length(aux_vec)]'/(length(aux_vec)+1)* (M_out(index_nz+1,11)-M_out(index_nz,11));
            end
        else   % correct if last node in the melting column has no melt (by projecting backwards)
            % Correct melt wt
            M_out(aux_vec,12) = M_out(aux_vec(1)-1,12);
            M_out(aux_vec,13:18) = M_out(aux_vec,13:18) - repmat(M_out(aux_vec,12),1,6).*M_out(aux_vec,13:18)./repmat(sum(M_out(aux_vec,13:18),2),1,6);
            % Correct melt vol
            M_out(aux_vec,19) = M_out(aux_vec(1)-1,19);
            M_out(aux_vec,20:25) = M_out(aux_vec,20:25) - repmat(M_out(aux_vec,19),1,6).*M_out(aux_vec,20:25)./repmat(sum(M_out(aux_vec,20:25),2),1,6);
            % Correct F
            M_out(aux_vec,11) = M_out(aux_vec,11)+ [1:length(aux_vec)]'*(M_out(aux_vec(1)-1,11)-M_out(aux_vec(1)-2,11));
        end
    else
        aux_vec = [];
    end
end

%% Smooth data - only if there is melt

if sum(M_out(:,12))>0
    
    smooth_factor = 3;
    M_out_smooth = M_out;
    M_out_smooth(M_out(:,11)>0,11) = smooth(M_out(M_out(:,11)>0,11),smooth_factor);  % F
    M_out_smooth(M_out(:,12)>0,12) = smooth(M_out(M_out(:,12)>0,12),smooth_factor);  % melt  - wt
    M_out_smooth(M_out(:,13)>0,13) = smooth(M_out(M_out(:,13)>0,13),smooth_factor);  % ol    - wt
    M_out_smooth(M_out(:,14)>0,14) = smooth(M_out(M_out(:,14)>0,14),smooth_factor);  % cpx   - wt
    M_out_smooth(M_out(:,15)>0,15) = smooth(M_out(M_out(:,15)>0,15),smooth_factor);  % opx   - wt
    M_out_smooth(M_out(:,16)>0,16) = smooth(M_out(M_out(:,16)>0,16),smooth_factor);  % grt   - wt
    M_out_smooth(M_out(:,17)>0,17) = smooth(M_out(M_out(:,17)>0,17),smooth_factor);  % sp    - wt
    M_out_smooth(M_out(:,18)>0,18) = smooth(M_out(M_out(:,18)>0,18),smooth_factor);  % pl    - wt
    M_out_smooth(:,12:18) = M_out_smooth(:,12:18)./repmat(sum(M_out_smooth(:,12:18),2),1,7)*100;
    
    M_out_smooth(M_out(:,19)>0,19) = smooth(M_out(M_out(:,19)>0,19),smooth_factor);  % melt  - vol
    M_out_smooth(M_out(:,20)>0,20) = smooth(M_out(M_out(:,20)>0,20),smooth_factor);  % ol    - vol
    M_out_smooth(M_out(:,21)>0,21) = smooth(M_out(M_out(:,21)>0,21),smooth_factor);  % cpx   - vol
    M_out_smooth(M_out(:,22)>0,22) = smooth(M_out(M_out(:,22)>0,22),smooth_factor);  % opx   - vol
    M_out_smooth(M_out(:,23)>0,23) = smooth(M_out(M_out(:,23)>0,23),smooth_factor);  % grt   - vol
    M_out_smooth(M_out(:,24)>0,24) = smooth(M_out(M_out(:,24)>0,24),smooth_factor);  % sp    - vol
    M_out_smooth(M_out(:,25)>0,25) = smooth(M_out(M_out(:,25)>0,25),smooth_factor);  % pl    - vol
    M_out_smooth(:,19:25) = M_out_smooth(:,19:25)./repmat(sum(M_out_smooth(:,19:25),2),1,7);
    
    
    %% F out of smoothed data
    F = 0*F_out(1,:);
    F_out_smooth = 0*F_out;
    ind = find(M_out_smooth(:,12)>0); ind_a = ind(1);
    for index = ind_a:size(M_out_smooth,1)
        F0 = F;
        X0 = [0 M_out_smooth(index-1,13:18)]/sum(M_out_smooth(index-1,13:18));
        if batch == 1  % not tested
            X0 = X(1:end)/sum(X(1:end));
        end
        X  = M_out_smooth(index,12:18)/sum(M_out_smooth(index,12:18));
%         if X(1)<1e-5
%             F = 0*F0;
%         else
            F = F0 + (1-F0(1)).*(X-X0);  % this is how F is computed!
%         end
        F_out_smooth(index,:) = F;
    end
    
    
else
    F_out_smooth = F_out;
    M_out_smooth = M_out;
end

%% save melt path
PTF        = [M_out(:,1),M_out(:,2)-273.15,M_out(:,11)];
PTF_smooth = [M_out_smooth(:,1),M_out_smooth(:,2)-273.15,M_out_smooth(:,11)];


%% pressure to grid
P_grid = [M_out(:,1); 0]*1e5;
rho_s  = M_out(end:-1:1,10);
rho_s(1) = rho_m;

DeltaP = P_grid(end-1:-1:1)-P_grid(end:-1:2);
dh = DeltaP./(g*rho_s);

grid_z_smooth = flipud(cumsum(dh));


% %% Refinement option
% 
% [M_out_smooth,F_out_smooth,grid_z_smooth,PTF_smooth] = refinement1D(M_out_smooth,F_out_smooth,ztop,parameters);

