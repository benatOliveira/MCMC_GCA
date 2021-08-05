function [K_j] = changePartitionCoefficient(K_j,Melt_comp,Cpx_comp,Grt_comp,P,T,TP_weight,grid_z)

% Normalize input compositions
Melt_comp  = Melt_comp./repmat(sum(Melt_comp,2),1,9)*100;
Cpx_comp   = Cpx_comp./repmat(sum(Cpx_comp,2),1,9)*100;
Grt_comp   = Grt_comp./repmat(sum(Grt_comp,2),1,9)*100;

% Check if there is any compositional GAP
if sum(~isnan(Cpx_comp(:,1)))>1
Cpx_comp(isnan(Cpx_comp(:,1)) & TP_weight(:,3)>0,:) = interp1(grid_z(~isnan(Cpx_comp(:,1))),Cpx_comp(~isnan(Cpx_comp(:,1)),:),grid_z(isnan(Cpx_comp(:,1)) & TP_weight(:,3)>0),'linear','extrap');
end
if sum(~isnan(Grt_comp(:,1)))>1
Grt_comp(isnan(Grt_comp(:,1)) & TP_weight(:,5)>0,:) = interp1(grid_z(~isnan(Grt_comp(:,1))),Grt_comp(~isnan(Grt_comp(:,1)),:),grid_z(isnan(Grt_comp(:,1)) & TP_weight(:,5)>0),'linear','extrap');
end

%% Cpx - Wood Blundy 1997

factor_Cpx = 12./((4*Cpx_comp(:,1)/60.084)+(4*Cpx_comp(:,9)/79.879)+(6*Cpx_comp(:,2)/101.961)+(2*Cpx_comp(:,3)/71.846)+(2*Cpx_comp(:,4)/70.937)+(2*Cpx_comp(:,5)/40.304)+(2*Cpx_comp(:,6)/56.077)+(2*Cpx_comp(:,7)/61.979)+(6*Cpx_comp(:,8)/151.99)+(2*0/94.196));

Si_Cpx = factor_Cpx.*Cpx_comp(:,1)/60.084;
Ti_Cpx = factor_Cpx.*Cpx_comp(:,9)/79.879;
Al_Cpx = 2*factor_Cpx.*Cpx_comp(:,2)/101.961;
Aliv_Cpx = 2-Si_Cpx;
Alvi_Cpx = Al_Cpx-Aliv_Cpx;
Fe_Cpx = factor_Cpx.*Cpx_comp(:,3)/71.846;
Mn_Cpx = factor_Cpx.*Cpx_comp(:,4)/70.937;
Mg_Cpx = factor_Cpx.*Cpx_comp(:,5)/40.304;
Ca_Cpx = factor_Cpx.*Cpx_comp(:,6)/56.077;
Na_Cpx = 2*factor_Cpx.*Cpx_comp(:,7)/61.979;
Cr_Cpx = 2*factor_Cpx.*Cpx_comp(:,8)/151.99;
K_Cpx  = 2*factor_Cpx*0/94.196;

Mg_Melt = Melt_comp(:,5)/40.304;
Fe_Melt = Melt_comp(:,3)/71.846;

Mg_number_Cpx = Mg_Cpx./(Mg_Cpx+Fe_Cpx);
Mg_number_Melt = Mg_Melt./(Mg_Melt+Fe_Melt);
XMg = Mg_number_Cpx.*(1-Alvi_Cpx-Ti_Cpx-Cr_Cpx);

%  1     2     3     4     5     6     7     8     9     10    11    12    13       
%  La    Ce    Pr    Nd    Sm    Eu    Gd    Tb    Dy    Ho    Er    Yb    Lu
r0 = [1.16	1.143	1.126	1.109	1.079	1.066	1.053	1.04	1.027	1.015	1.004	0.985	0.977];


E_Cpx = repmat(318.6+6.9*P/1e4-0.036*T,1,length(r0));
r_Cpx = repmat(0.974+0.067*Ca_Cpx-0.051*Alvi_Cpx,1,length(r0));
D0_Cpx = repmat((Mg_number_Melt./XMg).*exp((88750-65.644*(T)+0.705*(P/1e4)*10000-0.077*(P/1e4).^2*10000)./(8.314*(T))),1,length(r0));
D_Cpx = D0_Cpx.*exp((-910.17.*E_Cpx./(T)).*((1/3).*(repmat(r0,size(r_Cpx,1),1)-r_Cpx).^3+0.5.*r_Cpx.*(repmat(r0,size(r_Cpx,1),1)-r_Cpx).^2)); D_Cpx(1,:) = D_Cpx(2,:);
D_Cpx(isinf(D_Cpx)) = NaN;

%% Grt - Van Westrenen et. al. 2001

total_Grt = sum(Grt_comp,2);
oxygen_Grt = ((Grt_comp(:,5)/40.31)+(Grt_comp(:,6)/56.08)+2*(Grt_comp(:,1)/60.1)+3*(Grt_comp(:,2)/102)+(Grt_comp(:,3)/71.85)+3*(0/159.7)+(Grt_comp(:,4)/70.94)+2*(Grt_comp(:,9)/79.9)+(0/94.2)+(Grt_comp(:,7)/61.98)+5*(0/141.9)+3*(Grt_comp(:,8)/151.99))./total_Grt;
DMg_Grt = Grt_comp(:,5)./Melt_comp(:,5);
Mg_Grt = (Grt_comp(:,5)/40.31)./total_Grt.*(12./oxygen_Grt);
Ca_Grt = (Grt_comp(:,6)/56.08)./total_Grt.*(12./oxygen_Grt);	
Si_Grt = (Grt_comp(:,1)/60.1)./total_Grt.*(12./oxygen_Grt);		
Al_Grt = 2*(Grt_comp(:,2)/102)./total_Grt.*(12./oxygen_Grt);	
Fe_Grt = (Grt_comp(:,3)/71.85)./total_Grt.*(12./oxygen_Grt);		
Fe3_Grt	= 0*Fe_Grt;	
Mn_Grt = (Grt_comp(:,4)/70.94)./total_Grt.*(12./oxygen_Grt);	
Ti_Grt = (Grt_comp(:,9)/79.9)./total_Grt.*(12./oxygen_Grt);	
K_Grt	= 0*Ti_Grt;
Na_Grt = 2*(Grt_comp(:,7)/61.98)./total_Grt.*(12./oxygen_Grt);	
Pp_Grt	= 0*Na_Grt;	
Cr_Grt = 2*(Grt_comp(:,8)/151.99)./total_Grt.*(12./oxygen_Grt);	

XPyrope         = Mg_Grt./(Mg_Grt+Ca_Grt+Fe_Grt+Mn_Grt);
XAlmandine      = Fe_Grt./(Mg_Grt+Ca_Grt+Fe_Grt+Mn_Grt);
XSpessartine	= Mn_Grt./(Mg_Grt+Ca_Grt+Fe_Grt+Mn_Grt);
XAndradite      = Fe3_Grt./(Fe3_Grt+Cr_Grt+Al_Grt);
XUvarovite      = Cr_Grt./(Cr_Grt+Fe3_Grt+Al_Grt);
XGrossular      = Ca_Grt./(Mg_Grt+Ca_Grt+Fe_Grt+Mn_Grt)-XUvarovite-XAndradite;

r_Grt = repmat(0.9302*XPyrope+0.993*XGrossular+0.916*XAlmandine+0.946*XSpessartine+1.05*(XAndradite+XUvarovite)-0.005*(P/1e4-3),1,length(r0));
E_Grt = repmat(3.5*10^12*(1.38+r_Grt(:,1)).^(-26.7),1,length(r0));
gMg_Grt = exp(19000*(XGrossular).^2./(8.314*T));
RTlnDgD0 = (418+10.4*P/1e4-0.226*T)*1000;
D0_Grt = repmat(exp(RTlnDgD0./(8.314*T))./((DMg_Grt.^2).*(gMg_Grt.^2)),1,length(r0));
D_Grt = D0_Grt.*exp(-910.17.*E_Grt.*(0.5*r_Grt.*(repmat(r0,size(r_Grt,1),1)-r_Grt).^2+(1/3)*(repmat(r0,size(r_Grt,1),1)-r_Grt).^3)./(T)); D_Grt(1,:) = D_Grt(2,:);
D_Grt(isinf(D_Grt)) = NaN;


%% Rearrange data

for index_REE = 1:length(r0)
    REE_K =  cell2mat(K_j(:,index_REE));
    REE_K(~isnan(D_Cpx(:,index_REE)),2) = D_Cpx(~isnan(D_Cpx(:,index_REE)),index_REE);
    REE_K(~isnan(D_Grt(:,index_REE)),4) = D_Grt(~isnan(D_Grt(:,index_REE)),index_REE);
    K_j(:,index_REE) = mat2cell(REE_K,ones(size(REE_K,1),1),size(REE_K,2));
end
    


