%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fe,mg,ca] = create_sample(al,bounds0)
% ... creates major-element compositions 
% from global correlations + natural variability.
% Only uniform distributions are implemented. 

bounds(1,1)=bounds0(1,3); bounds(2,1)=bounds0(2,3);   % Al
bounds(1,2)=bounds0(1,18); bounds(2,2)=bounds0(2,18); % Fe
bounds(1,3)=bounds0(1,19); bounds(2,3)=bounds0(2,19); % Mg
bounds(1,4)=bounds0(1,20); bounds(2,4)=bounds0(2,20); % Ca

%% compute the new values for Fe, Mg and Ca...    
mean_mg = 49.369 + (-4.106*al) + (0.343*al^2);  % mean MgO
mean_ca= -0.164 + (0.906*al );                  % mean CaO

% new random MgO from uniform distributions
xlow= mean_mg - 2.4 ; xup = mean_mg + 2.4;
mg = xlow + (xup-xlow)*rand(1,1);

% new random CaO from uniform distributions
xlow= mean_ca - 0.9; xup=mean_ca + 0.9;
ca = xlow + (xup-xlow)*rand(1,1);

   while (ca < 0.1 || ca < bounds(1,4) || ca > bounds(2,4))
     ca = xlow + (xup-xlow)*rand(1,1);
   end
   
% new Fe distribution as function of Ca
 mean_fe= 7.4527 + (0.5689*ca) + (-0.0863*ca^2); % mean FeO
 xlow= mean_fe - (-0.1*ca + 1.6);
 xup = mean_fe + (-0.216*ca + 1.25);
    if(ca < 1)
      xlow=6.1;
      xup =8.9;
    end
    
    if(xlow < bounds(1,2)); xlow=bounds(1,2);end
    if(xup > bounds(2,2)) ; xup = bounds(2,2);end

    fe = xlow + (xup-xlow)*rand(1,1);
end