% Creates Figures for the RGR example

%% PRELIMINARY
% burn
burn = options.burn;

%REE
% mean of melt trace element composition
meltREE = melt_erupt(burn:end,1:13);
mean_meltREE = sum(meltREE)/size(meltREE,1);
% variance
var_meltREE = sum((meltREE-mean_meltREE).^2)/size(meltREE,1);
std_meltREE=[mu(1)*(13/100) mu(2)*(12/100) mu(3)*(11/100) mu(4)*(9/100) mu(5)*(7.5/100) mu(6)*(7.5/100) mu(7)*(6.3/100) mu(8)*(7/100) mu(9)*(6.26/100) mu(10)*(7./100) mu(11)*(8./100) mu(12)*(7.6/100) mu(13)*(14./100)];

%majors
% mean of melt trace element composition
meltMAJORS = melt_erupt(burn:end,14:22);
mean_meltMAJORS = sum(meltMAJORS)/size(meltMAJORS,1);
% variance
var_meltMAJORS = sum((meltMAJORS-mean_meltMAJORS).^2)/size(meltMAJORS,1);
std_meltMAJORS=[mu(14)*(3/100) mu(15)*(5.5/100) mu(16)*(8.6/100) mu(17)*(7.5/100) mu(18)*(15./100) mu(19)*(7.5/100) mu(20)*(11./100) mu(21)*(6/100) mu(22)*(9./100)];




% number of simulations after burning
nsimu = size(meltREE,1);%size(chainnew,1);


%% FIGURE 1
auxREEplot = 0.25;
% Define bounds of plots
minREE = [0.01  0.01  0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01];   % requires adjustment; problem dependent
maxREE = [190 170 170 120 90   90   80   80   80   80   80   50   50];           % requires adjustment; problem dependent

figure(1);
clf(1);
for indexREE = 1:13
scattercloud_1D(indexREE*ones(1,nsimu),meltREE(:,indexREE),100,2,'k.',flipud(hot(256)),indexREE-auxREEplot ,indexREE+auxREEplot,minREE(indexREE),maxREE(indexREE)); %axis square;xlim([minT maxT]);hold on;
hold on
end

xlim([0.5 13.5])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13],'XTickLabel',...
    {'La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Yb','Lu'});
errorbar(mu(1:13),std_meltREE,'o','color','b','CapSize',10);
scatter([1:13],mu(1:13),50,'o','filled','b','MarkerEdgeColor','k')
ylabel('Aggregated melts [ppm, CI normalized]')

%% FIGURE 2
std_meltOX = std_meltMAJORS;
meltOX = meltMAJORS;
auxOXplot = 0.25;
% Define bounds of plots
minOX = [35 10 5  5   5   0.01 0.01 ];   % requires adjustment; problem dependent
maxOX = [60 30 30 25  25  20   20   ];   % requires adjustment; problem dependent

% errorbar
mean_bar = [mu(14) mu(15) mu(16) mu(18) mu(19) mu(20) mu(21)];
std_bar  = [std_meltOX(1) std_meltOX(2) std_meltOX(3) std_meltOX(5) std_meltOX(6) std_meltOX(7) std_meltOX(8)];
figure(2);
clf(2);
for indexOX = 1:7
scattercloud_1D(indexOX*ones(1,nsimu),100*meltOX(:,indexOX),100,2,'k.',flipud(hot(256)),indexOX-auxOXplot ,indexOX+auxOXplot,minOX(indexOX),maxOX(indexOX)); %axis square;xlim([minT maxT]);hold on;
hold on
end
xlim([0.5 9.5])
ylim([0 55])
set(gca,'XTick',[1 2 3 4 5 6 7],'XTickLabel',...
    {'SiO_2','Al_2O_3','FeO','MgO','CaO','Na_2O','Cr_2O_3'});
errorbar(100*mean_bar,100*std_bar,'o','color','b','CapSize',10);
scatter([1:7],100*mean_bar,50,'o','filled','b','MarkerEdgeColor','k')
ylabel('Aggregated melts [wt%]')


%% FIGURE 3

% which data do you want to plot?
             % Tp,ztop,Al,Na,  La,Ce,Pr,Nd,Sm,Eu,Gd,Tb,Dy,Ho,Er,Yb,Lu,  FeO,MgO,CaO,C2rO3
aux_index = [  1  1    1  1    1  0  0  0  1  0  0  0  0  0  0  0  1    0   0   0   0];

chainnew1 = chain(burn:end,:); % set the burn-in (real dimension of the chain)
chainnew1(:,5:17) = chainnew1(:,5:17)/100;

% Define bounds of plots
minChain = 0.9*min(chainnew1) ;  
minChain(1:4) = [1330 50 3.0   0.15];
maxChain = 1.1*max(chainnew1) ;  
maxChain(1:4) = [1345 65 3.4 0.25];
label_aux = {'Tp','ztop','Al','Na','La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Yb','Lu','FeO','MgO','CaO','C2rO3'};
figure(3);
clf(3);
for index_Y = 1:length(aux_index)
    index_SUM_Y = sum(aux_index(1:index_Y));
    for index_X = 1:length(aux_index)
        index_SUM_X = sum(aux_index(1:index_X));
        if index_Y>index_X && aux_index(index_X)==1 && aux_index(index_Y)==1 
            subplot(sum(aux_index),sum(aux_index),(index_SUM_Y-1)*(sum(aux_index))+index_SUM_X)
            scattercloud(chainnew1(:,index_X),chainnew1(:,index_Y),60,2,'k.',flipud(hot(256)),minChain(index_X),maxChain(index_X),minChain(index_Y),maxChain(index_Y)); axis square;hold on;
%             xlabel('Tp [C]'); ylabel('ztop [km]')
            if index_SUM_X == 1
                ylabel(label_aux(index_Y))
            end
            if index_SUM_Y == sum(aux_index)
                xlabel(label_aux(index_X))
            end
        end
    end
end



%% FIGURE 4
auxREEplot = 0.25;
minREE = [0.5  0.5  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
maxREE = 2*[2 2 2 2 2   2   2   2   2   2   2   2   2];

meltREE_source = chain(burn:end,5:17)/100; 

%                  La      Ce      Pr      Nd      Sm      Eu      Gd      Tb    Dy      Ho      Er      Yb      Lu
REE_factor =      [237     613     92.8    457     148     56.3    199     36.1  246     54.6    160     161     24.6]/1000; % CI from McDonough & Sun (1995)
meltREE_source = meltREE_source./repmat(REE_factor,size(meltREE_source,1),1);

figure(4);
clf(4);
for indexREE = 1:13
scattercloud_1D(indexREE*ones(1,nsimu),meltREE_source(:,indexREE),250,2,'k.',flipud(hot(256)),indexREE-auxREEplot ,indexREE+auxREEplot,minREE(indexREE),maxREE(indexREE)); %axis square;xlim([minT maxT]);hold on;
hold on
end

xlim([0.5 13.5])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13],'XTickLabel',...
    {'La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Yb','Lu'});
ylabel('Source composition [ppm, CI normalized]')

%% FIGURE 5
% 1  2    3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18  19  20  21    22
% Tp ztop Al Na La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Yb Lu FeO MgO CaO C2rO3 Ti
Al2O3_ox = chain(burn:end,3);
Na2O_ox = chain(burn:end,4);
TiO2_ox = 0.13*ones(size(Na2O_ox));
MnO_ox  = 0.13*ones(size(Na2O_ox));
Cr2O3_ox = 0.57*ones(size(Na2O_ox));

[SiO2,Al2O3,FeOt,MnO,MgO,CaO,Na2O,Cr2O3,TiO2] = major_creator_new(Al2O3_ox,MnO_ox,Na2O_ox,Cr2O3_ox,TiO2_ox);

meltOX = [SiO2 Al2O3 FeOt MgO CaO Na2O Cr2O3];

auxOXplot = 0.25;
minOX = [40 2.5 5  35 2   0.1   0.1];   % requires adjustment; problem dependent
maxOX = [50 5 10 50   5  0.5  5 ];      % requires adjustment; problem dependent

figure(5);
clf(5);
for indexOX = 1:7
scattercloud_1D(indexOX*ones(1,nsimu),meltOX(:,indexOX),40,2,'k.',flipud(hot(256)),indexOX-auxOXplot ,indexOX+auxOXplot,minOX(indexOX),maxOX(indexOX)); %axis square;xlim([minT maxT]);hold on;
hold on
end
xlim([0.5 9.5])
ylim([0 55])
set(gca,'XTick',[1 2 3 4 5 6 7],'XTickLabel',...
    {'SiO_2','Al_2O_3','FeO','MgO','CaO','Na_2O','Cr_2O_3'});
ylabel('Source composition [wt%]')



