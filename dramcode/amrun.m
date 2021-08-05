function [results,chain,misfit,misfit2,melt_erupt,melt_path,s2chain,boundError,error]=amrun(model,data,params,options)
% AMRUN   Adaptive Metropolis MCMC run based on the algorithm of Haario
%         et al. (2001)
%
% This function generates a MCMC chain using adaptation for a model defined
% by a user supplied sum-of-squares function and with additive i.i.d. Gaussian
% errors for the observations. Based on the matlab DRAM code by M. Laine.
%
%
% input:
%
% model.ssfun    =  sum-of-squares function, ss=ssfun(par,data),
%                   that returns  -2*log(p(y|par))
% model.priorfun =  prior "sum-of-squares", priorfun(par,params),
%                   that returns -2*log(p(par)),
%                   default: inline('0','x','params')
%
% data           =  extra argument for ssfun (to pass the data etc.)
%
% params.par0    =  initial parameter vector (a row vector)
% params.sigma2  =  initial/prior value for the Gaussian error variance
% params.n0      =  precision of sigma2 as imaginative observations
%                   if n0<0, no sigma2 update
% params.n       =  number of actual observations (for sigma2 update)
% params.bounds  =  2*npar matrix of parameter bounds
%                   default: [-Inf,Inf]
%
% options.nsimu  =  length of the chain
% options.qcov   =  proposal covariance matrix
%
% parameters for AM
% options.adaptint  = 10;  % how often to adapt, if zero, no adaptation
%
% output:
%
% results  structure that contains some info about the run
% chain    nsimu*npar MCMC chain
% s2chain  sigma� chain (if generated)


%% get values from the input structs
nsimu  = getpar(options,'nsimu',10000);
%... initial parameter vector
par0   = getpar(params,'par0'); par0=par0(:)'; % row vector
%... number of parameters
npar   = length(par0);
%... 2*npar matrix of parameter bounds
bounds = getpar(params,'bounds',(ones(npar,2)*diag([-Inf,Inf]))');
%... parameters for underlying grid
grid = getpar(params,'grid');
%... parameters for stencil
stencil = getpar(params,'stencil');
%... sum-of-squares function, ssfun(par,data),  -2*log(p(y|theta))
ssfun  = getpar(model,'ssfun');
%... prior "sum-of-squares", -2*log(p(theta))
priorfun = getpar(model,'priorfun',inline('0','x','params'));

%... chemical data
ree_dat   = getpar(data,'mu');  % data vector
ndat      = length(ree_dat);    % size of data vector
cova_ree  = getpar(data,'Lam'); % get the inv of covariance matrix

%%% parameters for AM
%... how often to adapt, if zero, no adaptation
adaptint = getpar(options,'adaptint',100);
%... multiplier for first non-adaptive interval
int_fac  = getpar(options,'int_fac',2);
%... burn-in period
burn  = getpar(options,'burn',50);
%... scale for adapting the propsal
adascale = getpar(options,'adascale',2.4/sqrt(npar));
%... small factor for covariace update
qcovadj  = getpar(options,'qcovadj',1e-5);

%... precision of sigma2 as imaginative observations
%    if n0<0, no sigma2 update
n0  = getpar(params,'n0',-1);
%... initial/prior value for the Gaussian error variance
sigma2 = getpar(params,'sigma2',1);
%... number of observations (needed for sigma2 update)
if n0>=0, n = getpar(params,'n'); end
%... proposal covariance
qcov = getpar(options,'qcov');
%... how often do we save the chain
save_chain = getpar(options,'save_chain');

printint  = getpar(options,'printint',500);
verbosity = getpar(options,'verbosity',0);
plot1 = getpar(options,'plot',0);


R       = chol(qcov);         %  Cholesky factor of proposal covariance
chain   = zeros(nsimu,npar);  % we store the chain here
misfit  = zeros(nsimu,1);     % we store the misfit here
misfit2 = zeros(nsimu,length(par0)+2);   % we store the misfit of all trials here (with paramaters)

s20 = 0;
if n0>=0
    s2chain = zeros(nsimu,1);   % the sigma2 chain
    s20 = sigma2;
else
    s2chain = [];
end

oldpar       = par0(:)';  % first row of the chain
[oldss,melt_1,mpath,flag_error,fail] = feval(ssfun,oldpar,ree_dat,cova_ree,grid,stencil);% first sum-of-squares
oldprior     = feval(priorfun,oldpar,params);
outodbounds  = 0;

% accept your first guess
if fail == 0
    acce            = 1;             %  how many accepted moves
    chain(1,:)      = oldpar ;
    misfit(1,:)     = oldss ;
    melt_erupt(1,:) = melt_1 ;
    melt_path{1}    = mpath;
    oldmelt         = melt_1;
    oldpath         = mpath;
    error(1,:)      = flag_error;
    boundError(1,:) = [oldpar,outodbounds];
else
    fprintf('ERROR!!!!! initial try failed!!!!!!!\n')
    keyboard
end

if s20>0
    fprintf('ERROR!!!!! check sigma2 i dramrun!!!!!!!\n')
    s2chain(1,:) = sigma2;
end

if(plot1==1)
    subplot(2,2,1);
    scatter(chain(1,3),chain(1,4),'*g');axis square;hold on %axis([1200,1600,20,100]);
    subplot(2,2,2);
    scatter(chain(1,1),chain(1,2),'*g');axis square;hold on %
    subplot(2,2,3);
    plot(ree_dat(1:13),'*g'); xlim([0,14]);axis square;hold on %
    plot(melt_1(1:13),'or'); xlim([0,14]);axis square;hold on %
    
    subplot(2,2,4);
    plot(ree_dat(15:22),'*g');axis square;hold on %
    plot(melt_1(15:22),'or');pause(.0001);axis square;hold on %
end

% covariance update uses these to store previous values
chaincov = []; chainmean = []; wsum = []; lasti = 0;

countfail=0;
jca=10;


%%%%%%%%%%%%%%%%%%%%  Simulation starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for isimu=2:nsimu
    if any(isimu==save_chain:save_chain:nsimu)
        baseFileName = sprintf('workspace_isimu_%d',isimu);
        fullFileName = fullfile('outputs', baseFileName);
        save(fullFileName);
    end
    if isimu/printint == fix(isimu/printint) % info on every printint iteration
        fprintf('isimu=%d, %d%% done, accepted: %d%%\n',...
            isimu,fix(isimu/nsimu*100),fix((acce/isimu)*100));
    end
    
    newpar = oldpar+randn(1,npar)*R;     % a new proposal
    
    %... correct proposals (before adaptation) to account for correlations on major oxides
    if npar>17  % if
        if isimu <= adaptint*int_fac
            [newpar(18),newpar(19),newpar(20)]=create_sample(newpar(3),bounds); %FeO,MgO,CaO as functions of Al
        end
    end
    accept = 0;
    outodbounds = 0;
    
    % check bounds
    if any(newpar<bounds(1,:)) | any(newpar>bounds(2,:))
        newss = Inf;
        outodbounds = 1;
        fprintf('proposal out of bounds 1 \n');
        
        newprior = 0;
        alpha12 = 0;
    else % inside bounds, check if accepted
        if newpar(4)+0.05>0.6*newpar(3)-0.2 % not really inside bounds...
            newss = Inf;
            outodbounds = 1;
            fprintf('proposal out of bounds 2 \n');
            
            newprior = 0;
            alpha12 = 0;
        else
            [newss, melt_1,mpath,flag_error,fail]    = feval(ssfun,newpar,ree_dat,cova_ree,grid,stencil);   % sum-of-squares
            
            if fail == 1
                countfail = countfail+1;
            end
            newprior = feval(priorfun,newpar,params); % prior ss
            alpha12 = min(log(1),(-0.5*(newss-oldss)/sigma2 -0.5*(newprior-oldprior)));
            if log(rand) < alpha12 % we accept
                accept   = 1;
                %jca=jca+1;
                acce     = acce+1;
                oldpar   = newpar;
                oldss    = newss;
                oldprior = newprior;
                oldmelt = melt_1;
                oldpath = mpath;
            end
        end
    end
    if(plot1==1);if(newss > 500); miscol=500; else; miscol=newss;end;h = colorbar; set(h, 'ylim', [0 500])
        subplot(2,2,1); scatter(newpar(3),newpar(4),7.5,miscol,'filled');axis square; colormap(jet);colorbar,hold on %axis([1200,1600,20,100]);
        subplot(2,2,2); scatter(newpar(1),newpar(2),7.5,miscol,'filled');pause(.0001);axis square;hold on %axis([1200,1600,20,100]);
        if(accept ==1 && jca==10);subplot(2,2,3); plot(melt_1(1:13)); xlim([0,14]);hold on;axis square;hold on; jca=10; %end %axis([1200,1600,20,100]);
            subplot(2,2,4); plot(melt_1(15:22));pause(.0001);axis square;hold on ;end%axis([1200,1600,20,100]);
    end
    
    chain(isimu,:)      = oldpar;
    misfit(isimu,:)     = oldss;
    misfit2(isimu,:)    = [newpar,newss,accept];
    melt_erupt(isimu,:) = oldmelt;
    melt_path{isimu}    = oldpath;
    error(isimu,:)      = flag_error;
    boundError(isimu,:) = [newpar,outodbounds];
    
    if(plot1==1)
        if accept ==1
            subplot(2,2,1);
            scatter(chain(isimu,3),chain(isimu,4),'or','LineWidth',1.);hold on %; axis([0,1600,20,100]);
            subplot(2,2,2);
            scatter(chain(isimu,1),chain(isimu,2),'or','LineWidth',1.);pause(.0001);hold on %axis([1200,1600,20,100]);
            %subplot(2,2,3);
            %plot(melt_1); xlim([0,14]);hold on
        end
    end
    
    % update the error variance sigma2
    if s20 > 0
        fprintf('ERROR!!!!! check sigma2 i dramrun!!!!!!!\n')
        sigma2  = 1./gammar_mt(1,1,(n0+n)./2,2./(n0*s20+oldss));
        s2chain(isimu,:) = sigma2;
    end
    
    % check sigma... must be always 1!
    if sigma2 > 1 | sigma2 < 1
        fprintf('ERROR!!!!! check sigma2 i dramrun!!!!!!!\n')
    end
    
    if(isimu >= adaptint*int_fac )
        if adaptint>0 && fix(isimu/adaptint) == isimu/adaptint
            % adapt the proposal covariances
            if verbosity, fprintf('adapting\n'); end
            % update covariance and mean of the chain
            [chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1, ...
                chaincov,chainmean,wsum);
            %chaincov
            %chainmean
            %pause
            lasti = isimu;
            [Ra,is] = chol(chaincov + eye(npar)*qcovadj);
            if is % singular cmat
                fprintf('Warning cmat singular, not adapting\n');
            else
                R = Ra*adascale;
            end
        end
    end
    
end

%... calculate covariance and mean of the chain
[chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1, ...
    chaincov,chainmean,wsum);

%% save results
results.class = 'MCMC results';
results.accepted=acce./nsimu;              % acceptance ratio
results.mean = chainmean;
results.cov  = chaincov;
results.qcov = R'*R;
results.R = R;
results.nsimu = nsimu;
results.adascale = adascale;
results.adaptint = adaptint;
results.bounds = bounds;




function y=getpar(options,par,default)
%GETPAR get parameter value from a struct
% options   options struct
% par       parameter value to extract from the struct
% default   default value if par is not a member of the options struct

if isfield(options,par)
    y = getfield(options,par);
elseif nargin>2
    y = default;
else
    error(sprintf('Need value for option: %s',par));
end
