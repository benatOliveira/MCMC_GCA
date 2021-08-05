function [newss] = misfit(data_REE,INVCOV,xx)

% function to compute the misfit (sum-of-squares) of a model. Possible compositional data is
% defined below... user can change these as required.

% new local vectors for manipulations
data1=data_REE(1:13);
data2=data_REE(14:22);
xx1=xx(1:13);
xx2=xx(14:22);
INVCOV1 = INVCOV(1:13,1:13);
INVCOV2 = INVCOV(14:22,14:22);


%% CALCULATE  MISFIT

% remove values
index_remove1 = data1==0; % index detecting the values to be removed
xx1(index_remove1) = [];
data1(index_remove1) = [];
INVCOV1(index_remove1,:) = [];
INVCOV1(:,index_remove1) = [];

index_remove2 = data2==0; % index detecting the values to be removed
xx2(index_remove2) = [];
data2(index_remove2) = [];
INVCOV2(index_remove2,:) = [];
INVCOV2(:,index_remove2) = [];

% compute misfit + weightings
newss_trace   = (xx1-data1)*INVCOV1*(xx1-data1)' ;
newss_majors  = (xx2-data2)*INVCOV2*(xx2-data2)' ;

newss = newss_trace + newss_majors*5;
end
