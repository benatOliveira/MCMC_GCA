% directory = '<Directory path>';      % Full path of the directory to be searched in
directory = [pwd,filesep,'perplex-model/Variables_Table_extra_aux'];
filesAndFolders = dir(directory);     % Returns all the files and folders in the directory
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory  
namesFiles = extractfield(filesInDir,'name')';
numOfFiles = length(filesInDir);

Al2O3_vector = [0.2:0.2:4.6];
Na2O_vector = [0.01 0.05:0.05:0.7];
Tp_vector = [1200:5:1600];             % [C]
noFile = [];
aux = 0;

for index_TP = 1:length(Tp_vector)
for index_Al2O3 = 1:length(Al2O3_vector)
for index_Na2O = 1:length(Na2O_vector)
    if Na2O_vector(index_Na2O)>0.6*Al2O3_vector(index_Al2O3)
        break
    end
    baseFileName = (sprintf('input_trace_smooth_isentrope_Na2O_%d_Al2O3_%d_Tp_%d.mat',index_Na2O,index_Al2O3,index_TP));
    ind = find(ismember(namesFiles, baseFileName));
    if isempty(ind)
        disp(['File ' baseFileName ', not found'])
        noFile = [noFile; index_Na2O,index_Al2O3,index_TP];
        aux = aux + 1;
    end
end
end
end
baseAux = sprintf('%d',aux);
disp([ baseAux ' files not found'])
