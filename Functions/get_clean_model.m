function get_clean_model(f)
disp(f)

files = dir('cancer_models/*.xml'); %update according to model folder
name = files(f).name(1:end-4)

% change to COBRA path
disp('set cobra path')
COBRA_PATH = '/work/ankueken/Git/cobratoolbox/';
addpath(genpath(COBRA_PATH));

model = readCbModel(strcat(files(f).folder,'/',files(f).name));

disp('Done read model')
if ~isfield(model,'csense')
    model.csense = repmat('E',size(model.mets));
end

disp('Clean model')
model=removeRxns(model,model.rxns(find(all(model.S==0))));
model=removeMetabolites(model,model.mets(find(all(model.S'==0))));

[solo.x,solo.f,solo.stat,solo.output]=linprog(-model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub);

[mini,maxi] = linprog_FVA(model,0.001);
thr=1e-9;
BLK=model.rxns(find(abs(mini)<thr & abs(maxi)<thr));

model=removeRxns(model,BLK);

[sol.x,sol.f,sol.stat,sol.output]=linprog(-model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub);

clear BLK sol solo files

model=convertToIrreversible(model);
cd Functions/
save(strcat('../cancer_models/cancer_models_pre_processed/final_cancer_models/' ,name, '.mat'),"model")

% pathToR = '"C:\Users\Anika\AppData\Local\Programs\R\R-4.3.0\bin\Rscript.exe"';

system(strjoin({pathToR 'get_AY_matrix.r',strcat('../cancer_models/cancer_models_pre_processed/final_cancer_models/', name,'.mat')}));

cd ../cancer_models/cancer_models_pre_processed/final_cancer_models/

A=importdata(strcat(name,'.A'));
model.A=sparse(A.data);
model.complexes=A.textdata(2:end);
clear A
Y=importdata(strcat(name,'.Y'));
model.Y=sparse(Y.data);
clear Y 
cd ../../../

save(['Results_balanced_cancer/pre_balanced_cancer/' name '_pre_balanced.mat'],'-v7.3')

end
