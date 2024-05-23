clear

%% Requirements
% Cobra Toolbox Installation % (https://opencobra.github.io/cobratoolbox/stable/installation.html)
% F2C2 Installation (https://sourceforge.net/projects/f2c2/)
% R Installation (https://www.r-project.org/)
% change path to R installation in line 96
% add Cobra and F2C2 to your matlab path

%% the example 
% (1) removes blocked reactions
% (2) calaculate balanced complexes
% (3) calculate concordant complexes
% (4) calculate balancing potential
% output:
% model - model structure 
% T - table of model reactions
% Qi - the balancing potential
% Qi_triv - number of trivially balanced complexes
% Qi_non_triv - number of non-trivially balanced complexes
% Qi_non_triv_type_2 - number of non-trivially balanced complexes type II
% B - index of balanced complexes
% class_with_balanced - concordance modules include the module of balanced
% complexes


%%{
addpath('../Functions/')

name = 'tcacycle';
model = readCbModel([name '.xml']);

% further simplify the model to have smaler example
model.S([1:2 4 12 17],:) = 0;
model.S(:,[6 16]) = 0;
model.S(:,8) = model.S(:,8)-model.S(:,15);
model.S(:,[12 15])=0;
model.S(16,:)=0;
model.lb([5 9 11])=0;

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

model.lb(mini>=0)=0;

model=removeRxns(model,BLK);

[sol.x,sol.f,sol.stat,sol.output]=linprog(-model.c,model.S(model.csense=='L',:),model.b(model.csense=='L'),model.S(model.csense=='E',:),model.b(model.csense=='E'),model.lb,model.ub);

if abs(sol.f)<abs(solo.f)*0.5
    disp('No biomass due to removal of blocked reactions')
end

clear BLK sol solo 

% changeCobraSolver('glpk')
model=convertToIrreversible(model);

save(strcat(name,'.mat'),"model")

T=table(cell(length(model.rxns),1),cell(length(model.rxns),1),cell(length(model.rxns),1),'VariableNames',{'Reaction_abbreviation' 'Reaction_formula' 'Reversible'}); % create empty table


T.Reaction_abbreviation = model.rxns(1:length(model.rxns)); % reaction name abbreviation can be found in .rxns, long names in .rxnNames

for i=1:length(model.rxns) % write formula
    % find substrates and their stoichiometry
    substrates = strjoin(strcat(num2str(abs(model.S(model.S(:,i)<0,i))), '*', model.mets(model.S(:,i)<0)),' + ');
    % find products and their stoichiometry
    products = strjoin(strcat(num2str(model.S(model.S(:,i)>0,i)), '*', model.mets(model.S(:,i)>0)),' + ');
    % combine both 
    T.Reaction_formula{i} = strcat(substrates, ' -> ', products);

    % fill reversible column using if statement
    if (model.lb(i) < 0 & model.ub(i) > 0) % reversible reactions have positive upper bound and negative lower bound
        T.Reversible(i) = {'yes'};
    else
         T.Reversible(i) = {'no'};
    end
end
T

cd ../Functions/
pathToR = '"C:\Users\Anika\AppData\Local\Programs\R\R-4.3.0\bin\Rscript.exe"';

system(strjoin({pathToR 'get_AY_matrix.r',strcat('../Example/', name, '.mat')}));

cd ../Example/
% random model
A=importdata(strcat(name,'.A'));
model.A=A.data;
model.complexes=A.textdata(2:end);
Y=importdata(strcat(name,'.Y'));
model.Y=Y.data;
clear A Y pathToR mini maxi
save(strcat(name,'.mat'),"model")
%}

%%
model.rev = model.lb<0 & model.ub>0;

[fctable,BLK] = F2C2('glpk', ...
    CobraToF2C2(model));

%% find balanced complexes
[B, group] = find_balanced_complexes(model);

CC = sparse(zeros(length(model.complexes)));
CC(B,B) = 1;

%% species degree

% model.ub([5 9 16])=0;

species_degree_complexes = sum(model.Y~=0,2);
species_degree_reactions = sum(model.S~=0,2);

idx = find(species_degree_complexes==2);
for i=1:length(idx)
    tcc = find(model.Y(idx(i),:)~=0);
    if CC(tcc,tcc)~=1
        CC(tcc,tcc) = [0 -1;-1 0];
       % CC(logical(eye(size(CC)))) = 0;
    end
end

options = optimset('linprog');
options.Display = 'off';
        
disp('start FVA ...')
[mini_elementary,maxi_elementary] = linprog_FVA(model,0);

try
    disp('Start sampling ...')
    sample=[];
    for r=1:length(model.lb)
        clc
        disp('Start sampling ...')
        disp(size(sample,2))
        model_s=model;
        sv= mini_elementary(r) + (maxi_elementary(r)-mini_elementary(r)).*0.5;
        model_s.lb(r)=sv;
        model_s.ub(r)=sv;
        [X,~,ExitFlag]=linprog(model_s.c,[],[],model_s.S,model_s.b,model_s.lb,model_s.ub,options);
        if ExitFlag==1
            sample(:,end+1)=X;
        end
    end

    At_sample=nan(size(model.A,1));
    At1=round(model.A*sample,9);

    for i=1:size(model.A,1)
        disp(i/size(model.A,1))
        At2=repmat(At1(i,:),size(At1,1),1);
        At_sample(:,i)=std((At1+eps)./(At2+eps),0,2,'omitnan')./mean((At1+eps)./(At2+eps),2,'omitnan');
    end

    At_sample(logical(eye(size(At_sample))))=nan;
    [At_rows,At_cols]=find(abs(At_sample)<0.02);
    At = [At_rows At_cols];

    clear model_s At1 At2 At_rows At_cols r sv X ExitFlag 
catch
end

if ~exist('At')
    At=nchoosek(1:length(model.complexes),2);
end

CC_temp=CC;
CC_temp(B,:)=5;
CC_temp(:,B)=5;

[checked_rows,checked_cols]=find(triu(CC_temp,1)~=0);
At = setdiff(At,[checked_rows checked_cols; checked_cols checked_rows],'rows');

clear checked_rows checked_cols CC_temp

% find concordant complexes
coupling_pairs = find_coupling_complexes(model,At,group);


    for i = 1:size(coupling_pairs,1)
        if CC(coupling_pairs(i,1),coupling_pairs(i,2)) == 0
            CC(coupling_pairs(i,1),coupling_pairs(i,2)) = 2;
            CC(coupling_pairs(i,2),coupling_pairs(i,1)) = 2;
        end
    end

%% 
% *Concordance modules*
    
    % *Group mutually concordant complexes without balanced*

   [CP1(:,1),CP1(:,2)] = find(CC==-1);
   [CP2(:,1),CP2(:,2)] = find(CC==2);
   CP = [CP1; CP2];
   unclassified = unique(CP);

    class=[];
    while ~isempty(unclassified)
        i = unclassified(1);
        class{end+1} = i;
        [r,~]=find(ismember(CP,i));
        while length(unique(reshape(CP(r,:),[],1))) > length(class{end})
            class{end} = unique(reshape(CP(r,:),[],1));
            unclassified = setdiff(unclassified,class{end});
            i = class{end};
            [r,~]=find(ismember(CP,i));
        end
    end
    

    % *Group mutually concordant complexes with balanced*
   clear CP
   [CP(:,1),CP(:,2)] = find(CC~=0);
   unclassified = unique(CP);
    
    class_with_balanced=[];
    while ~isempty(unclassified)
        i = unclassified(1);
        class_with_balanced{end+1} = i;
        [r,~]=find(ismember(CP,i));
        while length(unique(reshape(CP(r,:),[],1))) > length(class_with_balanced{end})
            class_with_balanced{end} = unique(reshape(CP(r,:),[],1));
            unclassified = setdiff(unclassified,class_with_balanced{end});
            i = class_with_balanced{end};
            [r,~]=find(ismember(CP,i));
        end
    end

%% balancing

[Max_dc,Min_dc] = find_directionally_coupled_complexes(model,B,1,1);

Max_dc_no_f = Max_dc;
% do not consider concordant coupling pairs
for i=1:size(coupling_pairs,1)
    Max_dc_no_f(coupling_pairs(i,1),coupling_pairs(i,2))=Inf;
    Max_dc_no_f(coupling_pairs(i,2),coupling_pairs(i,1))=Inf;
end
[row_inx,col_inx]=find(Max_dc_no_f==0 & Min_dc==0);
directionally_coupled_pairs=[row_inx,col_inx]; % row fixed to 0, col is optimized

Qi = sum(Max_dc==0 & Min_dc==0,2); % Qi the number of complexes that become balanced
% complexes that are balanced complexes are not considered in that number
% trivial and non-trivial

for i=1:length(class_with_balanced)
    Max_dc_no_f(class_with_balanced{i},class_with_balanced{i}) = Inf;
end
Qi_non_triv = sum(Max_dc_no_f==0 & Min_dc==0,2);
% only non-trivial balancing

Qi_triv = Qi - Qi_non_triv;

% for the non-trivial balancing pairs we need to check the type
% (a) due to blocking of flux, (b) complex becomes a non-trivial balanced
% complex with in flux = out flux

[r_check, c_check] = find(Max_dc_no_f==0 & Min_dc==0);
Dc_non_triv_type_2 = Max_dc_no_f==0 & Min_dc==0; % all non-triv have entry 1
for r=1:length(r_check)
    objective = double(model.A(c_check(r),:)<0);
    [~,Rf] = linprog(-objective,[],[],[model.S; model.A(r_check(r),:)],[model.b;0],zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);
    
    if isempty(Rf) || abs(Rf)<1e-9
        Dc_non_triv_type_2(r_check(r),c_check(r)) = 0; % system is blocked - set value to 0 s.t. it does not count anymore
    end
end
Qi_non_triv_type_2 = sum(Dc_non_triv_type_2,2);

max(Qi_non_triv_type_2)
max(Qi_triv) % size concordance module -1 because of not counting the balancing complex itself

