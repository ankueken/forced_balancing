clear

%% organism specific GEMs
result_directory = 'Results_forced_balancing_cancer_GEMs';

%% calculate balancing potential
files = dir([result_directory '\*.mat']);

% for each file in directory calculate balancing potential
for f=1:length(files)
    disp(f)
    full_filename = [files(f).folder '\' files(f).name];
    load(full_filename)
    clearvars -except model B Max_dc Min_dc f Qi* files num_balanced zero_potential

    % What is the balancing potential?
    Qi{f} = sum(Max_dc==0 & Min_dc==0,2); 
    % Qi - the number of complexes that become balanced complexes after balancing an other complex Ci
    % balanced complexes are not considered in that number
    
    num_balanced(f,1) = length(B); % number of balanced complexes per model
   
    zero_potential{f} = (sum(Qi{f}==0)-length(B))/(length(Qi{f}-length(B)));
    clear B model Max_dc Min_dc
end
save('Results_balancing_potential_cancer_GEMs.mat')

%% Fit power-law and other distributions 

% use SFAnalysis
nbins=20;

for f=1:length(files)
    [BinCounts,BinEdges] = hist(Qi{f},1:2:length(Qi{f}),nbins);
    
    Bin.Count = BinCounts';
    Bin.Center = BinEdges';
    
    Bin_reduced.Center = log10(Bin.Center(Bin.Count>0)); % balancing potential > 0
    Bin_reduced.Count = Bin.Count(Bin.Count>0); % frequency
    Bin_reduced.Count = log10(Bin_reduced.Count/sum(Bin_reduced.Count));

     writetable(table(Bin.Center(Bin.Count>0),Bin.Count(Bin.Count>0),'VariableNames',{'xvalue','counts'}),...
         ['../SFAnalysis/example/degseqs/' files(f).name(1:end-3) 'txt'],'Delimiter',',')
  clear Bin*

end

%% compare cancer and healthy distribution parameters

[ks(1),p_ks(1)] = kstest2(Qi{1},Qi{2});
[ks(2),p_ks(2)] = kstest2(Qi{3},Qi{4});
[ks(3),p_ks(3)] = kstest2(Qi{5},Qi{6});
[ks(4),p_ks(4)] = kstest2(Qi{7},Qi{8});
[ks(5),p_ks(5)] = kstest2(Qi{9},Qi{10});
[ks(6),p_ks(6)] = kstest2(Qi{11},Qi{12});
[ks(7),p_ks(7)] = kstest2(Qi{13},Qi{14});
[ks(8),p_ks(8)] = kstest2(Qi{15},Qi{16});
[ks(9),p_ks(9)] = kstest2(Qi{17},Qi{18});

%% complexes whose balancing potential changes from healthy to cancer

complexes = cell(length(files),1);
for f=1:length(files)
    load(['Results_forced_balancing_cancer_GEMs\' files(f).name],'Max_dc','Min_dc','model','B')
    Qi{f} = sum(Max_dc==0 & Min_dc==0,2); 

    BALANCED{f} = B;
    MODELS{f} = model;
    complexes = cell(length(model.complexes),1);
    for i=1:length(model.complexes)
        complexes{i,1} = strjoin(strcat(num2str(model.Y(model.Y(:,i)~=0,i)), '*', model.mets(model.Y(:,i)~=0) ),'+');
    end
    COMPLEXES{f} = complexes;

    clear model B complexes
end

[~,inx{1},inx{2}]=intersect(COMPLEXES{1},COMPLEXES{2});
[~,inx{3},inx{4}]=intersect(COMPLEXES{3},COMPLEXES{4});
[~,inx{5},inx{6}]=intersect(COMPLEXES{5},COMPLEXES{6});
[~,inx{7},inx{8}]=intersect(COMPLEXES{7},COMPLEXES{8});
[~,inx{9},inx{10}]=intersect(COMPLEXES{9},COMPLEXES{10});
[~,inx{11},inx{12}]=intersect(COMPLEXES{11},COMPLEXES{12});
[~,inx{13},inx{14}]=intersect(COMPLEXES{13},COMPLEXES{14});
[~,inx{15},inx{16}]=intersect(COMPLEXES{15},COMPLEXES{16});
[~,inx{17},inx{18}]=intersect(COMPLEXES{17},COMPLEXES{18});

[sum(abs((Qi{1}(inx{1})+eps)./(eps+Qi{2}(inx{2}))-1)>=0.25)/length(inx{2})
    sum(abs((Qi{3}(inx{3})+eps)./(eps+Qi{4}(inx{4}))-1)>=0.25)/length(inx{4})
    sum(abs((Qi{5}(inx{5})+eps)./(eps+Qi{6}(inx{6}))-1)>=0.25)/length(inx{6})
    sum(abs((Qi{7}(inx{7})+eps)./(eps+Qi{8}(inx{8}))-1)>=0.25)/length(inx{8})
    sum(abs((Qi{9}(inx{9})+eps)./(eps+Qi{10}(inx{10}))-1)>=0.25)/length(inx{10})
    sum(abs((Qi{11}(inx{11})+eps)./(eps+Qi{12}(inx{12}))-1)>=0.25)/length(inx{12})
    sum(abs((Qi{13}(inx{13})+eps)./(eps+Qi{14}(inx{14}))-1)>=0.25)/length(inx{14})
    sum(abs((Qi{15}(inx{15})+eps)./(eps+Qi{16}(inx{16}))-1)>=0.25)/length(inx{16})
    sum(abs((Qi{17}(inx{17})+eps)./(eps+Qi{18}(inx{18}))-1)>=0.25)/length(inx{18})]

