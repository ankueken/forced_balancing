clear

%% organism specific GEMs
result_directory = 'Results_directionally_coupled_organism_GEMs';

%% calculate balancing potential
files = dir([result_directory '\*.mat']);

% for each file in directory calculate balancing potential
for f=1:length(files)
    disp(f)
    full_filename = [files(f).folder '\' files(f).name];
    load(full_filename,'Results_balanced', 'coupling_pairs', 'B', 'Max_dc*', 'Min_dc', 'class_with_balanced')
    model = Results_balanced.MODEL_r{1};
    clearvars -except model coupling_pairs B Max_dc* Min_dc f Qi* Bin files class_with_balanced num_balanced zero_potential

    % What is the balancing potential?
    Qi{f} = sum(Max_dc==0 & Min_dc==0,2); 
    % Qi - the number of complexes that become balanced complexes after balancing an other complex Ci
    % balanced complexes are not considered in that number
    % we distinguish trivial and non-trivial balancing
    
    for i=1:length(class_with_balanced) % concordance modules, balanced complex module included
        Max_dc_no_f(class_with_balanced{i},class_with_balanced{i}) = Inf; % remove balancing within concordance module 
    end
    Qi_non_triv{f} = sum(Max_dc_no_f==0 & Min_dc==0,2); % only non-trivial balancing
    
    Qi_triv{f} = Qi{f} - Qi_non_triv{f}; 

    num_balanced(f,1) = length(B) % number of balanced complexes per model
    
    % for the non-trivial balancing pairs we need to check the type
    % (a) due to blocking of flux, (b) complex becomes a non-trivial balanced
    % complex with in flux = out flux
    options = optimset('linprog');
    options.Display = 'off';
    
    [r_check, c_check] = find(Max_dc_no_f==0 & Min_dc==0);
    Dc_non_triv_type_2{f} = (Max_dc_no_f==0 & Min_dc==0); % 1 if non triv coupling
    for r=1:length(r_check)
        objective = double(model.A(c_check(r),:)<0);
        [~,Rf] = linprog(-objective,[],[],[model.S; model.A(r_check(r),:)],[model.b;0],zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);
                        
        if isempty(Rf) || abs(Rf)<1e-9
            Dc_non_triv_type_2{f}(r_check(r),c_check(r)) = 0; % system is blocked
        end
    end
    clear r_check c_check Rf objective r i options model coupling_pairs class_with_balanced Max_dc* Min_dc
    Qi_non_triv_type_2{f} = sum(Dc_non_triv_type_2{f},2);

    zero_potential{f} = (sum(Qi{f}==0)-length(B))/(length(Qi{f}-length(B)))
    clear B
end
save('Results_balancing_potential_organism_GEMs.mat')


%% plot distribution

label_names={'\it{A. niger} iMA871',...
'\it{A. thaliana} AraCore',...
'\it{C. reinhardtii} iCre1355',...
'\it{E. coli} iJO1366',...
'\it{M. acetivorans} iMB745',...
'\it{M. barkeri} iAF692',...
'\it{M. musculus}',...
'\it{M. tuberculosis} iNJ661m',...
'\it{N. pharaonis}',...
'\it{P. putida} iJN746',...
'\it{T. maritima} iLJ478',...
'\it{S. cerevesiae} Yeast8'};

Qi_percentage = cellfun(@(x,y) sum(x>0),Qi);
Qi_percentage = Qi_percentage./(cellfun(@length,Qi)-num_balanced');
Qi_non_triv_percentage = cellfun(@(x) sum(x>0),Qi_non_triv);
Qi_non_triv_percentage = Qi_non_triv_percentage./(cellfun(@length,Qi)-num_balanced');
Qi_non_triv_type_2_percentage = cellfun(@(x) sum(x>0),Qi_non_triv_type_2);
Qi_non_triv_type_2_percentage = Qi_non_triv_type_2_percentage./(cellfun(@length,Qi)-num_balanced');

colormap('parula')
subplot(1,2,1)
bar([Qi_percentage;Qi_non_triv_percentage;Qi_non_triv_type_2_percentage]'*100,'grouped')
legend('balancing potential', 'non-trivial balancing potential',...
    'non-trivial balancing potential of type II')
set(gca, 'XTickLabel',label_names,'Box','off')
ylabel({'Percent of model complexes'; 'with non-zero balancing potential'})
legend boxoff
ylim([0 100])
set(gca,'colororder',parula(3))

R=table(label_names',Qi_percentage'*100, Qi_non_triv_percentage'*100, Qi_non_triv_type_2_percentage'*100, 'VariableNames',{'model' 'Qi' 'Qi_non_triv' 'Qi_non_triv_type_II'})
