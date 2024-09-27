%% Calculate forced balancing of complexes

clear
addpath('Functions\')

% Results of concordant complexes obtained as described in
% "The hidden simplicity of metabolic networks is revealed by multi-reaction dependencies"
% Küken et al. 2021 Sci. Adv.
% using the functions available at https://github.com/ankueken/concordant_complexes 
files = dir("Results_concordant_cancer\*.mat"); 

for f_n=1:length(files) 
    clearvars -except files f_n

    % load results on previously calculated balanced and concordant complexes
    load(strcat('Results_concordant_cancer/',files(f_n).name));
    
    [Max_dc,Min_dc] = find_directionally_coupled_complexes(Results_balanced.MODEL_r{1},B,1,1);
    
    save(strcat('Results_forced_balancing_cancer_GEMs/',name,'.mat'))
end

Max_dc_no_f = Max_dc;
% do not consider concordant coupling pairs
for i=1:length(coupling_pairs)
    Max_dc_no_f(coupling_pairs(i,1),coupling_pairs(i,2))=Inf;
    Max_dc_no_f(coupling_pairs(i,2),coupling_pairs(i,1))=Inf;
end
[row_inx,col_inx]=find(Max_dc_no_f==0 & Min_dc==0);
directionally_coupled_pairs=[row_inx,col_inx]; % row fixed to 0, col is optimized
save(strcat('Results_forced_balancing_cancer_GEMs/',name,'.mat'))


