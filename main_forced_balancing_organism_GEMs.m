% Calculate forced balancing of complexes

clear
addpath('Functions\')

species_list={'A_niger_iMA871';'ArabidopsisCoreModel';'M_acetivorans_iMB745';
    'Ecoli2011-iJO1366';'M_musculus';'M_barkeri_iAF692';'T_maritima_iLJ478';
    'C_reinhardtii_iCre1355_auto';'M_tuberculosis_iNJ661m';'P_putida_iJN746';
    'N_pharaonis';'YeastGEM'};

for f_n=1:length(species_list) 
    clearvars -except species_list f_n

    % load results on previously calculated balanced and concordant complexes
    load(strcat('Results_concordant_organism_GEMs/',species_list{f_n},'.mat'));
    
    [Max_dc,Min_dc] = find_directionally_coupled_complexes(Results_balanced.MODEL_r{1},B,1,1);
    
    save(strcat('Results_forced_balancing_organism_GEMs/',name,'.mat'))
end

Max_dc_no_f = Max_dc;
% do not consider concordant coupling pairs
for i=1:length(coupling_pairs)
    Max_dc_no_f(coupling_pairs(i,1),coupling_pairs(i,2))=Inf;
    Max_dc_no_f(coupling_pairs(i,2),coupling_pairs(i,1))=Inf;
end
[row_inx,col_inx]=find(Max_dc_no_f==0 & Min_dc==0);
directionally_coupled_pairs=[row_inx,col_inx]; % row fixed to 0, col is optimized
save(strcat('Results_forced_balancing_organism_GEMs/',name,'.mat'))


