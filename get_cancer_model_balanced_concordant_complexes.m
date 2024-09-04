% calculate concordant complexes for cancer models
% cobra toolbox required here
clear

files = dir('cancer_models\*.xml');
addpath("Functions\")

for f=1:length(files)

    % step 1: preprocessing
    get_clean_model(f)

    % step 2: find balanced compelxes
    load(strcat('Results_balanced_cancer\pre_balanced_cancer\', files(f).name(1:end-4),'_pre_balanced.mat'));

    [B,group] = find_balanced_complexes(model);
    save(['Results_balanced_cancer\balanced_cancer\' files(f).name(1:end-4) '_balanced.mat'],'-v7.3')

    % step 3: find concordant complexes
    coupling_pairs = find_coupling_complexes(model,B,f,1e-9,group,'Results_concordant_cancer/', files(f).name(1:end-4));

    % to then inspect forced balancing --> see
    % main_forced_balancing_cancer_GEMs.m
end