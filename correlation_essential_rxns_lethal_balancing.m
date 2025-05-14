% check if number of differentially leathal balanced complexes is correlated to 

% difference in the number of essential reactions between healthy and
% cancer

clear
filesC = dir('Results_lethality/*cancer_balanced.mat');
filesH = dir('Results_lethality/*normal_balanced.mat');

for f=1:length(filesH)

    C=load(['Results_lethality/' filesC(f).name]);
    H=load(['Results_lethality/' filesH(f).name]);
    
    disp(filesC(f).name)
    disp(filesH(f).name)
    
    C.complexes = cell(length(C.model.complexes),1);
        for i=1:length(C.model.complexes)
            C.complexes{i,1} = strjoin(strcat(num2str(C.model.Y(C.model.Y(:,i)~=0,i)), '*', C.model.mets(C.model.Y(:,i)~=0) ),'+');
        end
    
    H.complexes = cell(length(H.model.complexes),1);
        for i=1:length(H.model.complexes)
            H.complexes{i,1} = strjoin(strcat(num2str(H.model.Y(H.model.Y(:,i)~=0,i)), '*', H.model.mets(H.model.Y(:,i)~=0) ),'+');
        end
    
    [joint_complexes, inx_joint_C, inx_joint_H] = intersect(C.complexes,H.complexes);
    length(joint_complexes)
    
    
    lethal_candidates(f) = length(find(C.Bio_after_balancing(inx_joint_C)==0 & H.Bio_after_balancing(inx_joint_H)>0.9*H.Bio_opt));
    
end


%% 

clear
filesC = dir('Results_lethality_rxns/*cancer_balanced.mat');
filesH = dir('Results_lethality_rxns/*normal_balanced.mat');

for f=1:length(filesH)

    C=load(['Results_lethality_rxns/' filesC(f).name]);
    H=load(['Results_lethality_rxns/' filesH(f).name]);
    
    num_essential_rxns_C(f,1) = length(find(C.Bio_after_balancing==0));
    num_essential_rxns_H(f,1) = length(find(H.Bio_after_balancing==0));

    [~,rxns_common_C_inx,rxns_common_H_inx] = intersect(C.model.rxns,H.model.rxns);

    num_essential_switching(f,1) = length(find(C.Bio_after_balancing(rxns_common_C_inx) == 0 & H.Bio_after_balancing(rxns_common_H_inx) > 0.9*H.Bio_opt));

end

num_cand_balancing = [250 352 391 89 797 548 454 534 343]';

corr(num_cand_balancing,num_essential_switching)

corr(num_cand_balancing, (num_essential_rxns_C./num_essential_rxns_H)*100)

corr(num_cand_balancing, num_essential_rxns_C-num_essential_rxns_H)
