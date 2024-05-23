clear
filesC = dir('Results_lethality/*cancer_balanced.mat');
filesH = dir('Results_lethality/*normal_balanced.mat');

%% find candidate complexes that are leathal in cancer but not in healthy tissue
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
        
    lethal_candidates = find(C.Bio_after_balancing(inx_joint_C)==0 & H.Bio_after_balancing(inx_joint_H)>0.9*H.Bio_opt);
        
    CANDIDATE_COMPLEXES{f,1} = C.complexes(inx_joint_C(lethal_candidates));
    [~,candidate_rxns{f,1}] = find(C.model.A(inx_joint_C(lethal_candidates),:)~=0);
    [~,candidate_rxns_H{f,1}] = find(H.model.A(inx_joint_H(lethal_candidates),:)~=0);
    CANDIDATE_EC{f,1} = C.model.rxnECNumbers(candidate_rxns{f,1});
    CANDIDATE_SYSTEM{f,1} = C.model.subSystems(candidate_rxns{f,1});

end

%% lethality - association to essential reactions?
clear
% changeCobraSolver('gurobi')
filesC = dir('Results_lethality/*cancer_balanced.mat');
filesH = dir('Results_lethality/*normal_balanced.mat');

filesCR = dir('Results_lethality_rxns/*cancer_balanced.mat');
filesHR = dir('Results_lethality_rxns/*normal_balanced.mat');

for f = 1:9
    C=load(['Results_lethality/' filesC(f).name]);
    H=load(['Results_lethality/' filesH(f).name]);
    
    disp(filesC(f).name)
    disp(filesH(f).name)
    
    C.complexes = cell(length(C.model.complexes),1);
        for i=1:length(C.model.complexes)
            C.complexes{i,1} = strjoin(strcat(num2str(C.model.Y(C.model.Y(:,i)~=0,i)), '*', C.model.mets(C.model.Y(:,i)~=0) ),'+');
        end
    C.complexes_long = cell(length(C.model.complexes),1);
        for i=1:length(C.model.complexes)
            C.complexes_long{i,1} = strjoin(strcat(num2str(C.model.Y(C.model.Y(:,i)~=0,i)), '*', C.model.metNames(C.model.Y(:,i)~=0) ),'+');
        end
    
    H.complexes = cell(length(H.model.complexes),1);
        for i=1:length(H.model.complexes)
            H.complexes{i,1} = strjoin(strcat(num2str(H.model.Y(H.model.Y(:,i)~=0,i)), '*', H.model.mets(H.model.Y(:,i)~=0) ),'+');
        end
    
    [joint_complexes, inx_joint_C, inx_joint_H] = intersect(C.complexes,H.complexes);
    length(joint_complexes)

    exchange_rxns = all(C.model.S>=0) + all(C.model.S<=0);
    [exchange_mets,~] = find(C.model.S(:,exchange_rxns~=0)~=0);
    exchange_mets=unique(exchange_mets);
        
    lethal_candidates = find(C.Bio_after_balancing(inx_joint_C)==0 & H.Bio_after_balancing(inx_joint_H)>0.9*H.Bio_opt);
        
    Crxns=load(['Results_lethality_rxns/' filesC(f).name]);
    Hrxns=load(['Results_lethality_rxns/' filesH(f).name]);
    
    disp(filesCR(f).name)
    disp(filesHR(f).name)

    concordant_model_C = load(['Results_concordant_cancer\' strrep(filesC(f).name,'balanced','concordant')]);
    concordant_model_H = load(['Results_concordant_cancer\' strrep(filesH(f).name,'balanced','concordant')]);
    CANDIDATE_COMPLEXES_no_essential{f,1} = [];
    for i = 1:length(lethal_candidates)

        modelC=concordant_model_C.Results_balanced.MODEL_r{1};
        modelC.S = [modelC.S; modelC.A(inx_joint_C(lethal_candidates(i)),:)];
        modelC.b = [modelC.b; 0];

        modelC.c(find(contains(modelC.rxns,'bio')))=1;

        sol_opt = optimizeCbModel(modelC);

        candidate_module = find(cellfun(@(x) any(x==inx_joint_C(lethal_candidates(i))),concordant_model_C.class_with_balanced));

        [~,candidate_rxns_C{f}] = find(Crxns.model.A(concordant_model_C.class_with_balanced{candidate_module},:)~=0);
        candidate_rxns_C{f} = unique(candidate_rxns_C{f});
        bio_candidate_rxns = Crxns.Bio_after_balancing(candidate_rxns_C{f});
        percentage_essential{f}(i) = sum(bio_candidate_rxns==0)/length(bio_candidate_rxns);

        if percentage_essential{f}(i) == 0 && isempty(intersect(exchange_mets,find(C.model.Y(:,inx_joint_C(lethal_candidates(i)))~=0)))
            CANDIDATE_COMPLEXES_no_essential{f,1}{end+1,1} = C.complexes_long(inx_joint_C(lethal_candidates(i)));
        end

        candidate_module = find(cellfun(@(x) any(x==inx_joint_H(lethal_candidates(i))),concordant_model_H.class_with_balanced));

        [~,candidate_rxns_H{f}] = find(Hrxns.model.A(concordant_model_H.class_with_balanced{candidate_module},:)~=0);
        candidate_rxns_H{f} = unique(candidate_rxns_H{f});
        
        bio_candidate_rxns_H = Hrxns.Bio_after_balancing(candidate_rxns_H{f});
        percentage_essential_H{f}(i) = sum(bio_candidate_rxns_H==0)/length(bio_candidate_rxns_H);
    end

    C1=sum(Crxns.Bio_after_balancing(candidate_rxns_C{f}) == 0);
    C2=sum(Crxns.Bio_after_balancing(candidate_rxns_C{f}) > 0);
    C3=sum(Crxns.Bio_after_balancing(setdiff(1:size(Crxns.model.A,2),candidate_rxns_C{f})) == 0);
    C4=sum(Crxns.Bio_after_balancing(setdiff(1:size(Crxns.model.A,2),candidate_rxns_C{f})) > 0);

    FT{f}=[C1 C2;C3 C4];
    [h{f},p{f}]=fishertest([C1 C2;C3 C4]);

    C1=sum(Hrxns.Bio_after_balancing(candidate_rxns_H{f}) == 0);
    C2=sum(Hrxns.Bio_after_balancing(candidate_rxns_H{f}) > 0);
    C3=sum(Hrxns.Bio_after_balancing(setdiff(1:size(Hrxns.model.A,2),candidate_rxns_H{f})) == 0);
    C4=sum(Hrxns.Bio_after_balancing(setdiff(1:size(Hrxns.model.A,2),candidate_rxns_H{f})) > 0);

    FT_H{f}=[C1 C2;C3 C4];
    [h_H{f},p_H{f}]=fishertest([C1 C2;C3 C4]);

end

% metabolites in candidate complexes (no essential reactions)
for i=1:9
    candidate_mets_no_essential{i}=[];
    for j=1:length(CANDIDATE_COMPLEXES_no_essential{i})
        if ~isempty(CANDIDATE_COMPLEXES_no_essential{i}{j}{1})
            CANDIDATE_COMPLEXES_no_essential{i}{j}{1} = strrep(CANDIDATE_COMPLEXES_no_essential{i}{j}{1},' ','');
            STRparts = strsplit(CANDIDATE_COMPLEXES_no_essential{i}{j}{1},'\+\d+','DelimiterType','RegularExpression');
            for k=1:length(STRparts)
                STRmets = strsplit(STRparts{k},'\d*\*','DelimiterType','RegularExpression');
                candidate_mets_no_essential{i}{end+1,1} = STRmets{2};
            end
        end
    end
    candidate_mets_no_essential{i} = unique(candidate_mets_no_essential{i});
end


% metabolites in candidate complexes composed of single metabolite (no essential reactions)
for i=1:9
    candidate_mets_single_no_essential{i}=[];
    for j=1:length(CANDIDATE_COMPLEXES_no_essential{i})
        if ~isempty(CANDIDATE_COMPLEXES_no_essential{i}{j}{1})
            CANDIDATE_COMPLEXES_no_essential{i}{j}{1} = strrep(CANDIDATE_COMPLEXES_no_essential{i}{j}{1},' ','');
            STRparts = strsplit(CANDIDATE_COMPLEXES_no_essential{i}{j}{1},'\+\d+','DelimiterType','RegularExpression');
            if length(STRparts)==1
                k=1;
                STRmets = strsplit(STRparts{k},'\d*\*','DelimiterType','RegularExpression');
                candidate_mets_single_no_essential{i}{end+1,1} = STRmets{2};
            end
        end
    end
    candidate_mets_single_no_essential{i} = unique(candidate_mets_single_no_essential{i});
end

entire_list = [candidate_mets_single_no_essential{1};candidate_mets_single_no_essential{2};candidate_mets_single_no_essential{3};candidate_mets_single_no_essential{4};...
    candidate_mets_single_no_essential{5};candidate_mets_single_no_essential{6};candidate_mets_single_no_essential{7};candidate_mets_single_no_essential{8};candidate_mets_single_no_essential{9}];
entire_list = unique(entire_list);

joint_mets_mat = zeros(length(entire_list),9);
for i=1:9
    if ~isempty(candidate_mets_single_no_essential{i})
        [~,inx] = intersect(entire_list,candidate_mets_single_no_essential{i});
    
        joint_mets_mat(inx,i) = 1;
    end
end

(cellfun(@(x) length(find(x==0)),percentage_essential)./cellfun(@length, percentage_essential))*100
cellfun(@mean,percentage_essential)
cellfun(@mean,percentage_essential_H)
cellfun(@median,percentage_essential)

