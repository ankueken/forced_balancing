% script to check how often forced balancing of a complex results into
% blocking of all input reactions, example E. coli
clear
load('Results_forced_balancing_organism_GEMs\Ecoli2011-iJO1366.mat','Results_balanced');
model = Results_balanced.MODEL_r{1}; clear Results_balanced

input_rxns = find(all(model.S>=0));
[input_mets,~] = find(model.S(:,input_rxns)>0);
output_rxns = find(all(model.S<=0));
[output_mets,~] = find(model.S(:,output_rxns)<0);
[~,~,output_mets_idx] = intersect(input_mets,output_mets);

options = optimset('linprog');
options.Display = 'off';
    
objective = zeros(size(model.S,2),1);
objective(input_rxns) = 1;
model.ub(output_rxns(output_mets_idx)) = 0;

for i=1:size(model.A,1) % for each complex force balancing and check maximum flux input reactions
    disp(i/size(model.A,1))

    [R.x,R.f_k,R.ExitFlag]=linprog(-objective,[],[],[model.S; model.A(i,:)],[model.b;0],zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);

    if ~isempty(R.f_k) && R.f_k<-1e-4
        Input_after_balancing(i,1) = 1;
    else
        Input_after_balancing(i,1) = 0;
    end

end





