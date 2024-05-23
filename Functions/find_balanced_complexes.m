function [B,group] = find_balanced_complexes(model)
% function to find balanced complexes
%
% B = find_balanced_complexes(model,threshold)
%
% Input
%   model: struct with at least following fields
%           .c      objective vectore
%           .S      stoichiometric matrix
%           .b      right-hand side vector
%           .A      complex-reaction matrix
%           .lb     lower bound on flux
%           .ub     upper bound on flux
%   threshold: value when to consider complex balanced (e.g. 1e-9)
%              (optional)
% 
% Output
%   B:  vector of indices of balanced complexes

options = optimset('linprog');
options.Display = 'off';

threshold=1e-9;

for i=1:size(model.A,1)
    disp(i/size(model.A,1))
    c = model.A(i,:)';
    
    [~,f,stat]=linprog(c*-1,[],[],model.S,model.b,model.lb,model.ub,options);
    if stat == 1
        Maximum(i) = f*-1;
    else
        Maximum(i) = nan;
    end
    
    [~,f,stat]=linprog(c,[],[],model.S,model.b,model.lb,model.ub,options);
    if stat == 1
        Minimum(i) = f;
    else
        Minimum(i) = nan;
    end
end

Minimum=round(Minimum*(1/threshold))/(1/threshold);
Maximum=round(Maximum*(1/threshold))/(1/threshold);

B = intersect(find(Maximum==0),find(Minimum==0));

group = repmat('U',length(Maximum),1);
group(intersect(find(Maximum==0),find(Minimum==0))) = 'B';
group(intersect(find(Maximum>0),find(Minimum>=0))) = 'P';
group(intersect(find(Maximum<=0),find(Minimum<0))) = 'N';

% save(['../Results/' files(f).name '_balanced_' num2str(start) '.mat'])

end