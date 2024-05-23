function [Maximum_c,Minimum_c] = find_directionally_coupled_complexes(model,B,NCPU,f_n)
% Function to test balancing capacity of a complex.
% 
% M_ij =  min / max A_i: * v
%           s.t.
%           AYv = Sv = b
%           A_j: * v = 0
%           lb <= v <= ub
%
% Input: model - model strcuture including complex-rxn matrix in model.A
%                and species-complex matrix in model.Y, 
%                other field required .S,.lb,.ub,.b
%        B     - set of balanced complexes, if not empty these will be
%                excluded from calculation of balancing effect to speed up
%                calculations
%        NCPU  - number of cores used in parallel loop
%
% Output: Maximum_c , Minimum_c - matrix M_ij with entry being the 
% maximum/minimum flux around a model complex j when i is balancing
%
if nargin < 3 || isempty(NCPU)
    NCPU = 1;
elseif NCPU > 1
    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= NCPU
        delete(p);
        parpool(NCPU);
    end
end

options = optimset('linprog');
options.Display = 'off';

%
Maximum_c=Inf(size(model.A,1),size(model.A,1));
Minimum_c=-Inf(size(model.A,1),size(model.A,1));

checked = nan(size(model.A,1));
checked(B,:) = 1;
checked(:,B) = 1;

for i=1:size(model.A,1)

    clc
    disp(f_n)
    disp(strcat('progress in check coupling:',{' '},num2str((i/size(model.A,1)*100)),'%'))
    
    if NCPU > 1
        params.Threads = 1;
        parfor j = 1:size(model.A,1)
            
            %    disp(strcat('check coupling file',{' '}, num2str(f), ':',{' '},num2str((i/size(model.A,1)*100)),'%'))
            
            % check if we tested the pair already
            if isnan(checked(i,j)) && i~=j
                
                checked(i,j)=1;
                
                objective=model.A(j,:);
                
                [~,f_k,ExitFlag]=linprog(-objective,[],[],[model.S; model.A(i,:)],[model.b; 0],zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);
                
                if ExitFlag == 1
                    
                    Maximum_c(i,j) = -f_k;
                    
                    [~,f_k,ExitFlag]=linprog(objective,[],[],[model.S; model.A(i,:)],[model.b;0],zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);
                    
                    if ExitFlag == 1
                        
                        Minimum_c(i,j) = f_k;
                    else
                        Minimum_c(i,j) = -Inf;
                    end
                    
                else
                    Maximum_c(i,j) = Inf;
                end
            end
        end
    else
        
        for j=1:size(model.A,1)
           
            % check if we tested the pair already
            if isnan(checked(i,j)) && i~=j
                
                checked(i,j)=1;
                
                objective=model.A(j,:);
                
                [R.x,R.f_k,R.ExitFlag]=linprog(-objective,[],[],[model.S; model.A(i,:)],[model.b;0],zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);
                
                if R.ExitFlag == 1
                    
                    Maximum_c(i,j) = -R.f_k;
                    
                    [R.x,R.f_k,R.ExitFlag]=linprog(objective,[],[],[model.S; model.A(i,:)],[model.b;0],zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);
                    
                    if R.ExitFlag == 1
                        
                        Minimum_c(i,j) = R.f_k;
                    else
                        Minimum_c(i,j) = -Inf;
                    end
                    
                else
                    Maximum_c(i,j) = Inf;
                end
            end
        end
end
end

end
