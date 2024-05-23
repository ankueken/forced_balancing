function C = find_coupling_complexes(model,At,group)

options = optimset('linprog');
options.Display = 'off';

Maximum_c_p=Inf(size(model.A,1),size(model.A,1));
Minimum_c_p=-Inf(size(model.A,1),size(model.A,1));
Maximum_c_n=Inf(size(model.A,1),size(model.A,1));
Minimum_c_n=-Inf(size(model.A,1),size(model.A,1));

At_rows = At(:,1); At_cols = At(:,2);

fprintf('start coupling calculation for %i pairs \n', length(At_rows))

for i=1:length(At_cols)
    
    % check if we tested the pair already
    if At_rows(i)~=At_cols(i) && group(At_cols(i))~='N' 
              
        objective=model.A(At_rows(i),:);
        
        % maximize 
        [R.x,R.f_k,R.ExitFlag]=linprog([-objective 1],[eye(size(model.S,2)) -model.ub; -eye(size(model.S,2)) model.lb],[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At_cols(i),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; 0],[ones(size(model.S,2),1)*1e9; 999],options);
        
        if R.ExitFlag == 1
         
            Maximum_c_p(At_rows(i),At_cols(i)) = (model.A(At_rows(i),:)*R.x(1:end-1))/(model.A(At_cols(i),:)*R.x(1:end-1));
            
            % minimize
            [R.x,R.f_k,R.ExitFlag]=linprog([objective 1],[eye(size(model.S,2)) -model.ub; -eye(size(model.S,2)) model.lb],[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At_cols(i),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; 0],[ones(size(model.S,2),1)*1e9; 999],options);
            
            if R.ExitFlag ==1
                
                Minimum_c_p(At_rows(i),At_cols(i)) = (model.A(At_rows(i),:)*R.x(1:end-1))/(model.A(At_cols(i),:)*R.x(1:end-1));
            else
                Minimum_c_p(At_rows(i),At_cols(i)) = -Inf;
            end
        else
            Minimum_c_p(At_rows(i),At_cols(i)) = -Inf;
        end
    end
    if At_rows(i)~=At_cols(i) && group(At_cols(i))~='P'
        
        objective=model.A(At_rows(i),:);
                        
        [R.x,R.f_k,R.ExitFlag]=linprog([-objective 1],[-eye(size(model.S,2)) model.ub; eye(size(model.S,2)) -model.lb],[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At_cols(i),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; -999],[ones(size(model.S,2),1)*1e9; 0],options);
        
        if R.ExitFlag == 1
           
            Maximum_c_n(At_rows(i),At_cols(i)) = (model.A(At_rows(i),:)*R.x(1:end-1))/(model.A(At_cols(i),:)*R.x(1:end-1));
            
            [R.x,R.f_k,R.ExitFlag]=linprog([objective 1],[-eye(size(model.S,2)) model.ub; eye(size(model.S,2)) -model.lb],[zeros(size(model.S,2),1);zeros(size(model.S,2),1)],[model.S zeros(size(model.S,1),1); model.A(At_cols(i),:) 0],[model.b;1],[-ones(size(model.S,2),1)*1e9; -999],[ones(size(model.S,2),1)*1e9; 0],options);
            
            if R.ExitFlag == 1
               
                Minimum_c_n(At_rows(i),At_cols(i)) = (model.A(At_rows(i),:)*R.x(1:end-1))/(model.A(At_cols(i),:)*R.x(1:end-1));
            else
                Minimum_c_n(At_rows(i),At_cols(i)) = -Inf;
            end
            
        else
            Maximum_c_n(At_rows(i),At_cols(i)) = Inf;
        end
    end
end

Maximum_c_p(logical(eye(size(Maximum_c_p,1)))) = Inf;
Maximum_c_n(logical(eye(size(Maximum_c_n,1)))) = Inf;

for i=1:size(Maximum_c_p,1)
    for j=1:size(Maximum_c_p,1)
        if Maximum_c_p(i,j)~=Inf && Maximum_c_p(j,i)~=Inf
            Maximum_c_p(i,j)=Inf;
        end
    end
end

for i=1:size(Maximum_c_n,1)
    for j=1:size(Maximum_c_n,1)
        if Maximum_c_n(i,j)~=Inf && Maximum_c_n(j,i)~=Inf
            Maximum_c_n(i,j)=Inf;
        end
    end
end

Set_p = find(group=='P');Set_n = find(group=='N');Set_u = find(group=='U');

% for positive nominator check Maximum_c_p == Minimum_c_p
% save(['A_thaliana_concordant_all_' num2str(start) '.mat'])


[row_idx_c,col_idx_c] = find(round(Maximum_c_p(:,Set_p),2) == round(Minimum_c_p(:,Set_p),2));

Minimum_Set_p=Minimum_c_p(:,Set_p);

val_c2=[];
for i=1:length(row_idx_c)
    val_c2(i,1)=Minimum_Set_p(row_idx_c(i),col_idx_c(i));
end

CP=[row_idx_c,Set_p(col_idx_c),round(val_c2)];

% for negative nominator check Maximum_c_n == Minimum_c_n

[row_idx_c,col_idx_c] = find(round(Maximum_c_n(:,Set_n),2) == round(Minimum_c_n(:,Set_n),2));

Minimum_Set_n=Minimum_c_n(:,Set_n);

val_c2=[];
for i=1:length(row_idx_c)
    val_c2(i,1)=Minimum_Set_n(row_idx_c(i),col_idx_c(i));
end

CN=[row_idx_c,Set_n(col_idx_c),round(val_c2)];

% for unassigned nominator check Maximum_c_n == Minimum_c_n && Maximum_c_p == Minimum_c_p

linear_idx = intersect(intersect(find(round(Maximum_c_p(:,Set_u),2) == round(Minimum_c_p(:,Set_u),2)),...
    find(round(Maximum_c_n(:,Set_u),2) == round(Minimum_c_n(:,Set_u),2))),...
    find(round(Maximum_c_p(:,Set_u),2) == round(Minimum_c_n(:,Set_u),2)));

[row_idx_c,col_idx_c] = ind2sub(size(Minimum_c_p(:,Set_u)),linear_idx);
Minimum_Set1=Minimum_c_n(:,Set_u);

val_c1=[];
for i=1:length(row_idx_c)
    val_c1(i,1)=Minimum_Set1(row_idx_c(i),col_idx_c(i));
end

CU=[row_idx_c,Set_u(col_idx_c),round(val_c1,2)];


C = [CU; CP; CN];

end
