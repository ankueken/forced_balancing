function lethality_check_essential_rxns(f)

files=dir('Results_forced_balancing_cancer_GEMs/*.mat');
load(['Results_forced_balancing_cancer_GEMs/' files(f).name],'model')

objective=zeros(size(model.c));
bio = find(contains(model.rxns,'bio','IgnoreCase',true));
if isempty(bio)
  disp('No biomass found')
else

  
  options = optimset('linprog');
  options.Display = 'off';
  objective(bio) = 1;
  
  [Bio.x,Bio.f_k,Bio.ExitFlag]=linprog(-objective,[],[],model.S,model.b,zeros(size(model.S,2),1),ones(size(model.S,2),1)*1e9,options);
  Bio_opt = -Bio.f_k;
  
    if Bio.f_k==0
      disp('No biomass')
    else
    
    Bio_after_balancing = zeros(size(model.A,1),1);
    
      for i=1:size(model.A,2) % for each rxn being blocked check maximum biomass 
      disp(i/size(model.A,2))
              
      ub = ones(size(model.S,2),1)*1e9;
      ub(i) = 0;
      [R.x,R.f_k,R.ExitFlag]=linprog(-objective,[],[],model.S,model.b,zeros(size(model.S,2),1),ub,options);
      
      if ~isempty(R.f_k)
        Bio_after_balancing(i,1) = -R.f_k;
      end       
    
    end      
          
end

save(['Results_lethality_rxns/' files(f).name])


end