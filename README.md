# forced_balancing
Code to reproduce results presented in the manuscipt of "Imposing multireaction dependencies as an innovative strategy to alter metabolic functions"
by Anika Küken, Damoun Langary, Angela Angeleska, Zoran Nikoloski

# Requirements
Matlab and R Installation
Cobra Toolbox Installation % (https://opencobra.github.io/cobratoolbox/stable/installation.html)
F2C2 Installation (https://sourceforge.net/projects/f2c2/)

We fitted the balancing potentials to a power law as well as four non-scale-free distributions, as described by Broido and Clauset.
(A. D. Broido, A. Clauset, Scale-free networks are rare. Nature Communications 2019 10:1 10, 1–10 (2019)) 

Run the example to see a full workflow including
  (1) removal of blocked reactions
  (2) calaculation of balanced complexes
  (3) calculation of concordant complexes
  (4) calculation of balancing potential
  
The output contains the following variables:
% model - model structure 
% T - table of model reactions
% Qi - the balancing potential
% Qi_triv - number of trivially balanced complexes
% Qi_non_triv - number of non-trivially balanced complexes
% Qi_non_triv_type_2 - number of non-trivially balanced complexes type II
% B - index of balanced complexes in model
% class_with_balanced - concordance modules include the module of balanced complexes

