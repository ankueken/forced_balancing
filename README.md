# forced_balancing
Code to reproduce results presented in the manuscipt of "Imposing multireaction dependencies as an innovative strategy to alter metabolic functions"
by A. Küken, D. Langary, A. Angeleska, Z. Nikoloski

# Requirements
Matlab and R Installation <br>
Cobra Toolbox Installation % (https://opencobra.github.io/cobratoolbox/stable/installation.html) <br>
F2C2 Installation (https://sourceforge.net/projects/f2c2/) <br>

We fitted the balancing potentials to a power law as well as four non-scale-free distributions, as described by Broido and Clauset.
(A. D. Broido, A. Clauset, Scale-free networks are rare. Nature Communications 2019 10:1 10, 1–10 (2019)) 

# Example
Run the example from folder Example\ to see a full workflow including: <br>
  (1) removal of blocked reactions<br>
  (2) calaculation of balanced complexes<br>
  (3) calculation of concordant complexes<br>
  (4) calculation of balancing potential<br>
  
The output contains the following variables:<br>
% model - model structure <br>
% T - table of model reactions<br>
% Qi - the balancing potential<br>
% Qi_triv - number of trivially balanced complexes<br>
% Qi_non_triv - number of non-trivially balanced complexes<br>
% Qi_non_triv_type_2 - number of non-trivially balanced complexes type II<br>
% B - index of balanced complexes in model<br>
% class_with_balanced - concordance modules include the module of balanced complexes

