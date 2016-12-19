function ele = check_ele(T,sip,ele_cutoff);

% CHECK_ELE	Checks elevation cutoff angle in SIP structure (greated by get_sip)
%		and returns a n_obs x n_sat array with:
%		- 0 when below cutoff
%		- 1 when above cutoff
%		This array can then be multiplied by the iec array to remove data
%		below the cutoff angle
%
%		Note that times have to be rounded because:
%		- time in sip structure comes from sp3 file
%		- time in T comes from rinex file

round_fact = 1000;
ele = [];

% message
disp(['Removing observations below ' num2str(ele_cutoff) ' deg. elevation';]);

% T = time vector in UT hours of current day
Tr = round(T.*round_fact);

% for each satellite
prns = fieldnames(sip);
for i = 1:length(prns)

  % extract the right field
  s = cell2mat(prns(i));         
  sv = getfield(sip,s);      

  % find the time when it is above the elevation cutoff angle
  I=find(sv(:,4)>ele_cutoff);
  Tok=round(sv(I,1) .* round_fact);

  % find corresponding time in full day time vector
  [C,Iok,IB] = intersect(Tr,Tok);

  % make E vector with 0 if ele below cutoff, 1 otherwise
  E = zeros(length(T),1);
  E(Iok) = 1;

  % increment output matrix
  ele = [ele E];

  % clean up
 % clear E C Iok IB s sv;
end
