function data = make_data(sr,sf,sl,sip,sv);

% MAKE_DATA	Merge sr, sf, sip, and lgu structures into a single structure 'data'
%                 sf = filtered IEC
%                 sr = unfiltered IEC
%                 sl = LG corrected for ambiguity bias but not for interfequency bias
%                 sip = 
%		Structure 'data' contains, for each prn:
%		   time_iec iec_filt iec_raw time_sip lat_sip lon_sip elev azim dist_epi azim_epi mag_los emf lgu
%		sv is a vector with the list of prn's to consider
%
%		data = make_data(sr,sf,sl,sip,sv);
%

% message
disp(['Merging iec and sip into single data structure']);

% initialize structure
data = [];

% get prn names
prn = char(fieldnames(sip));

% for each prn
for i=1:length(sv)

  % extract data
  eval(['tmp_raw = sr.PRN' num2str(sv(i)) ';']);
  eval(['tmp_fil = sf.PRN' num2str(sv(i)) ';']);
  eval(['tmp_sip = sip.PRN' num2str(sv(i)) ';']);
  eval(['tmp_lgu = sl.PRN' num2str(sv(i)) ';']);

  % extract time vector and round to the second
  tr = round(tmp_raw(:,1)*3600);
  tf = round(tmp_fil(:,1)*3600);
  ts = round(tmp_sip(:,1)*3600);
  tl = round(tmp_lgu(:,1)*3600);

  % find sip entries matching filtered iec (end of iec may not have sip values?)
  I = []; J = []; K = [];
  for j=1:length(tf)
    I = [I find(ts==tf(j))];
    J = [J find(tr==tf(j))];
    K = [K find(tl==tf(j))];
  end

  % rewrite tmp_fil
  tmp_fil = tmp_fil(1:length(I),:);

  % extract arrray from matching sip
  match_sip = tmp_sip(I,:);

  % extract vector from matching raw iec
  match_raw = tmp_raw(J,2);

  % extract vector from matching lgu
  match_lgu = tmp_lgu(K,2);

  % fill up data structure
  sf_sip = [tmp_fil match_raw match_sip match_lgu];
  eval(['data = setfield(data,''PRN' num2str(sv(i)) ''',sf_sip);']);

end
