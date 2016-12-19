function [SR,SF,SL,sv_new] = sort_traces(iec,iec_filt,lgu,t,sv,exc_sv);

% SORT_TRACES	Sorts iec, iec_filt, and lgu arrays and return structures:
%                 SF = filtered IEC
%                 SR = unfiltered IEC
%                 SL = LG corrected for ambiguity bias but not for interfequency bias
%		Also returns sv_new, new vector with prn numbers
%               Needs the list of 'live' (sv) and excluded (exc_sv) PRNs

% number of satellites and observations
nsat = size(iec,2);
nobs = size(iec,1);

% initialize structures
SF = [];
SR = [];
SL = [];

% write structures: [time data]
sv_new = [];
for i = 1:length(sv)

  % do it if this sv is not to be excluded
  if(isempty(find(sv(i)==exc_sv)))

    % get rid of zeros
    I = find(iec_filt(:,i));
    % fill out SF structure
    eval(['F' num2str(sv(i)) '= [t(I) iec_filt(I,i)];']);
    if (~isempty(eval(['F' num2str(sv(i))])))
       field = ['PRN' num2str(sv(i))];
       eval(['SF = setfield(SF,field,F' num2str(sv(i)) ');']);
       % use this opportunity to write sv vector
       sv_new = [sv_new sv(i)];
    end

    % get rid of zeros
    I = find(iec(:,i));
    % fill out SR structure
    eval(['R' num2str(sv(i)) '= [t(I) iec(I,i)];']);
    if (~isempty(eval(['R' num2str(sv(i))])))
       field = ['PRN' num2str(sv(i))];
       eval(['SR = setfield(SR,field,R' num2str(sv(i)) ');']);
    end

    % get rid of zeros
    I = find(lgu(:,i));
    % fill out SL structure
    eval(['L' num2str(sv(i)) '= [t(I) lgu(I,i)];']);
    if (~isempty(eval(['L' num2str(sv(i))])))
       field = ['PRN' num2str(sv(i))];
       eval(['SL = setfield(SL,field,L' num2str(sv(i)) ');']);
    end
  end
end

