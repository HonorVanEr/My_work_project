function data = swap_filt(data,cutoff,f_type,splint);

% SWAP_FILT	Calculate new filterd time series and replaces
%		the corresponding entry in the 'data' structure.
%		  cutoff = filter cut-off period in minutes
%		  f_type = filter type (h/b/l)
%		  splint = sampling interval in seconds
%
%		data = swap_filt(data,cutoff,f_type,splint);
%

% get field names
fields = fieldnames(data);

% message
disp(['Swapping filtered signal using cutoff=' num2str(cutoff) ',type=' f_type ' ...']);

% for each field
for i = 1:length(fields)

  % read data
  s = cell2mat(fields(i));
  eval (['P = data.' s ';']);

  % read time and iec
  t = P(:,1);
  iec = P(:,3);

  % filter each trace
  iec_filt = filt_alltraces2(iec,t,splint,cutoff,f_type);

  % fill out data array
  eval (['data.' s '(:,2) = iec_filt;']);

end
