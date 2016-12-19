function [time_tec,tec,tec_inv,tec_smooth,br] = inv_ion(data,tgd,ele_cut,file_name);

% INV_ION	Computes total ionospheric electron content by inverting
%               LG for IEC and receiver interfrequency bias.
%
%	[time_tec,tec,tec_inv,tec_smooth,br] = inv_ion(data)
%
%       Inputs:
%         data = data structure (from make_data)
%         tgd = vector of transmitter group delays (from get_tgd)
%         ele_cut = elevation cut-off angle (degrees)
%         file_name = output file name for TEC
%
%       Outputs:
%         time_tec = time vector for tec (hours)
%	  tec = integrated electron content from averaging (el/m^2) => BIASED!
%	  tec_inv = integrated electron content from inversion (TECU = 10^16 el/m^2)
%         tec_smooth = smoothed tec_inv (TECU)
%	  br = receiver interfrequency bias

% check that data is not empty
if (isempty(data))
  disp(['--> Empty data structure, quitting inv_ion']);
  time_tec = [];
  tec = []; tec_inv = []; tec_smooth = [];
  br = [];
  return
end

% define some constants (should be read from a file)
A = 40.3;           % ionosphere constant
f1 = 1.57542e9;     % L1 frequency
f2 = 1.2276e9;      % L2 frequency
c = 0.299792458e9;  % speed of light
l1 = c/f1;
l2 = c/f2;

decim = 3;          % decimation factor
data_min = 2;       % minimum number of data per epoch

% multiplication factors
K = (A * (f1^2-f2^2)) / (f1^2*f2*c);
F = (f1^2*f2*c) / (A * (f1^2-f2^2));

% get field names
names = char(fieldnames(data));

% initialize data matrix
M = [];

% for each satellite
disp(['--> Reading data structure...']);
j = 1;
PRNS = [];
for i = 1:size(names,1)

  % read data structure
  eval(['time_iec = data.' names(i,:) '(:,1);']); % time in hours UT
  eval(['iec = data.' names(i,:) '(:,3);']);      % slant electron content in el/m^2
  eval(['ele = data.' names(i,:) '(:,7);']);      % sat. elevation angle in degrees
  eval(['emf = data.' names(i,:) '(:,12);']);     % elevation mapping function
  eval(['lgu = data.' names(i,:) '(:,13);']);     % LG corrected from mean(PG-LG)

  % find prn number
  prn = str2num(names(i,4:length(names(i,:))));

  % find corresponding transmitter group delay
  I = find(tgd(:,1)==prn);

  % check if tgd found for each satellite
  if (isempty(I))
    disp(['WARNING: no tgd found for prn' num2str(prn)]);
  else
    % get tgd
    tgd_prn = tgd(I,2);

    % apply elevation cut-off
    J = find(ele > ele_cut);

    % check if there is remaining data above cut-off
    if (isempty(J) ~= 1)
      % calculate TEC: apply elevation mapping function and correct TGD
      tec_prn = ((lgu(J) - f2 * tgd_prn) * F) .* emf(J);
      % convert to TECU
      tec_prn = tec_prn .* 1e-16;
      % write data into single matrix
      I = j .* ones(length(J),1);
      svprn = prn .* ones(length(J),1);
      M = [M ; I time_iec(J) lgu(J) emf(J) iec(J) K./emf(J) tec_prn svprn];
      % increment list of prns
      PRNS = [PRNS; names(i,:)];
      % increment j
      j = j + 1;
    end
  end
end

% figure out data and unknowns
n_data = size(M,1);
t_vect = uniq(sort(M(:,2))); n_epochs = length(t_vect);
s_vect = uniq(sort(M(:,1))); n_prns = length(s_vect);
%n_unkn = n_epochs + n_prns + 1;
n_unkn = n_epochs + 1;
disp(['--> Found ' int2str(n_data) ' data and ' int2str(n_unkn) ' unknowns total']);

% decimate: helps solve rank deficiency problem (for A)...
disp(['--> Decimating data by ' num2str(decim) '...']);
j = 1;
for i=1:decim:length(t_vect)
  t_decim(j) = t_vect(i);
  j = j + 1;
end

% figure out live data and unknowns
disp(['--> Reordering data matrix (data_min = ' num2str(data_min) ')...']);
live_unkn = 0;
live_data = 0;
N = [];
A = [];
D = [];
dpe = [];
tec = []; time_tec = [];
i = 1; j = 1;
for ti=t_decim
  % find data for epoch t
  I = find(M(:,2)==ti);
  % need at least data_min data per epoch
  if (length(I) >= data_min)
    live_unkn = live_unkn + 1;
    live_data = live_data + length(I);
    epoch = i .* ones(length(I),1);      % epoch #: will be used to built design matrix
    lg_data = j:j-1+length(I);           % data #: will be used to built design matrix
    N = [N; M(I,:) lg_data' epoch];      % reorder M into N
    D = [D ; M(I,3)];                    % read corresponding LGU data
    dpe = [dpe ; length(I)];             % keep track of how many data per epoch
    tec = [tec ; mean(M(I,7))];          % weighted average of tec for all prns
    time_tec = [time_tec; ti];           % keep track of time
    i = i + 1;
    j = j + length(lg_data);
  end
end
disp(['--> Live: ' int2str(live_data) ' data and ' int2str(live_unkn) ' unknowns']);

% write design matrix and data vector
disp(['--> Writing design matrix and data vector...']);
SUBA = sparse(N(:,9),N(:,10),N(:,6)*1e16,live_data,live_unkn);
m = size(SUBA,1);
F2 = -f2 .* ones(m,1);
A = [SUBA F2];

% write covariance matrix
%C = zeros(n_unkn);

% append condition equations for smoothing
%smoothing = 

% solve least squares using pseudoinverse
disp(['--> Running least squares...']);
[n_data_live n_unkn_live] = size(A);
X = A \ D;

% extract unknowns
tec_inv = X(1:length(X)-1);
br = X(length(X));
disp(['--> Receiver interfrequency bias = ' num2str(br*1e9) ' ns']);

% produce smooth TEC
win_length = 900;  % window length in seconds
stp = 1;           % step in epochs
splint = 30;       % sampling interval of data in seconds
opt = 3;           % no detrend, no polynomial removal
tec_smooth = smooth_alltraces(tec_inv,time_tec,win_length,stp,splint,opt);

% write into file
fid = fopen(file_name,'w');
fprintf(fid, '%f %f %f %f\n', [time_tec tec tec_inv tec_smooth]');
fclose(fid);

% plot result
%figure;
%I = find(tec_smooth);
%plot(time_tec(I),tec_smooth(I),'k');
%xlabel('time (hours UT)');
%ylabel('Total Electron Content (10^{16} el/m^2)');
 
