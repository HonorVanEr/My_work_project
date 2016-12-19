function plot_trace(data,field,T_eq);

% PLOT_TRACE	Plots a single trace from 'field' given in structure 'data'
%
%		plot_trace(data,field,T_eq);
%

% check if field exists
if isfield(data,field)

  % initialize figure
  f = figure('Position',[296 100 560 725],'paperorientation','portrait','paperunits','inches','paperposition',[0 0 8.5 11],'papertype','usletter');

  % read data
  eval(['D = data.' field ';']);
  time_iec = D(:,1);
  iec_filt = D(:,2);
  iec_raw = D(:,3);
  time_sip = D(:,4);
  lat_sip = D(:,5);
  lon_sip = D(:,6);
  elev = D(:,7);
  azim = D(:,8);
  dist_epi = D(:,9);
  azim_epi = D(:,10);
  mag_los = D(:,11);

  % plot
  subplot(5,1,1); hold on; plot(time_iec,iec_raw);
  xlabel(['time (hr UT)']); ylabel(['raw TEC (el/m^2)']);

  subplot(5,1,2); hold on; plot(time_iec,iec_filt);
  xlabel(['time (hr UT)']); ylabel(['filt. TEC (el/m^2)']);

  subplot(5,1,3); hold on; plot(time_iec,elev);
  xlabel(['time (hr UT)']); ylabel(['elev. (deg.)']);

  subplot(5,1,4); hold on; plot(time_iec,azim);
  xlabel(['time (hr UT)']); ylabel(['azim. (deg.)']);

  subplot(5,1,5); hold on; plot(time_iec,mag_los);
  xlabel(['time (hr UT)']); ylabel(['mag/los angle (deg.)']);

else
  disp(['ERROR: field ' field ' not found in structure ' inputname(1) ]);
end
