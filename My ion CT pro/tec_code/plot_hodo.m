function [T,D,H,LAT,LON] = plot_hodo(data,latlim,lonlim);

% PLOT_HODO	Plot iec time series in 'true' [distance,time] domain
%               'data' is structure containing, for each subset, time and iec in columns 1 and 2
%		Note that data can be a 'site' structure or a 'prn' structure
%
%		Returns time (T), distance (D), hilbert transformed signal (H), latitude (LAT)
%		and longitude (LON) of the SIP.
%
%               [T,D,H,LAT,LON] = plot_hodo(data,latlim,lonlim);


%f=figure; hold on;
%xlabel('distance (km)');
%ylabel('local time (hours)');

% threshold to get rid of time series with abnormally large iec values
threshold = 1e16;

% initialize quantities to be gridded
T = []; D = []; H = []; LAT = []; LON = [];

names = char(fieldnames(data));

% for each field
for i=1:size(names,1)

  % read data structure
  eval(['time_iec = data.' names(i,:) '(:,1);']);
  eval(['iec = data.' names(i,:) '(:,2);']);
  eval(['lat = data.' names(i,:) '(:,5);']);
  eval(['lon = data.' names(i,:) '(:,6);']);
  eval(['dist_epi = data.' names(i,:) '(:,9);']);
  eval(['azim_epi = data.' names(i,:) '(:,10);']);

  if (max(iec) < threshold)

    % compute hilbert transfrom to get signal enveloppe
    hil = hilbert(iec);
    env = sqrt(real(hil).^2+imag(hil).^2);
    env = env*2;

    % append to time, distance, enveloppe
    T = [T;time_iec];
    D = [D;dist_epi];
    H = [H;env];
    LAT = [LAT;lat];
    LON = [LON;lon];

    %plot(dist_epi/1000, time_iec);
    %plot3(dist_epi/1000, time_iec, env);
    %plot(azim_epi*180/pi, time_iec);

  end

end

% find mins and maxs
[tmin,I] = min(T);
[tmax,J] = max(T);
dmin = D(I); lonmin = LON(I); latmin = LAT(I);
dmax = D(J); lonmax = LON(J); latmax = LAT(J);

% create map file for GMT
ll_name = ['llh_' inputname(1) '.dat'];
fid = fopen(ll_name,'w');
fprintf(fid, '%f %f %f\n', [LAT LON H]');

% create hodogram file for GMT
dt_name = ['dth_' inputname(1) '.dat'];
fid = fopen(dt_name,'w');
fprintf(fid, '%f %f %f\n', [D/1000 T H]');

% grid data over [time,distance] domain
% select near-field area
%max_dist = 200000;
%I = find (D < max_dist);
%T = T(I); D = D(I); H = H(I);
%step_dist = (max(D)-min(D))/100;
%step_time = 24/100;
%ti = [min(T):step_time:max(T)];
%di = [min(D):step_dist:max(D)];
%[DI,TI] = meshgrid(di,ti);
%HI = griddata(D,T,H,DI,TI,'cubic');
%pcolor(DI,TI,HI); shading('interp');

% prepare GMT run files
the_time = datestr(now);
run_name = ['run.' inputname(1) ];
fid = fopen(run_name,'w');

% PLOT MAP WITH GMT
fprintf(fid,'#! /bin/csh -f\n');
fprintf(fid,'# Created by plot_hodo on %s\n\n',the_time);
fprintf(fid,'gmtset BASEMAP_TYPE PLAIN\n');
fprintf(fid,'gmtset DEGREE_FORMAT 5\n');
fprintf(fid,'gmtset ANOT_FONT_SIZE 12\n');
fprintf(fid,'gmtset LABEL_FONT_SIZE 14\n');
fprintf(fid,'gmtset HEADER_FONT_SIZE 16\n\n');
fprintf(fid,'set lim = `minmax %s -C`\n',ll_name);
fprintf(fid,'set ncolor = `echo $lim[5] $lim[6] 200 | awk ''{print ($2-$1)/$3}''`\n');
fprintf(fid,'makecpt -T$lim[5]/3.5e15/$ncolor -Z >! color.cpt\n');
fprintf(fid,'set R = -R%.1f/%.1f/%.1f/%.1f\n',lonlim(1),lonlim(2),latlim(1),latlim(2));
fprintf(fid,'set S = -JM4i\n');
fprintf(fid,'set title = %s\n', inputname(1));
fprintf(fid,'pscoast $R $S -B1a4:.$title":WNse" -N1/1ta -W1/0 -P -K -Y5 >! %s.ps\n',inputname(1));
fprintf(fid,'psxy %s $R $S -: -Sc0.02 -Ccolor.cpt -O -K >> %s.ps\n',ll_name,inputname(1));
fprintf(fid,'pstext << eof $R $S -W0 -O -K >> %s.ps\n',inputname(1));
fprintf(fid,'%f %f 8 0 1 CM start\n%f %f 8 0 1 CM end\neof\n\n',lonmin,latmin,lonmax,latmax);

% PLOT HODOGRAM WITH GMT
%sort dth.dat -k 3 -o dth_sort.dat
fprintf(fid,'set lim = `minmax %s -C`\n',dt_name);
fprintf(fid,'set ncolor = `echo $lim[5] $lim[6] 200 | awk ''{print ($2-$1)/$3}''`\n');
fprintf(fid,'makecpt -T$lim[5]/3.5e15/$ncolor -Z >! color.cpt\n');
fprintf(fid,'set R = `minmax %s -I10/1`\n',dt_name);
fprintf(fid,'set S = -JX4i/4i\n');
fprintf(fid,'psxy %s $R $S -Ba400f50:"Distance (km)":/a1f0.2:"Time (hour UT)":nSWe -Sc0.02 -Ccolor.cpt -O -Y-4.5 >> %s.ps\n',dt_name,inputname(1));
fprintf(fid,'pstext << eof $R $S -W0 -O -K >> %s.ps\n',inputname(1));
fprintf(fid,'%f %f 8 0 1 CM start\n%f %f 8 0 1 CM end\neof\n\n',dmin,tmin,dmax,tmax);

fprintf(fid,'pstoimg %s.ps -antialias -crop a -density 150 -out %s.png -quiet -type png\n',inputname(1),inputname(1));
disp(['Created GMT file ' run_name ]);

%surface dth.dat -Gdth.grd -R0/400/0/24 -I0.24/4000
%grd2cpt dth.grd -Z > color.cpt
%grdimage dth.grd -R0/2400/0/24 -Ccolor.cpt -Jx0.25/0.00002 -B2a4/50000a100000::WenS -P >! junk.ps
