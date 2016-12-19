function get_tec(driver_file);
%driver_file = 'driver.038'
% GET_TEC	Main MATLAB program for computing ionospheric Total Electron Content
%               (TEC) from GPS phase and code data.
%
%               get_tec needs:
%                 - GPS observation data in rinex format
%                 - GPS navigation data in rinex format
%                 - precise orbits in sp3 format
%                 - a driver file that contains processing parameters
%
%               get_tec processes GPS data site by site and produces, for each site,
%               a structure containing:
%                   time_iec iec_filt iec_raw time_sip lat_sip lon_sip elev ...
%                   azim dist_epi azim_epi mag_los emf lgu
%
%               optionally, get_tec:
%                  - saves the site structures in .mat files
%                  - rearranges the site structures and produces one structure per
%                    satellite, same format as above
%                  - calculates absolute TEC and produces one output file per site
%                  - produces output files that can be plotted with GMT as
%                    hodograms (time - distance plots) and SIP maps
%
%		Compile readrinex inside Matlab first: 'mex readrinex.c'


% load add-on's by Thomas Dautermann
% these are files required for compute_position and and clean_data
% addpath add_on

% in case get_tec is not used as a function
%driver_file = 'driver.270';

disp('get_tec Version March 16th, 2007')

% remove matlab warnings
warning off MATLAB:polyfit:RepeatedPointsOrRescale;
warning off MATLAB:break_outside_of_loop;
warning off MATLAB:divideByZero;
warning off MATLAB:singularMatrix;

% this should be in driver file:
% filter_data='n';		% 'y' creates filtered trace as 2nd column in data structure,'n' uses raw IEC as 2nd colum
% cutoff =35;	          	% filter cut-off period in minutes
% f_type = 'h';			% filter type (h/b/l)
% make_hodo = 'n';		% prepare hodo data
% reorg_prn = 'n';                % reorganize data in structures by prn
% compute_tec = 'n';              % calculate absolute TEC
% save_as_mat = 'n';		% save site data structures in .mat files
% save_prn_asmat= 'n';		% save PRN data structure in mat files, requires reorg_prn='y'
% clean_data='n'; 		% filters oscillations in pseudorange and phase
% correct_IFB='n';		% after calculating tec (compute_tec=='y') correct traces for IFB
% write_tecfile='n';		% writes total TEC into textfile  site_day_year.tec, requires compute_tec='y'

% WGS-84 PARAMETERS
smj = 6378137.0;         % semi-major axis, m
smn = 6356752.314;       % semi-minor axis, m
fla = 1.0/298.257222101; % flattening 扁率
ecc = 2*fla - fla^2;     % eccentricity 离心率

% IEC Parameters
A = 40.3;           % ionosphere constant
f1 = 1.57542e9;     % L1 frequency
f2 = 1.2276e9;      % L2 frequency
c = 0.299792458e9;  % speed of light
K = (A * (f1^2-f2^2)) / (f1^2*f2*c);

% read driver file
fid = fopen(driver_file,'r');
if (fid==-1)
  disp(['ERROR: Cannot find driver file ' file]);
  return;
else
  n_rnxfiles = 0;
  while 1
    line = fgetl(fid);
    if ~isstr(line), break, end    %isstr可用ischar代替.
    if ~isempty(findstr(line,'comment')), comment = sscanf(line, '%s',1);, end
    if ~isempty(findstr(line,'time of event')), T_eq = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'latitude of event')), lat_epi = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'longitude of event')), lon_epi = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'latitude range')), latlim = sscanf(line, '%f %f');, end
    if ~isempty(findstr(line,'longitude range')), lonlim = sscanf(line, '%f %f');, end
    if ~isempty(findstr(line,'date')), date_epi = sscanf(line, '%d');, end
    if ~isempty(findstr(line,'elevation cutoff')), ele_cut = sscanf(line, '%d');, end
    if ~isempty(findstr(line,'orbit file')), sp3file = sscanf(line, '%s',1);, end
    if ~isempty(findstr(line,'ephemerides file')), efile = sscanf(line, '%s',1);, end
    if ~isempty(findstr(line,'ionex file')), ionexfile = sscanf(line, '%s',1);, end
    if ~isempty(findstr(line,'magnetic inclination')), inc = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'magnetic declination')), dec = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'ionospheric height')), alt_ion = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'maximum distance')), d_max = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'step for distance to source')), d_step = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'time window')), tmp = sscanf(line, '%f %f'); Tmin=tmp(1); Tmax=tmp(2);, end
    if ~isempty(findstr(line,'sampling interval')), splint = sscanf(line, '%f');, end
    if ~isempty(findstr(line,'filter data')), filter_data = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'filter type')), f_type = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'filter limit')), cutoff = sscanf(line, '%f %f');, end
    if ~isempty(findstr(line,'make hodo')), make_hodo = sscanf(line, '%c',1);, end   
    if ~isempty(findstr(line,'reorganize prn')), reorg_prn = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'compute tec')), compute_tec = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'save as mat')), save_as_mat = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'save prn as mat')), save_prn_as_mat = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'clean data')), clean_data = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'correct IFB')), correct_IFB = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'write TECfile')), write_tecfile = sscanf(line, '%c',1);, end
    if ~isempty(findstr(line,'data file'))
      rnxfile = sscanf(line, '%s',1);
      if (exist(rnxfile) ~= 0)
        n_rnxfiles = n_rnxfiles + 1;
        rnxlist(n_rnxfiles,:) = rnxfile;
      else
        disp(['WARNING: Cannot find rinex file ' rnxfile ', skipping...']);
      end
    end
  end
end
cutoff=1./(cutoff'.*1e-3)./60;

% message
disp(['======================']);
disp(['Driver file: ' driver_file]);
disp([comment ' ' num2str(date_epi)]);
disp(['Event at ' num2str(T_eq) ' TU, lat=' num2str(lat_epi) ', lon=' num2str(lon_epi)]);
disp(['Ionospheric height: ' num2str(alt_ion) ' meters']);
disp(['Mag. inclination/declination: ' num2str(inc) '/' num2str(dec) ' degrees']);
disp(['Elevation cutoff angle: ' num2str(ele_cut) ' degrees']);
disp(['Time window: ' num2str(Tmin) ' to ' num2str(Tmax) ' hours UT']);
disp(['Found ' num2str(n_rnxfiles) ' rinex file(s)']);
disp(['----------------------']);

% read sp3file and interpolate every splint seconds

% in old .sp3 files the satellites are identified by P 1, P 2,..
% in new .sp3 files the satellites are identified by PG01, PG02,...

% for  old sp3 use int_sp3(sp3file,splint)
% for  new  use int_sp3_new(sp3file,splint);
[SAT sp3_list] = int_sp3(sp3file,splint);

% read ephemerides file: get transmitter group delay
%tgd = get_tgd(efile);
tgd=get_tgd_ionex(ionexfile);
disp(['----------------------']);

%break

% process rinex files sequentially
SITE_LIST = []; SITE_COORD=[];
site_list = []; lat_sites = []; lon_sites = [];
for (irnx = 1:n_rnxfiles)

  % get rinex file name
  rnx_file = rnxlist(irnx,:);

  % read rinex data file
  [Observables,epochs,sv,apcoords]=readrinex(rnx_file);
  disp(['A priori Coordinates (XYZ,meters) = ' num2str(apcoords') ''])
  
  % check that data is actually sampled at splint for the whole file
  check_splint = check_sample(epochs,splint);
  if check_splint ~= 0
    disp(['WARNING: ' rnx_file ' not sampled at ' num2str(splint) ' sec, skipping...']);

  % check if station position is not zero (SHOULD BE COMPUTED)
  elseif sum(apcoords) == 0
    disp(['WARNING: ' rnx_file ' XYZ position not defined, skipping...']);

  % if ok, continue
  else

    % check for sv inconsistencies
    sp3_sv = [];
    exc_sv = [];
    for isv = 1:length(sv)
      if (find(sp3_list == sv(isv)))
        sp3_sv = [sp3_sv sv(isv)];
      else
        disp(['WARNING: PRN' num2str(sv(isv)) ' not in sp3 file']);
        exc_sv = [exc_sv sv(isv)];
      end
    end

    % get day of year from rinex file
    day = str2num(rnx_file(5:7));

    % get site name from rinex file name
    nm_site = rnx_file(1:4);
    SITE_LIST = [SITE_LIST ; nm_site];
    % get date from rinex file name
    this_date = [num2str(rnx_file(5:7)) '_' rnx_file(10:11)];

    % get coordinates of GPS site
    W=xyz2wgs([0 apcoords']);
    lat_site = W(3); lon_site = W(2); alt_site = 0;
    SITE_COORD=[SITE_COORD; [W(3) W(2)]];

    % compute distance and azim from event to site on WGS84 ellipsoid
    rng = distance(lat_epi,lon_epi,lat_site,lon_site,[smj ecc]);
    azi = azimuth(lat_epi,lon_epi,lat_site,lon_site,[smj ecc]);

    % message
    disp(['Processing GPS site ' nm_site ' [' num2str(irnx) '/' num2str(n_rnxfiles) '], lat=' num2str(lat_site) ' lon=' num2str(lon_site)])

    % organize observables into separate arrays
    L1=Observables.L1; L2=Observables.L2;
    P2=Observables.P2; C1=Observables.C1;
    if (isfield(Observables,'P1')) P1=Observables.P1; end;

    % replace NaN by zeros
    L1(isnan(L1)) = 0; L2(isnan(L2)) = 0;
    P2(isnan(P2)) = 0; C1(isnan(C1)) = 0;
    if (isfield(Observables,'P1'))
       P1(isnan(P1)) = 0;
    else
       P1 = [];
    end;

    % if P1 not defined in rinex or P1 defined but empty, then no P1 data => use C1 instead
    if (isempty(find(P1)))
      P1 = C1;
      disp(['WARNING: site ' nm_site ' does not have P1 data --> using C1']);
    end
    %filter pr and phase oscillations
    if clean_data=='y'
    	 disp('Cleaning Observables')
	 [L1,L2,C1,P1,P2]=clean_obs(L1,L2,C1,P1,P2);
    end %if clean_data

    % compute IEC and time derivative
    [IEC,IEC_DOT,LGU] = get_ion(L1,L2,C1,P1,P2,sv,tgd);

    % compute time vector (seconds of current day)
    T = epochs(:,4) + epochs(:,5)/60 + epochs(:,6)/3600;
  
    % select time interval
    I = find(T>=Tmin & T<=Tmax);
    t = T(I); iec = IEC(I,:); iec_dot = IEC_DOT(I,:); lgu = LGU(I,:);

    % replace NaN by zeros
    iec(isnan(iec)) = 0;
    iec(isnan(iec_dot)) = 0;
    iec(isnan(lgu)) = 0;

    % number of satellites and observations
    nsat = size(iec,2);
    nobs = size(iec,1);

    % get SIP information and create sip structure (only for svs present in sp3 file)
    sip = get_sip(SAT,sp3_sv,alt_ion,apcoords,Tmin,Tmax,lat_epi,lon_epi,dec,inc);

    % check data segmentation and keep longest segment only
    %iec_seg = chk_seg(iec);

    % correct for cycle slips, threshold = 1e20 el/m2
    %iec_slp = chk_slip(iec_seg,1e20);

    % filter trace for each satellite
    if 	filter_data=='y'
	%  filt_alltraces   has zeros at beginning and end 
	%  filt_alltraces2  has NaNs at beginning and end
    	iec_filt = filt_alltraces(iec,t,splint,cutoff,f_type);
    else
        iec_filt=iec;
    end

    % smooth trace for each satellite
    %iec_smooth = smooth_alltraces(iec_seg,t,500,1);

    % remove data below elevation cutoff angle
    n_bef = length(find(iec));
    ele = check_ele(t,sip,ele_cut);
    n_aft = length(find(ele));
    disp(['n_data: before = ' num2str(n_bef) ' / after = ' num2str(n_aft) ;]);
    tmp = iec .* ele; iec = tmp;
    tmp = iec_dot .* ele; iec_dot = tmp;
    tmp = iec_filt .* ele; iec_filt = tmp;
    tmp = lgu .* ele; lgu = tmp;

    % sort traces and create sr, sf, and sl structures (for all svs, but discard excluded)
    [sr,sf,sl,sv] = sort_traces(iec,iec_filt,lgu,t,sv,exc_sv);
    
    % at this point sv = final prn list

    % make final data structure
    data = make_data(sr,sf,sl,sip,sv);
    % data
    if filter_data=='y'
       data=swap_filt(data,cutoff,f_type,splint);
    end
    % correct data for receiver interfrequency bias (added 3-16-2007, thomas dautermann)
    ifb=ifbFromstruct(data,tgd);
    data=applyIfb(data,ifb);

    % assign final variable names
    eval(['data_' nm_site ' = data;']); eval(['sv_' nm_site ' = sv;']);

    % save data + lat/lon for each site in .mat file for future use
    if (save_as_mat == 'y')
      disp(['Saving TEC data as .mat file']);
      eval(['save ' nm_site '_' num2str(day) ' data_' nm_site ';']);
    end

    % keep site lat/lon
    %eval(['lat_' nm_site ' = lat_site;']); eval(['lon_' nm_site ' = lon_site;']);
    site_list = [site_list ; nm_site];
    lat_sites = [lat_sites; lat_site];
    lon_sites = [lon_sites; lon_site];

    % clean up
    %clear data iec iec_filt sip sr sf sl sv t;

  end
  disp(['----------------------']);
end

% calculate absolute TEC and write into file
if (compute_tec == 'y')
  disp([' ']);
  for j=1:size(SITE_LIST,1)
    site = SITE_LIST(j,:);
    disp(['Computing absolute TEC for ' site  '...']);
    eval(['data = data_' site ';']);
    file_name = [site '_' this_date '.tec'];
    [time_tec,tec,tec_inv,tec_smooth,br(j)] = inv_ion(data,tgd,ele_cut,file_name);
    eval(['tec_' site ' = [time_tec tec_smooth tec_inv];']);
    disp(['----------------------']);
   % data=setfield(data,'IFB',br(j));
    eval([' IFB_' site '=br(j);']);
    if write_tecfile=='y'
    	fid = fopen(file_name,'w');
    	fprintf(fid, '%f %f %f %f\n', [time_tec tec tec_inv tec_smooth]');
    	fclose(fid);
    end
  end
end

%correct traces for IFB -added by  Thomas Dautermann 3-24-06
if correct_IFB=='y'
disp('Applying IFB to IEC traces')
	for j=1:size(SITE_LIST,1)
		site = SITE_LIST(j,:);
		eval(['tmp_data = data_' site ';'])
			for i = 1:length(sp3_list)
				eval(['ok = isfield(tmp_data,''PRN' num2str(sp3_list(i)) ''');']);
				if (ok == 1)
					eval(['tmp_data.PRN' num2str(sp3_list(i)) '(:,3)=tmp_data.PRN' num2str(sp3_list(i)) ...
					'(:,3)+br(j)*f2/K;'])
				end
			end %for i
		eval(['data_' site  '=tmp_data;'])
	end %for j
end %if correct_IFB


% reorganize data in structures by prn number
if (reorg_prn == 'y')
  disp([' ']);disp(['Reorganizing data by prn number']);
  for i = 1:length(sp3_list)
    eval(['data_prn' num2str(sp3_list(i)) '=[];']);
    for j=1:size(SITE_LIST,1)
      site = switch1stLetter(SITE_LIST(j,:)); % 1st letter of field cannot be a number (added for Japanese sites)
      eval(['tmp_data = data_' site ';']);
      eval(['ok = isfield(tmp_data,''PRN' num2str(sp3_list(i)) ''');']);
      if (ok == 1)
        eval(['data_prn' num2str(sp3_list(i))  ' = setfield(data_prn' num2str(sp3_list(i)) ',site,tmp_data.PRN' num2str(sp3_list(i)) ');']);
	   if (save_prn_as_mat=='y')
	        eval(['save data_prn' num2str(sp3_list(i)) ' data_prn' num2str(sp3_list(i))])
 	   end
	end
    end
  end
  disp(['----------------------']);
end

% prepare data for hodo plot
if (make_hodo == 'y')
  disp([' ']);disp(['Creating dist-time-iec matrix for entire data set'])
  T=[]; D=[]; H=[];
  for i = 1:length(sp3_list)
    % check if this satellite is defined
    ok = eval(['exist(''data_prn' num2str(sp3_list(i)) ''');']);
    if (ok == 1)
      % check if structure for that satellite is empty
      ok = eval(['size(data_prn' num2str(sp3_list(i)) ',1);']);
        if (ok ~= 0)
          eval(['[T_tmp,D_tmp,H_tmp,LAT_tmp,LON_tmp] = plot_hodo(data_prn' num2str(sp3_list(i)) ',latlim,lonlim);']);
        end
      T = [T ; T_tmp]; D = [D ; D_tmp]; H = [H ; H_tmp];
    end
  end
  disp(['----------------------']);

  % prepare site_list to be used by plot_sip
  %eval(['site_list = fieldnames(data_prn' num2str(sp3_list(1)) ');']);
  %tmp = cell2mat(site_list);
  %site_list = tmp;
end

% save list of sites and their lat/lon for use by plot_sip
if (save_as_mat == 'y')
  disp(['Saving site names and coordinates']);
  save lat_sites lat_sites;
  save lon_sites lon_sites;
  save site_list site_list;
  disp(['----------------------']);
end
 plot_trace(data,'PRN12',T_eq);
 %[T,D,H,LAT,LON] = plot_hodo(data,latlim,lonlim);
%plot_data(data,factor,T_eq,threshold,plot_txt,decal);
disp(['======================']);

