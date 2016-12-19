function iec_filt = filt_alltraces2(iec,t,splint,cutoff,f_type);

% FILT_ALLTRACES2	Filters all traces, segment by segment
% has NaN instead of zeros at the beginning and end
%                         iec = array of iec data (MUST be one segment only!)
%                         t = corresponding time vector
%			  splint = sampling interval in seconds
%			  cutoff = filter cut-off period in minutes
%			  f_type = filter type (h/b/l)

% number of satellites and observations
nobs = size(iec,1);
nsat = size(iec,2);

% message
% disp(['Filtering ' num2str(nsat) ' traces, splint=' num2str(splint) ', cutoff=' num2str(cutoff) ', f_type=' f_type ;]);

% initialize output array (changed from zeros to NaN by thomas dautermann 3-31-2006)
iec_filt =NaN*zeros(nobs,nsat);

% for each satellite
for i = 1:nsat
   % find non-zeros
   I = find(iec(:,i));
   if (~isempty(I))
      % detect segments
      x = [0 ; iec(:,i)./(iec(:,i)+eps^2)];
      y = [iec(:,i)./(iec(:,i)+eps^2) ; 0];
      z = x - y;
      J = find(z);
      %disp(['Found ' num2str(length(J)+1) ' segments;']);
      % if there is at least one segment
      if (~isempty(J))
         % for each segment
         for j = 1:2:length(J)
            % assign segment
            P = [t(J(j):J(j+1)-1) iec(J(j):J(j+1)-1,i)];
            if (~isempty(P))
               % filter segment
               F = []; P1 = []; P2 = []; P3 = [];
               [P1,P2,P3,F,H,h,g,pyy,w,mag] = my_filt(P,splint,cutoff,f_type);
               if (~isempty(F))
                  hmin = find(t==H(1));
                  hmax = find(t==H(length(H)));
                  iec_filt(hmin:hmax,i) = F;
               end
            end
         end
      end
   end
end 
