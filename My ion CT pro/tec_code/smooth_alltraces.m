function iec_smooth = smooth_alltraces(iec_new,t,win_length,stp,splint,opt);

% SMOOTH_ALLTRACES	Smooth all traces using average over sliding window
%                         iec_new = array of iec data (MUST be one segment only!)
%                         t = corresponding time vector
%			  win_length = window length in seconds
%			  stp = step in epochs
%                         splint = sampling interval of data in seconds
%                         opt = 1 for detrend and remove polynomial
%                               2 for detrend only
%                               3 for no detrend and no polynomial

nobs = size(iec_new,1);
nsat = size(iec_new,2);
iec_smooth = zeros(nobs,nsat);

% smooth trace for each satellite
for i = 1:nsat
   % find non-zeros
   I = find(iec_new(:,i));
   if (~isempty(I))
      P = [t(I) iec_new(I,i)];
      if (~isempty(P))
         % detrend then fit and retrieve a polynomial
         if (opt == 1)
           P1 = detrend(P(:,2));
           p=polyfit(P(:,1),P1(:,1),4);
           f=polyval(p,P(:,1));
           P2 = P1(:,1)-f;
         end
         % detrend only
         if (opt ==2) P2 = detrend(P(:,2)); end
         % no detrend, no poly
         if (opt == 3)
           P2 = P(:,2);
         end
         % smooth using sliding window
         wdt = win_length/splint; % epochs
         k = 1;
         j = 1;
         while (k < (length(P2)-wdt))
            F(j) = mean(P2(k:k+wdt));
            H(j) = t(k);
            k = k + stp;
            j = j + 1;
         end
         if (exist('F'))
            hmin = find(t==H(1));
            hmax = find(t==H(length(H)));
            iec_smooth(hmin:hmax,i) = F';
         end
      end
   end
end

