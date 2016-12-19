function [P1,P2,P3,F,H,h,g,pyy,w,mag] = my_filt(P,splint,cutoff,f_type);

% MY_FILT	Filter a times series using a 4th order Butterworth zero-phase filter.
%		- P = trace to be processed, N x 2 (at least) array with (time,amplitude)
%		- splint = sampling interval in seconds
%		- cutoff = cut-off period in minutes, can be [period_max period_min] for a pass-band filter
%		- f_type = filter type (h/b/l)
%
%               Output:
%		- P1 = detrended signal
%		- P2 = detrended - polynomial
%		- P3 = P2 * taper
%		- F  = filtered signal
%		- H  = time vector
%
%		[P1,P2,P3,F,H,h,g,pyy,w,mag] = my_filt(P,splint,cutoff,f_type);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the filter
cutoff=cutoff.*60;
nyqfrq=1.0/(2*splint);
Wn=(cutoff.*nyqfrq).^(-1);
ord=5;
if f_type == 'h'
        [b,a] = butter(ord,Wn,'high');
else
        [b,a] = butter(ord,Wn);
end

clear P1 P2 P3 F H h g pyy w mag;
P1=[]; P2=[]; P3=[]; F=[]; H=[]; h=[]; g=[]; pyy=[]; w=[]; mag=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A! ca peut planter si la taille de P est < a taille de b ou a
% This can crash if the size of P is < than the size of b or a
[m1,n1] = size(a);
[m2,n2] = size(P);
if (m2 > 2*n1+1)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% detrend the raw signal
	P1 = detrend(P(:,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fit and retrieve a polynomial
        [p,s] = polyfit(P(:,1),P1,4);
        f = polyval(p,P(:,1));
        P2 = P1-f;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% create cos taper
	[m n]=size(P(:,1));
	%w=hanning(m);
	% other taper: only the end and beginning (10%) are tapered
	m_10 = round(10*m/100);
	w1 = .5*(1 - cos(pi*(1:m_10)'/(m_10+1)));
	w2 = ones(m-(2*m_10),1);
	w3 = flipud(w1);
	w = [w1;w2;w3];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% mult. taper*signal
	P3 = P2 .* w;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% zero-phase filter of the massaged signal
	%F=filter(b,a,P3(:,1));    
	% Data must have length more than 3 times filter order
        nb = length(b);
        na = length(a);
        nfilt = max(nb,na);
        nfact = 3*(nfilt-1);  % length of edge transients

	if (length(P3) > nfact)
		F=filtfilt(b,a,P3);    
	else
		F=-999;P1=0;P2=0;P3=0;H=0;h=0;g=0;pyy=0;w=0;mag=0;
		return;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% cut the end and beginning
	P1=P1(m_10:m-m_10);
	P2=P2(m_10:m-m_10);
	F=F(m_10:m-m_10);
	H=P(m_10:m-m_10,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fourrier transform of the tapered signal
	window=1024;
	window_2=window/2;
	y=fft(P3,window);
	% compute power spectral density
	n=length(y);
	pyy=y.*conj(y)/n;
	% graph the result
	h=nyqfrq*(0:window_2)/window_2;
	pyy(window_2+2:window)=[];
	pyy(2:window_2) = 2*pyy(2:window_2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% frequency response of the filter
	[g,w] = freqz(b,a,128);
	mag = abs(g); phase = angle(g);
	% convert xaxis to Hz
	w = (w/pi)*nyqfrq;
end

% plot the frequency response  of the filter
%figure; 
%semilogy(w,mag); 
%loglog(w,mag); 
%plot(w,mag,':');  
%hold on;
%title('filter response');
%xlabel('frequency (Hz)');
%ylabel('amplitude (dB)');

