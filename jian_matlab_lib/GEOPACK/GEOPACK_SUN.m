function [GST,SLONG,SRASN,SDEC] = GEOPACK_SUN(IYEAR,IDAY,IHOUR,MIN,ISEC)
% function [GST,SLONG,SRASN,SDEC] = GEOPACK_SUN(IYEAR,IDAY,IHOUR,MIN,ISEC)
%      SUBROUTINE SUN (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
% C
% C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
% C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
% C
% C-------  INPUT PARAMETERS:
% C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
% C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
% C
% C-------  OUTPUT PARAMETERS:
% C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
% C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
% C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
% C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
% C
% C  LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)
% C
% C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
% C

%      DOUBLE PRECISION DJ,FDAY
RAD = 57.295779513;
% C
if (IYEAR < 1901) | (IYEAR > 2099), return; end
FDAY=(IHOUR*3600+MIN*60+ISEC)/86400.D0;
DJ=365*(IYEAR-1900)+floor((IYEAR-1901)/4)+IDAY-0.5D0+FDAY;
T=DJ/36525.;
VL=rem(279.696678+0.9856473354*DJ,360.D0);
GST=rem(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD;
G=rem(358.475845+0.985600267*DJ,360.D0)/RAD;
SLONG=(VL+(1.91946-0.004789*T)*sin(G)+0.020094*sin(2.*G))/RAD;
if(SLONG > 6.2831853), SLONG=SLONG-6.2831853; end
if (SLONG < 0.), SLONG=SLONG+6.2831853; end
OBLIQ=(23.45229-0.0130125*T)/RAD;
SOB=sin(OBLIQ);
SLP=SLONG-9.924E-5;
% C
% C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
% C   THE ORBITAL MOTION OF THE EARTH
% C
SIND=SOB*sin(SLP);
COSD=sqrt(1.-SIND^2);
SC=SIND/COSD;
SDEC=atan(SC);
SRASN=3.141592654-atan2(cos(OBLIQ)/SOB*SC,-cos(SLP)/COSD);
%       RETURN
%       END
% end of function SUN
% C
% C================================================================================
% c
