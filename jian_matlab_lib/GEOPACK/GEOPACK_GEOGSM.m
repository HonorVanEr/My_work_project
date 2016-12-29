function [A,B,C] = GEOPACK_GEOGSM (X,Y,Z,J)
% function [A,B,C] = GEOPACK_GEOGSM (X,Y,Z,J)
%      SUBROUTINE GEOGSM (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,J)
% C
% C CONVERTS GEOGRAPHIC (GEO) TO GEOCENTRIC SOLAR MAGNETOSPHERIC (GSM) COORDINATES
% C   OR VICA VERSA.
% C
% C                   J>0                   J<0
% C----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSM,YGSM,ZGSM
% C---- OUTPUT:    XGSM,YGSM,ZGSM      XGEO,YGEO,ZGEO
% C
% C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEOGSM IN TWO CASES:
% C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
% C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE BEEN CHANGED
% C
% C     LAST MODIFICATION: MARCH 31, 2003
% C
% C     AUTHOR:  N. A. TSYGANENKO
% C

%COMMON /GEOPACK1/AA(17),A11,A21,A31,A12,A22,A32,A13,A23,A33,D,B(8)
% C
if J < 0,
    [A,B,C] = GEOPACK_GSM2GEO(X,Y,Z);
else
    [A,B,C] = GEOPACK_GEO2GSM(X,Y,Z);
end

function [XGSM,YGSM,ZGSM] = GEOPACK_GEO2GSM(XGEO,YGEO,ZGEO);
global GEOPACK1;
XGSM=GEOPACK1.A11*XGEO+GEOPACK1.A12*YGEO+GEOPACK1.A13*ZGEO;
YGSM=GEOPACK1.A21*XGEO+GEOPACK1.A22*YGEO+GEOPACK1.A23*ZGEO;
ZGSM=GEOPACK1.A31*XGEO+GEOPACK1.A32*YGEO+GEOPACK1.A33*ZGEO;

function [XGEO,YGEO,ZGEO] = GEOPACK_GSM2GEO(XGSM,YGSM,ZGSM);
global GEOPACK1;
XGEO=GEOPACK1.A11*XGSM+GEOPACK1.A21*YGSM+GEOPACK1.A31*ZGSM;
YGEO=GEOPACK1.A12*XGSM+GEOPACK1.A22*YGSM+GEOPACK1.A32*ZGSM;
ZGEO=GEOPACK1.A13*XGSM+GEOPACK1.A23*YGSM+GEOPACK1.A33*ZGSM;

% end of GEOGSM
% C
% C=====================================================================================
% C

