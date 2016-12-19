function mag_los = get_mag_los(dec,inc,apcoords,S);

% GET_MAG_LOS  compute angle between magnetic field direction and line-of-sight
%              for a given site
%
%              dec = magnetic declination (degrees)
%              inc = magnetic inclination (degrees)
%              apcoords = ECEF coordinates of GPS site (XYZ, meters)
%                         (created by readrinex)
%              S = matrix of satellite positions (ECEF, meters)

% compute magnetic direction in ECEF frame at GPS site location
xyz_mag = elleza(dec,inc,1,apcoords); % UNIT VECTOR

% compute line-of-sight direction in ECEF frame at GPS location
ONES = ones(size(S,1),1);
A = [ONES .* apcoords(1) ONES.*apcoords(2) ONES.*apcoords(3) ];
XYZ_LOS = S - A;

% convert to unit vector
rang = sqrt(XYZ_LOS(:,1).^2+XYZ_LOS(:,2).^2+XYZ_LOS(:,3).^2);
XYZ_LOS_U = [XYZ_LOS(:,1)./rang XYZ_LOS(:,2)./rang XYZ_LOS(:,3)./rang];

% compute scalar product between magnetic direction and line-of-sight
XYZ_MAG = [ONES .* xyz_mag(1) ONES.*xyz_mag(2) ONES.*xyz_mag(3) ];
D = dot(XYZ_LOS_U,XYZ_MAG,2);

% convert cosine into mag/los angle
mag_los = acos(D) .* (180/pi);

