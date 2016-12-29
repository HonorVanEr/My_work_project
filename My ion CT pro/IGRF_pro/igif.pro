

Re=6371.2
hight=(indgen(40)*5+205.0)/Re+1

hight_km=(indgen(40)*5+205.0)
theta=INTARR(N_ELEMENTS(hight))+42.17
phi=INTARR(N_ELEMENTS(hight))+128

geopack_recalc,2007,198,6,30,00
GEOPACK_IGRF_GEO, hight, theta, phi, br, btheta, bphi,/DEGREE

B=SQRT(br^2+btheta^2+bphi^2)
;print,br, btheta, bphi,B

igrf={hight_km:hight_km,Br:br,B:B}


write_ascii_cmdline,igrf,'igrf_hight.txt'




end