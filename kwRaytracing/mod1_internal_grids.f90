module internal_grids
    integer nex,nec,nro,nphi
    parameter (nex=2**13) !Must be power of two for FFTs
    parameter (nec=100)  !Must be small for speed
    parameter (nro=500,nphi=500)
    real Emax,Emin,dloge,earx(0:nex),earc(0:nec)
    double precision aprev,mu0prev,pem1(nro,nphi),re1(nro,nphi)
    complex FTbbodx(4*nex)
    logical firstcall
    data firstcall/.true./
end module internal_grids
