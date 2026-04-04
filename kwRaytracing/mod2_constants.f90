
!***********************************************************************
!This is code of YNOGK used for ray-traycing with General Relativity
!***********************************************************************

!************************************************************
module constants
    !********************************************************************
    !*    This module defines many constants often uesd in our code.
    !*    One can use these constants through a command "use constants" in their
    !*    own subroutines or functions.
    !********************************************************************
    implicit none
    Double precision  infinity,pi,dtors,sixteen,twopi,zero,one,two,three,four,six,half,&
    half2,mh,hbar,pho_v,plankc,five,dtor,eight
    parameter(infinity=1.D40,dtors=asin(1.D0)*2.D0/180.D0, &
    sixteen=16.D0, twopi=4.D0*dasin(1.D0), pi = dasin(1.D0)*2.D0)!3.141592653589793D0
    PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0, six=6.D0, half=0.5D0, half2=0.25D0, &
    mh=1.6726231D-24, hbar = 1.0545887D-27, plankc=6.626178D-27, pho_v=2.99792458D10, five=5.D0,&
    dtor=asin(1.D0)*2.D0/180.D0, eight=8.D0)
    !********************************************************************************************
end module constants
!********************************************************************************************
