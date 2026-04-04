! module blcoordinate
!   use blc_core_transforms, only: mutp, radiustp
!   use blc_disk_mapping,    only: pemdisk, pemdisk_all, mu2p, radius, mucos
!   use blc_integrals,       only: phi, ynogk, intrpart, inttpart
!   use blc_schwartz,        only: phyt_schwatz, schwatz_int, mu2p_schwartz
!   use blc_radius_extras,   only: r2p, rms, rmb, rph
!   use blc_geodesic_tools,  only: metricg, initialdirection, lambdaq, center_of_image
!   use blc_geodesic_tools,  only: p_total
!   implicit none
!   public
! end module blcoordinate

!*******************************************************************************
module blcoordinate
    !*******************************************************************************
    !*     PURPOSE: This module aims on computing 4 Boyer-Lindquist coordinates (r,\theta,\phi,t)
    !*              and affine parameter \sigam.
    !*     ACCURACY:   Machine.
    !*     AUTHOR:     Yang & Wang (2012)
    !*     DATE WRITTEN:  4 Jan 2012
    !***********************************************************************
    use constants
    use rootsfinding
    use ellfunction
    implicit none

    contains
    !********************************************************************************************
    !   SUBROUTINE YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
    !                     radi,mu,phi,time,sigma)
    SUBROUTINE YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,scal,&
        radi,mu,phi,time,sigma)
        !********************************************************************************************
        !*     PURPOSE:  Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine parameter
        !*               \sigma as functions of parameter p, i.e. functions r(p), \mu(p), \phi(p), t(p)
        !*               and \sigma(p). Cf. discussions in Yang & Wang (2012).
        !*     INPUTS:   p--------------independent variable, which must be nonnegative.
        !*               f1234----------array of p_1, p_2, p_3, p_4, which are the components of four-
        !*                              momentum of a photon measured under the LNRF frame. This array
        !*                              can be computed by subroutine lambdaq(...), see below.
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of observer or initialposition of photon.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  radi-----------value of function r(p).
        !*               mu-------------value of function \mu(p).
        !*               phi------------value of function \phi(p).
        !*               time-----------value of function t(p).
        !*               sigma----------value of function \sigma(p).
        !*               tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
        !*                              respectively.
        !*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
        !*                              respectively.
        !*     ROUTINES CALLED: INTRPART, INTTPART.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  5 Jan 2012
        !*     REVISIONS: ******************************************
        IMPLICIT NONE
        DOUBLE PRECISION f1234(4),lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,scal,&
        zero,one,two,three,four,phi_r,time_r,aff_r,phi_t,time_t,&
        mu_cos,r_coord,radi,mu,time,phi,sigma,p,Rab
        PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0)
        LOGICAL rotate,err
        INTEGER tm1,tm2,tr1,tr2

        !************************************************************************************
        ! call integrat_r_part to evaluate t_r,\phi_r,\sigma_r, and function r(p) (here is r_coord).
        ! call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
        ! scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
        call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,&
        scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
        ! call integrat_theta_part to evaluate t_\mu,\phi_\mu,\sigma_\mu, and function \mu(p) (here is mu_cos).
        ! call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,&
        !                 scal,phi_t,time_t,mu_cos,tm1,tm2)
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,lambdaBar,&
        scal,phi_t,time_t,mu_cos,tm1,tm2)
        radi=r_coord
        mu=mu_cos
        !time coordinate value, equation (74) of Yang & Wang (2012).
        time=time_r+time_t

        !time   write(*,*)time, time_r, time_t
        !affine parameter value, equation (74) of Yang & Wang (2012).
        !write(*,*)'ynogk=',aff_r,time_t,time_r,time_t
        sigma=aff_r+time_t
        !phi coordinate value.
        rotate=.false.
        err=.false.
        !write(*,*)'phi2=',p,phi_r,phi_t,tm1,tm2,time_r,time_t,lambda,f1234(3)
        IF(ABS(muobs).NE.ONE)THEN
            ! equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_t)
            IF(f1234(3).EQ.zero)THEN
                phi=phi+(tm1+tm2)*PI
            ENDIF
            phi=DMOD(phi,twopi)
            IF(phi.LT.zero)THEN
                phi=phi+twopi
            ENDIF
        ELSE
            ! equation (74) of Yang & Wang (2012).
            phi=-(phi_t+phi_r+(tm1+tm2)*PI)

            Rab=dsqrt(f1234(3)**two+f1234(2)**two)
            IF(phi.NE.zero)THEN
                rotate=.TRUE.
            ENDIF
            IF(Rab.NE.zero)THEN
                ! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=muobs*phi+asin(f1234(2)/Rab)
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=muobs*phi+PI-asin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=muobs*phi+PI-asin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=muobs*phi+twopi+asin(f1234(2)/Rab)
                endif
            ELSE
                phi=zero
            ENDIF
            IF(rotate)THEN
                phi=Mod(phi,twopi)
                IF(phi.LT.zero)THEN
                    phi=phi+twopi
                ENDIF
            ENDIF
        ENDIF
        RETURN
    END SUBROUTINE YNOGK

    !============================================================================================
    Function mucos(p,f12343,f12342,lambda,q,sinobs,muobs,a_spin,scal)
        !============================================================================================
        !*     PURPOSE:  Computes function \mu(p) defined by equation (32) in Yang & Wang (2012). That is
        !*               \mu(p)=b0/(4*\wp(p+PI0;g_2,g_3)-b1)+\mu_tp1. \wp(p+PI0;g_2,g_3) is the Weierstrass'
        !*               elliptic function.
        !*     INPUTS:   p--------------independent variable, which must be nonnegative.
        !*               f12342---------p_2, \theta component of four momentum of a photon measured under a LNRF.
        !*               f12343---------p_3, \phi component of four momentum of a photon measured under a LNRF..
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  mucos----------\mu coordinate of photon corresponding to a given p.
        !*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision mucos,f12343,f12342,p,sinobs,muobs,a_spin,q,lambda,mu_tp,&
        zero,b0,b1,b2,b3,g2,g3,tinf,fzero,a4,b4,AA,BB,two,scal,four,&
        mu_tp2,one,three,integ4(4),rff_p
        Double precision f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1
        complex*16 dd(3)
        integer ::  reals,p4,index_p4(4),del,cases,count_num=1
        logical :: mobseqmtp
        save  f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1,&
        mu_tp,b0,b1,b2,b3,g2,g3,dd,fzero,count_num,AA,BB,del
        parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0,three=3.D0)

        10    continue
        If(count_num.eq.1)then
            f12343_1=f12343
            f12342_1=f12342
            lambda_1=lambda
            q_1=q
            muobs_1=muobs
            sinobs_1=sinobs
            a_spin_1=a_spin
            scal_1=scal
            !*****************************************************************************************************
            If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
                mucos=muobs     !this is because that mu==1 for ever,this because that Theta_mu=-a^2(1-mu^2)
                return          !so,mu must =+1 or -1 for ever. q=-a^2, X=lambda/sin(theta)=0
            endif               !so Theta_mu=q+a^2mu^2-X^2mu^4=-a^2(1-mu^2)
            ! spin is zero.
            If(a_spin.eq.zero)then
                if(q.gt.zero)then
                    AA=sqrt((lambda**two+q)/q)
                    BB=sqrt(q)
                    If(f12342.lt.zero)then
                        mucos=sin(asin(muobs*AA)+p*BB*AA)/AA
                    else
                        If(f12342.eq.zero)then
                            mucos=cos(p*AA*BB)*muobs
                        else
                            mucos=sin(asin(muobs*AA)-p*AA*BB)/AA
                        endif
                    endif
                else
                    mucos=muobs
                endif
            else
                ! Equatorial plane motion.
                If(muobs.eq.zero.and.q.eq.zero)then
                    mucos=zero
                    return
                endif
                call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp)
                !   call mutp(f12342,f12343,sinobs,muobs,a_spin,lambdaBar,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp)
                a4=zero
                b4=one
                p4=0
                ! equations (26)-(29) in Yang & Wang (2012).
                b0=-four*a_spin**2*mu_tp**3+two*mu_tp*(a_spin**2-lambda**2-q)
                b1=-two*a_spin**2*mu_tp**2+one/three*(a_spin**2-lambda**2-q)
                b2=-four/three*a_spin**2*mu_tp
                b3=-a_spin**2
                ! equation (31) in Yang & Wang (2012).
                g2=three/four*(b1**2-b0*b2)
                g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)

                call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
                index_p4(1)=0
                cases=1
                ! equation (33) in Yang & Wang (2012).
                If(muobs.ne.mu_tp)then
                    tinf=b0/(four*(muobs-mu_tp))+b1/four
                    call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ4,cases)
                    fzero=integ4(1)
                else
                    fzero=zero
                endif
                If(f12342.lt.zero)then
                    fzero=-fzero
                endif
                ! equation (32) in Yang & Wang (2012).
                mucos=mu_tp+b0/(four*weierstrassP(p+fzero,g2,g3,dd,del)-b1)
                !write(*,*)'mu=',weierstrassP(p+fzero,g2,g3,dd,del),b0,mu_tp,b1!tinf,infinity,g2,g3,a4,b4,p4,fzero
                ! If muobs eq 0,q eq 0,and mu_tp eq 0,so b0 eq 0,
                ! so mucos eq mu_tp eq 0.
                count_num=count_num+1
            endif
        !**************************************************************************
        else
            If(f12343.eq.f12343_1.and.f12342.eq.f12342_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
                sinobs.eq.sinobs_1.and.muobs.eq.muobs_1.and.a_spin.eq.a_spin_1&
                .and.scal.eq.scal_1)then
                !******************************************************************
                If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
                    mucos=muobs                 !this is because that mu==1 for ever,this because that Theta_mu=-a^2(1-mu^2)
                    return                 !so,mu must =+1 or -1 for ever. q=-a^2, X=lambda/sin(theta)=0
                endif
                If(a_spin.eq.zero)then
                    if(q.gt.zero)then
                        If(f12342.lt.zero)then
                            mucos=sin(asin(muobs*AA)+p*BB*AA)/AA
                        else
                            If(f12342.eq.zero)then
                                mucos=cos(p*AA*BB)*muobs
                            else
                                mucos=sin(asin(muobs*AA)-p*AA*BB)/AA
                            endif
                        endif
                    else
                        mucos=muobs
                    endif
                else
                    If(muobs.eq.zero.and.q.eq.zero)then
                        mucos=zero
                        return
                    endif
                    ! equation (32) in Yang & Wang (2012).
                    mucos=mu_tp+b0/(four*weierstrassP(p+fzero,g2,g3,dd,del)-b1)
                !write(*,*)'mu=',mucos,weierstrassP(p+fzero,g2,g3,dd,del),mu_tp,b0,b1
                endif
            !********************************************************************
            else
                count_num=1
                goto 10
            endif
        endif
        return
    end Function mucos

    !********************************************************************************************
    !   Function radius(p,f1234r,lambda,q,a_spin,robs,scal)
    Function radius(p,f1234r,lambda,q,a_spin,lambdaBar,robs,scal)
        !============================================================================================
        !*     PURPOSE:  Computes function r(p) defined by equation (41) and (49) in Yang & Wang (2012). That is
        !*               r(p)=b0/(4*\wp(p+PIr;g_2,g_3)-b1)+r_tp1. \wp(p+PIr;g_2,g_3) is the Weierstrass'
        !*               elliptic function; Or r=r_+, r=r_-.
        !*     INPUTS:   p--------------independent variable, which must be nonnegative.
        !*               f1234r---------p_1, r components of four momentum of a photon measured under a LNRF.
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  radius---------radial coordinate of photon corresponding to a given p.
        !*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision radius,p,a_spin,lambdaBar,rhorizon,q,lambda,scal,zero,integ4(4),&
        cc,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,robs,cr,dr,integ04(4),&
        u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,a4,b4,one,two,four,sqt3,&
        integ14(4),three,six,nine,r_tp1,r_tp2,f1234r,tp2,t_inf,PI0_total,&
        PI0_inf_obs,PI0_obs_hori,PI01,PI0_total_2,rff_p
        Double precision f1234r_1,lambda_1,q_1,a_spin_1,lambdaBar_1,robs_1,scal_1
        parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
        complex*16 bb(1:4),dd(3)
        integer ::  reals,cases_int,del,index_p4(4),cases,count_num=1
        logical :: robs_eq_rtp,indrhorizon
        save  f1234r_1,lambda_1,q_1,a_spin_1,lambdaBar_1,robs_1,scal_1,r_tp1,r_tp2,reals,&
        robs_eq_rtp,indrhorizon,cases,bb,rhorizon,b0,b1,b2,b3,g2,g3,dd,del,cc,tinf,tp2,&
        thorizon,tinf1,PI0,u,w,v,L1,L2,m2,t_inf,pinf,a4,b4,PI0_total,PI0_inf_obs,PI0_obs_hori,&
        PI0_total_2

        20   continue
        If(count_num.eq.1)then
            f1234r_1=f1234r
            lambda_1=lambda
            q_1=q
            a_spin_1=a_spin
            lambdaBar_1=lambdaBar
            robs_1=robs
            scal_1=scal
            !*********************************************************************************************
            !   rhorizon=one+sqrt(one-a_spin**2)
            ! new horizon with lambdaBar (throat radius)
            rhorizon = (one + lambdaBar**2) + sqrt((one + lambdaBar**2)**2 - a_spin**2)

            a4=zero
            b4=one
            cc=a_spin**2-lambda**2-q
            robs_eq_rtp=.false.
            indrhorizon=.false.
            !   call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
            !                      reals,robs_eq_rtp,indrhorizon,cases,bb)
            call radiustp(f1234r,a_spin,lambdaBar,robs,lambda,q,r_tp1,r_tp2,&
            reals,robs_eq_rtp,indrhorizon,cases,bb)
            If(reals.ne.0)then
                ! equations (35)-(38) in Yang & Wang (2012).
                b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
                b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
                b2=four/three*r_tp1
                b3=one
                g2=three/four*(b1**2-b0*b2)
                g3=one/16.D0*(3*b0*b1*b2-2*b1**3-b0**2*b3)
                ! equation (39) in Yang & Wang (2012).
                If(robs-r_tp1.ne.zero)then
                    tinf=b0/four/(robs-r_tp1)+b1/four
                else
                    tinf=infinity
                endif
                If(rhorizon-r_tp1.ne.zero)then
                    thorizon=b1/four+b0/four/(rhorizon-r_tp1)
                else
                    thorizon=infinity
                endif
                tp2=b0/four/(r_tp2-r_tp1)+b1/four
                tinf1=b1/four

                call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
                index_p4(1)=0
                cases_int=1
                ! equation (42) in Yang & Wang (2012).
                call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                PI0=integ04(1)
                select case(cases)
                case(1)
                    If(.not.indrhorizon)then
                        If(f1234r.lt.zero)then
                            call weierstrass_int_j3(tinf1,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            PI0_total=PI0+integ14(1)
                            If(p.lt.PI0_total)then
                                ! equation (41) in Yang & Wang (2012).
                                radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)
                            else
                                radius=infinity  !Goto infinity, far away.
                            endif
                        else
                            call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            PI0_inf_obs=integ04(1)
                            If(p.lt.PI0_inf_obs)then
                                ! equation (41) in Yang & Wang (2012).
                                radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)
                            else
                                radius=infinity !Goto infinity, far away.
                            endif
                        endif
                    else
                        If(f1234r.lt.zero)then
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            PI0_obs_hori=integ04(1)
                            If(p.lt.PI0_obs_hori)then
                                ! equation (41) in Yang & Wang (2012).
                                radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)
                            else
                                radius=rhorizon !Fall into black hole.
                            endif
                        else
                            call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            PI0_inf_obs=integ04(1)
                            If(p.lt.PI0_inf_obs)then
                                ! equation (41) in Yang & Wang (2012).
                                radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)
                            else
                                radius=infinity !Goto infinity, far away.
                            endif
                        endif
                    endif
                case(2)
                    If(.not.indrhorizon)then
                        If(f1234r.lt.zero)then
                            PI01=-PI0
                        else
                            PI01=PI0
                        endif
                        ! equation (41) in Yang & Wang (2012).
                        radius=r_tp1+b0/(four*weierstrassP(p+PI01,g2,g3,dd,del)-b1)
                    else
                        If(f1234r.le.zero)then
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            PI0_obs_hori = integ14(1)
                            If(p.lt.PI0_obs_hori)then
                                ! equation (41) in Yang & Wang (2012).
                                radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)
                            else
                                radius=rhorizon !Fall into black hole.
                            endif
                        else
                            call weierstrass_int_J3(tp2,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            call weierstrass_int_J3(tp2,tinf,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            PI0_total_2=integ14(1)+integ4(1)
                            If(p.lt.PI0_total_2)then
                                ! equation (41) in Yang & Wang (2012).
                                radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)
                            else
                                radius=rhorizon !Fall into black hole.
                            endif
                        endif
                    endif
                end select
                If(a_spin.eq.zero)then
                    If(cc.eq.zero)then
                        If(f1234r.lt.zero)then
                            If(p.lt.one/rhorizon-one/robs)then
                                radius=robs/(robs*p+one)
                            else
                                radius=rhorizon
                            endif
                        else
                            If(p.lt.one/robs)then
                                radius=robs/(one-robs*p)
                            else
                                radius=infinity
                            endif
                        endif
                    endif
                    If(cc.eq.-27.D0)then
                        sqt3=sqrt(three)
                        If(f1234r.lt.zero)then
                            cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                            dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                            If(p.ne.zero)then
                                radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                            else
                                radius=robs!infinity
                            endif
                        else
                            cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                            dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                            PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                            -Log(one+two/sqt3)/three/sqt3
                            If(p.lt.PI0)then
                                radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                            else
                                radius=infinity
                            endif
                        endif
                    endif
                endif
            else
                u=real(bb(4))
                w=abs(aimag(bb(4)))
                v=abs(aimag(bb(2)))
                If(u.ne.zero)then
                    ! equation (45) in Yang & Wang (2012).
                    L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                    L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                    ! equation (46) in Yang & Wang (2012).
                    thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
                    ! equation (48) in Yang & Wang (2012).
                    m2=(L1-L2)/L1
                    tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                    t_inf=sqrt((L1-one)/(L1-L2))
                    ! equation (50) in Yang & Wang (2012).
                    pinf=EllipticF(tinf,m2)/w/sqrt(L1)
                    call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                    If(f1234r.lt.zero)then
                        PI0=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))
                        if(p.lt.PI0)then
                            ! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                            radius=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                        else
                            radius=rhorizon
                        endif
                    else
                        PI0=EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                        if(p.lt.PI0)then
                            ! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                            radius=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                        else
                            radius=infinity
                        endif
                    endif
                else
                    If(f1234r.lt.zero)then
                        if(p.lt.(atan(robs/w)-atan(rhorizon/w))/w)then
                            radius=w*tan(atan(robs/w)-p*w)
                        else
                            radius=rhorizon
                        endif
                    else
                        if(p.lt.(PI/two-atan(robs/w))/w)then
                            radius=w*tan(atan(robs/w)+p*w)
                        else
                            radius=infinity
                        endif
                    endif
                endif
            endif
            count_num=count_num+1
        else
            If(f1234r.eq.f1234r_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
                a_spin.eq.a_spin_1.and.robs.eq.robs_1.and.scal.eq.scal_1)then
                !***************************************************************************************************
                If(reals.ne.0)then
                    index_p4(1)=0
                    cases_int=1
                    select case(cases)
                    case(1)
                        If(.not.indrhorizon)then
                            If(f1234r.lt.zero)then
                                If(p.lt.PI0_total)then
                                    ! equation (41) in Yang & Wang (2012).
                                    radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)
                                else
                                    radius=infinity  !Goto infinity, far away.
                                endif
                            else
                                If(p.lt.PI0_inf_obs)then
                                    ! equation (41) in Yang & Wang (2012).
                                    radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)
                                else
                                    radius=infinity !Goto infinity, far away.
                                endif
                            endif
                        else
                            If(f1234r.lt.zero)then
                                If(p.lt.PI0_obs_hori)then
                                    ! equation (41) in Yang & Wang (2012).
                                    radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)
                                else
                                    radius=rhorizon !Fall into black hole.
                                endif
                            else
                                If(p.lt.PI0_inf_obs)then
                                    ! equation (41) in Yang & Wang (2012).
                                    radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)
                                else
                                    radius=infinity !Goto infinity, far away.
                                endif
                            endif
                        endif
                    case(2)
                        If(.not.indrhorizon)then
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0
                            endif
                            ! equation (41) in Yang & Wang (2012).
                            radius=r_tp1+b0/(four*weierstrassP(p+PI01,g2,g3,dd,del)-b1)
                        else
                            If(f1234r.le.zero)then
                                If(p.lt.PI0_obs_hori)then
                                    ! equation (41) in Yang & Wang (2012).
                                    radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)
                                else
                                    radius=rhorizon !Fall into black hole.
                                endif
                            else
                                If(p.lt.PI0_total_2)then
                                    ! equation (41) in Yang & Wang (2012).
                                    radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)
                                else
                                    radius=rhorizon !Fall into black hole.
                                endif
                            endif
                        endif
                    end select
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                If(p.lt.one/rhorizon-one/robs)then
                                    radius=robs/(robs*p+one)
                                else
                                    radius=rhorizon
                                endif
                            else
                                If(p.lt.one/robs)then
                                    radius=robs/(one-robs*p)
                                else
                                    radius=infinity
                                endif
                            endif
                        endif
                        If(cc.eq.-27.D0)then
                            sqt3=sqrt(three)
                            If(f1234r.lt.zero)then
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                                If(p.ne.zero)then
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=robs!infinity
                                endif
                            else
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                                PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                -Log(one+two/sqt3)/three/sqt3
                                If(p.lt.PI0)then
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=infinity
                                endif
                            endif
                        endif
                    endif
                else
                    If(u.ne.zero)then
                        call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                        If(f1234r.lt.zero)then
                            if(p.lt.PI0)then
                                ! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                                radius=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                            else
                                radius=rhorizon
                            endif
                        else
                            if(p.lt.PI0)then
                                ! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                                radius=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                            else
                                radius=infinity
                            endif
                        endif
                    else
                        If(f1234r.lt.zero)then
                            if(p.lt.(atan(robs/w)-atan(rhorizon/w))/w)then
                                radius=w*tan(atan(robs/w)-p*w)
                            else
                                radius=rhorizon
                            endif
                        else
                            if(p.lt.(PI/two-atan(robs/w))/w)then
                                radius=w*tan(atan(robs/w)+p*w)
                            else
                                radius=infinity
                            endif
                        endif
                    endif
                endif
            !***************************************************************************************************
            else
                count_num=1
                goto  20
            endif
        endif
        return
    End function radius

    !********************************************************************************************
    !   Function phi(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)
    Function phi(p,f1234,lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,scal)
        !********************************************************************************************
        !*     PURPOSE:  Computes function \phi(p).
        !*     INPUTS:   p--------------independent variable, which must be nonnegative.
        !*               f1234----------array of p_1, p_2, p_3, p_4, which are the components of four-
        !*                              momentum of a photon measured under the LNRF frame. This array
        !*                              can be computed by subroutine lambdaq(...), see below.
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  phi------------value of function \phi(p).
        !*     ROUTINES CALLED: INTRPART, INTTPART.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        Implicit none
        Double precision p,sinobs,muobs,a_spin,lambdaBar,phi,phi_r,phi_c,twopi,f1234(4),timet,&
        two,Rab,robs,scal,zero,one,lambda,q,&
        time_r,aff_r,mu_cos,r_coord
        parameter (zero=0.D0,one=1.D0,two=2.D0)
        logical :: rotate,err
        integer  tm1,tm2,tr1,tr2

        twopi=two*PI
        rotate=.false.
        err=.false.
        IF(a_spin.EQ.ZERO.and.(f1234(2).NE.zero.or.muobs.NE.zero))THEN
            !When spin a=0, and the motion of photon is not confined in the equatorial plane, then \phi_r = 0.
            phi_r=zero
        ELSE
            ! call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
            !                     scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
            call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,&
            scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
        ENDIF
        ! Call integrat_theta_part to evaluate \phi_\mu.
        ! call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal,phi_c,timet,mu_cos,tm1,tm2)
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,lambdaBar,scal,phi_c,timet,mu_cos,tm1,tm2)
        !write(*,*)'ss=', p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,scal
        !write(*,*)'phi2=',phi_r,phi_c,tm1,tm2
        If(abs(muobs).ne.one)then
            ! When the observer is on the equatorial plane, and p_\theta (f1234(2)) = 0, then the photon is
            ! confined in the equatorial plane.
            If(muobs.eq.zero.and.f1234(2).eq.zero)then
                phi_c=zero
            endif
            !Equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_c)
            If(f1234(3).eq.zero)then
                phi=phi+(tm1+tm2)*PI
            endif
            phi=dMod(phi,twopi)
            If(phi.lt.zero)then
                phi=phi+twopi
            Endif
        else
            !Equation (74) of Yang & Wang (2012).
            phi=-(phi_c+phi_r+(tm1+tm2)*PI)
            Rab=sqrt(f1234(3)**two+f1234(2)**two)
            If(phi.ne.zero)then
                rotate=.true.
            endif
            If(Rab.ne.zero)then
                ! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=phi+dasin(f1234(2)/Rab)
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=phi+twopi+dasin(f1234(2)/Rab)
                endif
            else
                phi=zero
            endif
            If(rotate)then
                phi=Mod(phi,twopi)
                If(phi.lt.zero)then
                    phi=phi+twopi
                Endif
            Endif
        endif
        return
    End Function phi

    !********************************************************************************************
    ! !  SUBROUTINE GEOKERR(p_int,rp,mup,varble,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
    ! !                     tr1,tr2,tm1,tm2,radi,mu,time,phi,sigma)
    ! SUBROUTINE GEOKERR(p_int,rp,mup,varble,f1234,lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,scal,&
    !     tr1,tr2,tm1,tm2,radi,mu,time,phi,sigma)
    !     !********************************************************************************************
    !     !*     PURPOSE:  Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine parameter
    !     !*               \sigma as functions of parameter p, i.e. functions r(p), \mu(p), \phi(p), t(p)
    !     !*               and \sigma(p). Cf. discussions in Yang & Wang (2012).
    !     !*     INPUTS:   p_int----------this parameter will be taken as independent variable, if
    !     !*                              varble='p', which must be nonnegative.
    !     !*               rp-------------this parameter will be taken as independent variable, if
    !     !*                              varble='r'.
    !     !*               mup------------this parameter will be taken as independent variable, if
    !     !*                              varble='mu'.
    !     !*               varble---------Tell the routine which parameter to be as independent variable,
    !     !*                              r, mu or p.
    !     !*               f12342---------array of f_1, f_2, f_3, f_4, which was defined by equation (102)-(105)
    !     !*                              in Yang & Wang (2012).
    !     !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
    !     !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
    !     !*                              \theta_{obs} is the inclination angle of the observer.
    !     !*               a_spin---------spin of black hole, on interval (-1,1).
    !     !*               robs-----------radial coordinate of observer or initialposition of photon.
    !     !*               scal-----------a dimentionless parameter to control the size of the images.
    !     !*                              Which is usually be set to 1.D0.
    !     !*               tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
    !     !*                              respectively. If varble='mu', these two parameter must be provided.
    !     !*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
    !     !*                              respectively. If varble='r', these two parameter must be provided.
    !     !*     OUTPUTS:  radi-----------value of function r(p).
    !     !*               mu-------------value of function \mu(p).
    !     !*               phi------------value of function \phi(p).
    !     !*               time-----------value of function t(p).
    !     !*               sigma----------value of function \sigma(p).
    !     !*               tm1,tm2--------number of times of the photon meets turning points \mu_tp1 and \mu_tp2
    !     !*                              respectively.
    !     !*               tr1,tr2--------number of times of the photon meets turning points r_tp1 and r_tp2
    !     !*                              respectively.
    !     !*     ROUTINES CALLED: INTRPART, INTTPART.
    !     !*     ACCURACY:   Machine.
    !     !*     AUTHOR:     Yang & Wang (2012)
    !     !*     DATE WRITTEN:  5 Jan 2012
    !     !*     REVISIONS: ******************************************
    !     IMPLICIT NONE
    !     DOUBLE PRECISION f1234(4),lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,scal,radi,mu,time,phi,sigma,&
    !     zero,one,two,three,four,phi_r,time_r,aff_r,phi_t,time_t,p,Rab,&
    !     rp,mup,p_int,mu_cos,r_coord
    !     CHARACTER varble
    !     PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0)
    !     LOGICAL rotate,err
    !     INTEGER tm1,tm2,tr1,tr2

    !     SELECT CASE(varble)
    !     CASE('r')
    !         radi=rp
    !         ! p=r2p(f1234(1),rp,lambda,q,a_spin,robs,scal,tr1,tr2)
    !         p=r2p(f1234(1),rp,lambda,q,a_spin,lambdaBar,robs,scal,tr1,tr2)
    !         mu=mucos(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)
    !     CASE('mu')
    !         mu=mup
    !         p=mu2p(f1234(3),f1234(2),lambda,q,mup,sinobs,muobs,a_spin,tm1,tm2,scal)
    !         ! p=mu2p(f1234(3),f1234(2),lambda,q,mup,sinobs,muobs,a_spin,lambdaBar,tm1,tm2,scal)
    !         ! radi=radius(p,f1234(1),lambda,q,a_spin,robs,scal)
    !         radi=radius(p,f1234(1),lambda,q,a_spin,lambdaBar,robs,scal)
    !     CASE('p')
    !         p=p_int
    !     END SELECT

    !     !************************************************************************************
    !     ! Call integrate_r_part to evaluate t_r,\phi_r,\sigma_r, and function r(p)=r_coord.
    !     ! call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
    !     !                       scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
    !     call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,&
    !     scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
    !     ! Call integrate_theta_part to evaluate t_\mu,\phi_\mu,\sigma_\mu, and function \mu(p)=mu_cos.
    !     ! call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal,&
    !     !                       phi_t,time_t,mu_cos,tm1,tm2)
    !     call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,lambdaBar,scal,&
    !     phi_t,time_t,mu_cos,tm1,tm2)
    !     radi=r_coord
    !     mu=mu_cos
    !     !time coordinate value **************************************************************
    !     time=time_r+time_t
    !     !affine parameter value *************************************************************
    !     sigma=aff_r+time_t
    !     !phi coordinate value ***************************************************************
    !     rotate=.false.
    !     err=.false.
    !     IF(ABS(muobs).NE.ONE)THEN
    !         ! Equation (74) of Yang & Wang (2012).
    !         phi=-(phi_r+phi_t)
    !         IF(f1234(3).EQ.zero)THEN
    !             phi=phi+(tm1+tm2)*PI
    !         ENDIF
    !         phi=DMOD(phi,twopi)
    !         IF(phi.LT.zero)THEN
    !             phi=phi+twopi
    !         ENDIF
    !     ELSE
    !         ! Equation (74) of Yang & Wang (2012).
    !         phi=-(phi_t+phi_r+(tm1+tm2)*PI)
    !         Rab=dsqrt(f1234(3)**two+f1234(2)**two)
    !         IF(phi.NE.zero)THEN
    !             rotate=.TRUE.
    !         ENDIF
    !         IF(Rab.NE.zero)THEN
    !             ! a muobs was multiplied to control the rotate direction
    !             if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
    !                 phi=muobs*phi+asin(f1234(2)/Rab)
    !             endif
    !             if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
    !                 phi=muobs*phi+PI-asin(f1234(2)/Rab)
    !             endif
    !             if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
    !                 phi=muobs*phi+PI-asin(f1234(2)/Rab)
    !             endif
    !             if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
    !                 phi=muobs*phi+twopi+asin(f1234(2)/Rab)
    !             endif
    !         ELSE
    !             phi=zero
    !         ENDIF
    !         IF(rotate)THEN
    !             phi=Mod(phi,twopi)
    !             IF(phi.LT.zero)THEN
    !                 phi=phi+twopi
    !             ENDIF
    !         ENDIF
    !     ENDIF
    !     RETURN
    ! END SUBROUTINE GEOKERR

    !**************************************************
    Function rms(a_spin)
        !**************************************************
        !*     PURPOSE: Computes inner most stable circular orbit r_{ms}.
        !*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
        !*     OUTPUTS:  radius of inner most stable circular orbit: r_{ms}
        !*     ROUTINES CALLED: root4
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: *******************************************************
        implicit none
        Double precision rms,a_spin,b,c,d,e
        complex*16 rt(1:4)
        integer  reals,i
        If(a_spin.eq.0.D0)then
            rms=6.D0
            return
        endif
        b=0.D0
        c=-6.D0
        d=8.D0*a_spin
        e=-3.D0*a_spin**2
        ! Bardeen et al. (1972)
        call root4(b,c,d,e,rt(1),rt(2),rt(3),rt(4),reals)
        Do i=4,1,-1
            If(aimag(rt(i)).eq.0.D0)then
                rms=real(rt(i))**2
                return
            endif
        enddo
    end function rms

    !***********************************************************
    Function rph(a_spin)
        !***********************************************************
        !*     PURPOSE: Computes photon orbit of circluar orbits: r_{ph}.
        !*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
        !*     OUTPUTS:  radius of photon orbit: r_{ph}
        !*     ROUTINES CALLED: NONE
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision rph,a_spin
        ! Bardeen et al. (1972)
        rph=2.D0*(1.D0+cos(2.D0/3.D0*acos(-a_spin)))
    End function  rph

    !*************************************************************
    Function rmb(a_spin)
        !*************************************************************
        !*     PURPOSE: Computes marginally bound orbit of circluar orbits: r_{mb}.
        !*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
        !*     OUTPUTS:  radius of marginally bound orbit: r_{mb}
        !*     ROUTINES CALLED: NONE
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision rmb,a_spin
        ! Bardeen et al. (1972)
        rmb=2.D0-a_spin+2.D0*sqrt(1.D0-a_spin)
    End function  rmb

    !********************************************************************************************
    subroutine mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)
        !********************************************************************************************
        !*     PURPOSE: Returns the coordinates of turning points \mu_tp1 and \mu_tp2 of poloidal motion, judges
        !*                whether the initial poloidal angle \theta_{obs} is one of turning points, if
        !*                it is true then mobseqmtp=.TRUE..
        !*     INPUTS:   f12342---------p_2, the \theta component of four momentum of the photon measured
        !*                              under the LNRF, see equation (84) in Yang & Wang (2012).
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*     OUTPUTS:  mu_tp1, mu_tp2----the turning points, between which the poloidal motion of
        !*                                 the photon was confined, and mu_tp2 <= mu_tp1.
        !*               reals------number of real roots of equation \Theta_\mu(\mu)=0.
        !*               mobseqmtp---If mobseqmtp=.TRUE., then muobs equals to be one of the turning points.
        !*     ROUTINES CALLED: NONE
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision sinobs,muobs,a_spin,lambda,q,zero,one,two,four,&
        mu_tp1,mu_tp2,delta,mutemp,f12342,f12343
        integer  reals
        logical :: mobseqmtp
        parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0)

        mobseqmtp=.false.
        If(a_spin .eq. zero)then
            If(f12342.ne.zero)then
                mu_tp1=sqrt(q/(lambda**two+q))
                mu_tp2=-mu_tp1
            else
                mu_tp1=abs(muobs)
                mu_tp2=-mu_tp1
                mobseqmtp=.true.
            endif
            reals=2
        ELSE
            If(lambda.ne.zero)then
                delta=(a_spin**two-lambda**two-q)**two+four*a_spin**two*q
                mu_tp1=dsqrt( dabs((dsqrt(delta)-(lambda**two+q-a_spin**two))/two) )/dabs(a_spin)
                If(dsqrt(delta)+(lambda**two+q-a_spin**two).le.zero)then
                    mu_tp2=dsqrt(-(dsqrt(delta)+(lambda**two+q-a_spin**two))/two)/dabs(a_spin)
                    If(f12342.eq.zero)then
                        If(abs(muobs-mu_tp1).le.1.D-4)then
                            mu_tp1=dabs(muobs)
                        else
                            mu_tp2=dabs(muobs)
                        endif
                        mobseqmtp=.true.
                    endif
                    reals=4
                else
                    If(f12342.ne.zero)then
                        mu_tp2=-mu_tp1
                    else
                        mu_tp1=dabs(muobs)
                        mu_tp2=-mu_tp1
                        mobseqmtp=.true.
                    endif
                    reals=2
                endif
            else
                If(abs(muobs).ne.one)then
                    If(q.le.zero)then
                        If(f12342.ne.zero)then
                            mu_tp2=dsqrt(-q)/dabs(a_spin)
                        else
                            mu_tp2=dabs(muobs)!a=B=zero.
                            mobseqmtp=.true.
                        endif
                        reals=4
                    else
                        mu_tp2=-one
                        reals=2
                    endif
                    mu_tp1=one
                else
                    mu_tp1=one
                    If(q.le.zero.and.f12342*f12342+f12343*f12343.ne.zero)then
                        mu_tp2=dsqrt(-q)/dabs(a_spin)
                        reals=4
                    else
                        mu_tp2=-one
                        reals=2
                    endif
                endif
            endif
        ENDIF
        If(abs(muobs).eq.one)mobseqmtp=.true.
        If(muobs.lt.zero.and.reals.eq.4)then
            mutemp=mu_tp1
            mu_tp1=-mu_tp2
            mu_tp2=-mutemp
        endif
        return
    end subroutine mutp

    !============================================================================================
    !   Subroutine radiustp(f12341,a_spin,robs,lambda,q,r_tp1,&
    !                     r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
    Subroutine radiustp(f12341,a_spin,lambdaBar,robs,lambda,q,r_tp1,&
        r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
        !********************************************************************************************
        !*     PURPOSE: Returns the coordinates of turning points r_tp1 and r_tp2 of radial motion, judges
        !*                whether the initial radius robs is one of turning points, if
        !*                it is true then robs_eq_rtp=.TRUE.. And if r_tp1 less or equal r_horizon,
        !*                then indrhorizon=.TRUE. Where r_horizon is the radius of the event horizon.
        !*     INPUTS:   f12341---------p_r, the r component of four momentum of the photon measured
        !*                              under the LNRF, see equation (83) in Yang & Wang (2012).
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of the observer.
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*     OUTPUTS:  r_tp1, r_tp2----the turning points, between which the radial motion of
        !*                                 the photon was confined, and r_tp2 >= r_tp1.
        !*               bb(1:4)----roots of equation R(r)=0.
        !*               reals------number of real roots of equation R(r)=0.
        !*               robs_eq_rtp---If robs_eq_rtp=.TRUE., then robs equal to be one of turning points.
        !*               cases-------If r_tp2=infinity, then cases=1, else cases=2.
        !*               indrhorizon----if r_tp1 less or equals r_horizon, indrhorizon=.TRUE..
        !*     ROUTINES CALLED: root4
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision f12341,a_spin,lambdaBar,robs,lambda,q,r_tp1,r_tp2,&
        zero,one,two,four,b1,c1,d1,e1,r1(2),rthroat,rhorizon
        integer  reals,i,j,cases
        logical :: robs_eq_rtp,indrhorizon
        complex*16 bb(4)
        parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0)

        !   rhorizon=one+sqrt(one-a_spin**two)
        rthroat = one + lambdaBar**two + sqrt((one + lambdaBar**two)**two - a_spin**two)
        rhorizon = rthroat

        robs_eq_rtp=.false.
        indrhorizon=.false.
        b1=zero
        c1=a_spin**two-lambda**two-q
        ! Only-g_rr-modified: use pure Kerr quartic R_K(r)=0 for turning points.
        ! The wormhole correction enters only through the radial integral speed factor
        ! sqrt(DeltaK/DeltaWH), NOT through the quartic.
        d1=two*(q+(lambda-a_spin)**two)
        e1=-q*a_spin**two
        call root4(b1,c1,d1,e1,bb(1),bb(2),bb(3),bb(4),reals)
        SELECT CASE(reals)
        CASE(4)
            IF(f12341.eq.zero)THEN
                IF(dabs(robs-real(bb(4))) .LE. 1.D-4)THEN
                    r_tp1=robs
                    r_tp2=infinity
                    cases=1
                ENDIF
                IF(dabs(robs-real(bb(2))) .LE. 1.D-4)THEN
                    r_tp1=robs
                    r_tp2=real(bb(3))
                    cases=2
                ENDIF
                IF(dabs(robs-real(bb(3))) .LE. 1.D-4)THEN
                    r_tp1=real(bb(2))
                    r_tp2=robs
                    cases=2
                ENDIF
                IF(dabs(robs-real(bb(1))) .LE. 1.D-4)THEN
                    r_tp1=-infinity
                    r_tp2=robs
                    cases=3
                    write(*,*)'radiustp(): wrong! 4 roots, cases = 3'
                    stop
                ENDIF
                robs_eq_rtp = .TRUE.
            ELSE
                If( robs.ge.real(bb(4)) )then
                    r_tp1=real(bb(4))
                    r_tp2=infinity
                    cases=1
                else
                    If( (robs.ge.real(bb(2)) .and. robs.le.real(bb(3))) )then
                        r_tp1=real(bb(2))
                        r_tp2=real(bb(3))
                        cases=2
                    else
                        IF( real(bb(1)) .GT. rhorizon .AND.  robs .LE. real(bb(1))  )THEN
                            write(*,*)'radiustp(): wrong! 4 roots,cases = 3'
                            stop
                            r_tp2=real(bb(1))
                            r_tp1=-infinity
                        ELSE
                            write(*,*)'radiustp(): wrong! 4 roots',robs,bb
                            stop
                        ENDIF
                    endif
                endif
            ENDIF
        CASE(2)
            j=1
            Do  i=1,4
                If (aimag(bb(i)).eq.zero) then
                    r1(j)=real(bb(i))
                    j=j+1
                endif
            Enddo
            If( robs.ge.r1(2) )then
                r_tp1=r1(2)
                r_tp2=infinity
                cases=1
            else
                If( r1(1).ge.rhorizon .and. robs.le.r1(1) )then
                    write(*,*)'radiustp(): wrong! 2 roots, cases = 3'
                    stop
                endif
            endif
            IF(f12341.eq.zero)THEN
                IF(dabs(robs-r1(2)) .LE. 1.D-4)THEN
                    r_tp1=robs
                    r_tp2=infinity
                ENDIF
                IF(dabs(robs-r1(1)) .LE. 1.D-4)THEN
                    r_tp1=-infinity
                    r_tp2=robs
                    write(*,*)'radiustp(): wrong! 2 roots, cases = 3'
                    stop
                ENDIF
                robs_eq_rtp=.TRUE.
            ENDIF
        CASE(0)
            r_tp1=zero
            r_tp2=infinity
            cases=1
        END SELECT

        IF(rhorizon.ge.r_tp1 .and. rhorizon.le.r_tp2)then
            indrhorizon=.true.
        Endif
    End Subroutine radiustp

    !********************************************************************************************
    Function mu2p(f12343,f12342,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal)
        !********************************************************************************************
        !*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute
        !*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
        !*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
        !*               component of 4 momentum of photon.
        !*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon
        !*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
        !*               f12343---------p_\phi, which is the \phi component of four momentum of a photon
        !*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
        !*                              respectively.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*               mu-------------\mu coordinate of photon.
        !*     OUTPUTS:  value of \mu part of integral of (24).
        !*     ROUTINES CALLED: mu2p_schwartz, mutp, root3, weierstrass_int_J3
        !*     ACCURACY:   Machine.
        !*     DATE WRITTEN:  4 Jan 2012
        !*     AUTHOR:     Yang & Wang (2012)
        !*     REVISIONS: ******************************************
        implicit none
        Double precision mu2p,f12342,f12343,mu,sinobs,muobs,a_spin,lambda,q,mu_tp,tposition,tp2,four,&
        b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,a4,b4,two,mu_tp2,&
        scal,zero,one,integ4(4),three,rff_p
        parameter (zero=0.D0,two=2.D0,four=4.D0,one=1.D0,three=3.D0)
        integer  t1,t2,reals,p4,index_p4(4),del,cases
        complex*16 dd(3)
        logical :: mobseqmtp

        If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
            mu2p=zero!-one
            return
        endif
        If(a_spin.eq.zero)then
            call mu2p_schwartz(f12343,f12342,lambda,q,mu,sinobs,muobs,t1,t2,mu2p,scal)
            return
        endif

        a4=zero
        b4=one
        p4=0
        mobseqmtp=.false.
        call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp)
        ! equatorial plane motion.
        If(mu_tp.eq.zero)then
            mu2p=zero
            return
        endif
        ! equations (26)-(29) in Yang & Wang (2012).
        b0=-four*a_spin**2*mu_tp**3+two*mu_tp*(a_spin**2-lambda**2-q)
        b1=-two*a_spin**2*mu_tp**2+one/three*(a_spin**2-lambda**2-q)
        b2=-four/three*a_spin**2*mu_tp
        b3=-a_spin**2
        g2=three/four*(b1**2-b0*b2)
        g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)
        ! equation (30) in Yang & Wang (2012).
        If(abs(mu-mu_tp).ne.zero)then
            tposition=b0/(four*(mu-mu_tp))+b1/four
        else
            tposition=infinity
        endif
        If(muobs.ne.mu_tp)then
            tinf=b0/four/(muobs-mu_tp)+b1/four
        else
            tinf=infinity
        endif

        call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
        index_p4(1)=0
        cases=1

        If(mu.gt.mu_tp.or.mu.lt.mu_tp2)then
            mu2p=-one
            return
        endif
        ! equation (30) in Yang & Wang (2012).
        tp2=b0/four/(mu_tp2-mu_tp)+b1/four
        If(t1.eq.0)then
            p1=zero
        else
            ! equation (53) in Yang & Wang (2012).
            call weierstrass_int_J3(tposition,infinity,dd,del,a4,b4,index_p4,rff_p,integ4,cases)
            p1=integ4(1)
        endif
        If(t2.eq.0)then
            p2=zero
        else
            call weierstrass_int_J3(tp2,tposition,dd,del,a4,b4,index_p4,rff_p,integ4,cases)
            p2=integ4(1)
        endif
        call weierstrass_int_J3(tinf,tposition,dd,del,a4,b4,index_p4,rff_p,integ4,cases)
        pp=integ4(1)

        ! equation (54) in Yang & Wang (2012).
        If(mobseqmtp)then
            If(muobs.eq.mu_tp)then
                mu2p=-pp+two*(t1*p1+t2*p2)
            else
                mu2p=pp+two*(t1*p1+t2*p2)
            endif
        else
            If(f12342.lt.zero)then
                mu2p=pp+two*(t1*p1+t2*p2)
            endif
            If(f12342.gt.zero)then
                mu2p=-pp+two*(t1*p1+t2*p2)
            endif
        endif
        return
    end Function mu2p

    !============================================================================================
    subroutine mu2p_schwartz(f12343,f12342,lambda,q,mu,sinobs,muobs,t1,t2,mu2p,scal)
        !********************************************************************************************
        !*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute
        !*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
        !*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
        !*               component of 4 momentum of photon.
        !*               And black hole spin is zero.
        !*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon
        !*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
        !*               f12343---------p_\phi, which is the \phi component of four momentum of a photon
        !*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
        !*                              respectively.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*               mu-------------\mu coordinate of photon.
        !*     OUTPUTS:  mu2p-----------value of \mu part of integral of (24).
        !*     ROUTINES CALLED: NONE.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision f12343,f12342,mu,sinobs,muobs,mu2p,pp,p1,p2,BB,two,&
        lambda,q,scal,zero,one,mu_tp,mu_tp2
        integer  t1,t2
        parameter(two=2.D0,zero=0.D0,one=1.D0)
        logical :: mobseqmtp
        If(f12343.eq.zero.and.f12342.eq.zero)then !this Theta=q(1-mu^2),so if B=0,then q=0.
            mu2p=-two            !so Theta_mu=0 for ever.But we do not need to
            return  !consider it,for q=0,so the next part means that it will return
        endif !zero value.
        mobseqmtp=.false.
        If(q.gt.zero)then
            BB=sqrt(q)
            If(f12342.ne.zero)then
                mu_tp=sqrt(q/(lambda**two+q))
                mu_tp2=-mu_tp
            else
                mu_tp=muobs
                mu_tp2=-mu_tp
                mobseqmtp=.true.
            endif
            If(abs(muobs).eq.one)mobseqmtp=.true.
            pp=(asin(mu/mu_tp)-asin(muobs/mu_tp))*mu_tp/BB
            If(t1.eq.0)then
                p1=zero
            else
                p1=(PI/two-asin(mu/mu_tp))*mu_tp/BB
            endif
            If(t2.eq.0)then
                p2=zero
            else
                p2=(asin(mu/mu_tp)+PI/two)*mu_tp/BB
            endif
            If(mobseqmtp)then
                If(muobs.eq.mu_tp)then
                    mu2p=-pp+two*(t1*p1+t2*p2)
                else
                    mu2p=pp+two*(t1*p1+t2*p2)
                endif
            else
                mu2p=sign(one,-f12342)*pp+two*(t1*p1+t2*p2)
            endif
        else
            mu2p=zero
        endif
        return
    end subroutine mu2p_schwartz
    !********************************************************************************************
    !   Function r2p(f1234r,rend,lambda,q,a_spin,robs,scal,t1,t2)
    Function r2p(f1234r,rend,lambda,q,a_spin,lambdaBar,robs,scal,t1,t2)
        !============================================================================================
        !*     PURPOSE:  Computes the value of parameter p from radial coordinate. In other words, to compute
        !*               the r part of integral of equation (24), using formula (58) in Yang & Wang (2012).
        !*               (58) is: p=-sign(p_r)*p_0+2*t1*p_2+2*t2*p_2. where p_r is initial radial
        !*               component of 4 momentum of photon.
        !*     INPUTS:   f1234r---------f_1, which was defined by equation (106) in Yang & Wang (2012).
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of the observer.
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*               t1,t2----------Number of photon meets the turning points r_tp1 and r_tp2
        !*                              respectively in radial motion.
        !*     OUTPUTS:  r2p------------value of r part of integral (24) in Yang & Wang (2012).
        !*     ROUTINES CALLED: radiustp, root3, weierstrass_int_j3, EllipticF.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  4 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision r2p,a_spin,lambdaBar,rhorizon,rthroat,q,lambda,scal,zero,integ4(4),&
        cc,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,robs,integ04(4),&
        u,v,w,L1,L2,thorizon,m2,pinf,a4,b4,one,two,four,sqrt3,&
        integ14(4),three,six,nine,r_tp1,r_tp2,f1234r,tp2,tp,t_inf,&
        pp,p1,p2,rend,rff_p
        parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
        complex*16 bb(1:4),dd(3)
        integer  reals,cases_int,del,index_p4(4),cases,t1,t2
        logical :: robs_eq_rtp,indrhorizon

        !   rhorizon=one+sqrt(one-a_spin**2)
        rthroat = (one + lambdaBar**two) + sqrt((one + lambdaBar**two)**two - a_spin**two)
        rhorizon = rthroat
        a4=zero
        b4=one
        cc=a_spin**2-lambda**2-q
        robs_eq_rtp=.false.
        indrhorizon=.false.
        call radiustp(f1234r,a_spin,lambdaBar,robs,lambda,q,r_tp1,&
        r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
        If(reals.ne.0)then
            If(rend.lt.r_tp1.or.rend.gt.r_tp2)then
                r2p=-one
                return
            endif
            ! equations (35)-(38) in Yang & Wang (2012). Pure Kerr quartic coefficients.
            b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
            b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
            b2=four/three*r_tp1
            b3=one
            g2=three/four*(b1**2-b0*b2)
            g3=one/16.D0*(3*b0*b1*b2-2*b1**3-b0**2*b3)
            ! equation (39) in Yang & Wang (2012).
            If(robs-r_tp1.ne.zero)then
                tinf=b0/four/(robs-r_tp1)+b1/four
            else
                tinf=infinity
            endif
            If(rhorizon-r_tp1.ne.zero)then
                thorizon=b1/four+b0/four/(rhorizon-r_tp1)
            else
                thorizon=infinity
            endif
            If(rend-r_tp1.ne.zero)then
                tp=b1/four+b0/four/(rend-r_tp1)
            else
                tp=infinity
            endif
            tp2=b0/four/(r_tp2-r_tp1)+b1/four
            tinf1=b1/four

            call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
            index_p4(1)=0
            cases_int=1
            ! equation (42) in Yang & Wang (2012).
            call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
            PI0=integ04(1)
            select case(cases)
            case(1)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then
                        call weierstrass_int_j3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                        pp=integ14(1)
                        If(t1.ne.zero)then
                            ! equation (57) in Yang & Wang (2012).
                            call weierstrass_int_j3(tp,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            p1=integ14(1)
                        else
                            p1=zero  !Goto infinity, far away.
                        endif
                        ! equation (58) in Yang & Wang (2012).
                        r2p=pp+two*p1*t1
                    else
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        pp=integ04(1)
                        r2p=-pp
                    endif
                else
                    If(f1234r.lt.zero)then
                        If(rend.le.rhorizon)then
                            tp=thorizon
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            r2p=integ04(1)
                        else
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            r2p=integ04(1)
                        endif
                    else
                        If(rend.lt.infinity)then
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            r2p=-pp
                        else
                            call weierstrass_int_J3(tinf,tinf1,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            r2p=-pp
                        endif
                    endif
                endif
            case(2)
                If(.not.indrhorizon)then
                    ! equation (57) in Yang & Wang (2012).
                    call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                    pp=integ4(1)
                    If(t1.eq.zero)then
                        p1=zero
                    else
                        call weierstrass_int_J3(tp,infinity,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        p1=integ4(1)
                    endif
                    If(t2.eq.zero)then
                        p2=zero
                    else
                        call weierstrass_int_J3(tp2,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        p2=integ4(1)
                    endif
                    If(f1234r.ne.zero)then
                        r2p=sign(one,-f1234r)*pp+two*(t1*p1+t2*p2)
                    else
                        ! equation (58) in Yang & Wang (2012).
                        If(robs.eq.r_tp1)then
                            r2p=-pp+two*(t1*p1+t2*p2)
                        else
                            r2p=pp+two*(t1*p1+t2*p2)
                        endif
                    endif
                else
                    If(f1234r.le.zero)then
                        If(rend.le.rhorizon)then
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            pp=integ4(1)
                        else
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            pp=integ4(1)
                        endif
                    else
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        pp=integ4(1)
                        If(t2.eq.zero)then
                            p2=zero
                        else
                            call weierstrass_int_J3(tp2,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            p2=integ4(1)
                        endif
                        ! equation (58) in Yang & Wang (2012).
                        r2p=-pp+two*t2*p2
                    endif
                endif
            end select
            If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        If(rend.le.rhorizon)then
                            r2p=one/rhorizon-one/robs
                        else
                            r2p=one/rend-one/robs
                        endif
                    else
                        If(rend.lt.infinity)then
                            r2p=one/robs-one/rend
                        else
                            r2p=one/robs
                        endif
                    endif
                endif
                If(cc.eq.-27.D0)then
                    sqrt3=sqrt(three)
                    If(f1234r.lt.zero)then
                        If(rend.gt.rhorizon)then
                            r2p=-log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/(sqrt3))/&
                            (robs-three)))/(three*sqrt3)+&
                            log(abs((sqrt(rend*(rend+6.D0))+(three+two*rend)/&
                            (sqrt3))/(rend-three)))/(three*sqrt3)
                        else
                            r2p=-log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/&
                            (sqrt3))/(robs-three)))/(three*sqrt3)+&
                            log(abs((sqrt(rhorizon*(rhorizon+6.D0))+(three+two*rhorizon)&
                            /(sqrt3))/(rhorizon-three)))/(three*sqrt3)
                        endif
                    else
                        If(rend.lt.infinity)then
                            r2p=-log(abs((sqrt(rend*(rend+6.D0))+(three+two*rend)/(sqrt3))/&
                            (rend-three)))/(three*sqrt3)+log(abs((sqrt(robs*(robs+6.D0))+&
                            (three+two*robs)/(sqrt3))/(robs-three)))/(three*sqrt3)
                        else
                            r2p=-log(one+two/sqrt3)/three/sqrt3+&
                            log(abs((sqrt(rend*(rend+6.D0))+(three+two*rend)/&
                            (sqrt3))/(rend-three)))/(three*sqrt3)
                        endif
                    endif
                endif
            endif
        else
            ! equation (44) in Yang & Wang (2012).
            u=real(bb(4))
            w=abs(aimag(bb(4)))
            v=abs(aimag(bb(2)))
            If(u.ne.zero)then
                ! equation (45) in Yang & Wang (2012).
                L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                ! equation (46) in Yang & Wang (2012).
                thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
                tp=sqrt((L1-one)/(L1-L2))*(rend-u*(L1+one)/(L1-one))/sqrt((rend-u)**2+w**2)
                ! equation (48) in Yang & Wang (2012).
                m2=(L1-L2)/L1
                tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                t_inf=sqrt((L1-one)/(L1-L2))
                ! equation (50) in Yang & Wang (2012).
                pinf=EllipticF(tinf,m2)/w/sqrt(L1)
                If(f1234r.lt.zero)then
                    If(rend.le.rhorizon)then
                        r2p=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))
                    else
                        r2p=pinf-EllipticF(tp,m2)/(w*sqrt(L1))
                    endif
                else
                    If(rend.lt.infinity)then
                        r2p=EllipticF(tp,m2)/(w*sqrt(L1))-pinf
                    else
                        r2p=EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                    endif
                endif
            else
                If(f1234r.lt.zero)then
                    If(rend.le.rhorizon)then
                        r2p=(atan(robs/w)-atan(rhorizon/w))/w
                    else
                        r2p=(atan(robs/w)-atan(rend/w))/w
                    endif
                else
                    if(rend.lt.infinity)then
                        r2p=(atan(rend/w)-atan(robs/w))/w
                    else
                        r2p=(PI/two-atan(robs/w))/w
                    endif
                endif
            endif
        endif
        return
    End function r2p

    !********************************************************************************************
    !   SUBROUTINE INTTPART(p,f12343,f12342,lambda,q,sinobs,muobs,a_spin,scal,phyt,timet,mucos,t1,t2)
    SUBROUTINE INTTPART(p,f12343,f12342,lambda,q,sinobs,muobs,a_spin,lambdaBar,scal,phyt,timet,mucos,t1,t2)
        !********************************************************************************************
        !*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, t and affine parameter \sigma,
        !*               expressed by equation (71) and (72) in Yang & Wang (2012).
        !*     INPUTS:   p--------------independent variable, which must be nonnegative.
        !*               f12342---------p_\theta, which is the \theta component of four momentum of a photon
        !*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
        !*               f12343---------p_\phi, which is the \phi component of four momentum of a photon
        !*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  phyt-----------value of integral \phi_\theta expressed by equation (72) in
        !*                              Yang & Wang (2012).
        !*               timet----------value of integral \t_\theta expressed by equation (71) in
        !*                              Yang & Wang (2012). And \sigma_\theta=time_\theta.
        !*               mucos----------value of function \mu(p).
        !*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
        !*                              respectively.
        !*     ROUTINES CALLED: mutp, root3, weierstrass_int_J3
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  5 Jan 2012
        !*     REVISIONS: ******************************************
        USE constants
        IMPLICIT NONE
        Double precision phyt,timet,f12343,f12342,p,sinobs,muobs,a_spin,lambdaBar,lambda,q,mu_tp1,tp2,tmu,&
        b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,Wmup,Wmum,tplus,tminus,&
        a4,b4,mu_tp2,scal,integ4(4),integ(4),rff_p,&
        integ14(4),PI0,integ04(4),mu2p,PI01,h,p1_t,p2_t,pp_t,p1_phi,&
        p2_phi,pp_phi,mucos,p_mt1_mt2,&
        PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p
        Double precision f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1
        integer ::  t1,t2,i,j,reals,index_p4(4),del,cases_int,count_num=1
        complex*16 dd(3)
        logical :: mobseqmtp
        save  f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1,a4,b4,mu_tp1,mu_tp2,reals,&
        mobseqmtp,b0,b1,b2,b3,g2,g3,dd,del,PI0,Wmup,Wmum,tplus,tminus,tp2,tinf,h,p_mt1_mt2,&
        PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p,PI01

        30      continue
        IF(count_num.eq.1)then
            f12343_1=f12343
            f12342_1=f12342
            lambda_1=lambda
            q_1=q
            muobs_1=muobs
            sinobs_1=sinobs
            a_spin_1=a_spin
            scal_1=scal
            t1=0
            t2=0
            !************************************************************************
            If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
                mucos = sign(one,muobs)
                timet=zero             !this is because that mu==1 for ever
                phyt=zero              !this is because that mu==1 for ever,this
                count_num=count_num+1
                return        !because that Theta_mu=-a^2(1-mu^2), so,mu must =+1 or -1 for ever.
            endif
            If(muobs.eq.zero.and.(abs(lambda).lt.abs(a_spin)).and.q.eq.zero)then
                timet=zero
                phyt=zero
                mucos=zero
                count_num=count_num+1
                return
            endif
            mobseqmtp=.false.
            call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)
            If(mu_tp1.eq.zero)then
                !photons are confined in the equatorial plane, so the integrations about \theta are valished.
                timet=zero
                phyt=zero
                mucos=zero
                count_num=count_num+1
                return
            endif
            !**************************************************************************
            If(a_spin.eq.zero)then
                timet=zero
                CALL phyt_schwatz(p,f12343,f12342,lambda,q,sinobs,muobs,scal,phyt,mucos,t1,t2)
                count_num=count_num+1
                return
            endif
            a4=zero
            b4=one
            ! equations (26)-(29) in Yang & Wang (2012).
            b0=-four*a_spin**2*mu_tp1**3+two*mu_tp1*(a_spin**2-lambda**2-q)
            b1=-two*a_spin**2*mu_tp1**2+one/three*(a_spin**2-lambda**2-q)
            b2=-four/three*a_spin**2*mu_tp1
            b3=-a_spin**2
            g2=three/four*(b1**2-b0*b2)
            g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)
            call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
            ! equation (30) in Yang & Wang (2012).
            If(muobs.ne.mu_tp1)then
                tinf=b0/four/(muobs-mu_tp1)+b1/four
            else
                tinf=infinity
            endif
            If(mu_tp1-one.ne.zero)then
                ! equation (64) in Yang & Wang (2012).
                Wmum=b0/(eight*(-one-mu_tp1)**2)
                Wmup=b0/(eight*(one-mu_tp1)**2)
                tminus=b0/four/(-one-mu_tp1)+b1/four
                tplus=b0/four/(one-mu_tp1)+b1/four
            endif
            index_p4(1)=0
            cases_int=1
            call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
            PI0=integ04(1)
            ! equation (34) in Yang & Wang (2012).
            If(f12342.lt.zero)then
                PI01=-PI0
            else
                PI01=PI0
            endif
            tmu=weierstrassP(p+PI01,g2,g3,dd,del)
            ! equation (32) in Yang & Wang (2012).
            mucos = mu_tp1+b0/(four*tmu-b1)
            tp2=b0/four/(mu_tp2-mu_tp1)+b1/four
            h=-b1/four
            !to get number of turning points of t1 and t2.
            !111111111**********************************************************************************
            call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
            call weierstrass_int_J3(tinf,tmu,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
            !write(*,*)tp2,tinf,tmu,mu_tp1,mu_tp2
            ! equation (51) in Yang & Wang (2012).
            p_mt1_mt2=integ14(1)
            PI2_p=p_mt1_mt2-PI0
            pp=integ4(1)
            p1=PI0-pp
            p2=p_mt1_mt2-p1
            PI1_phi=zero
            PI2_phi=zero
            PI1_time=zero
            PI2_time=zero
            Do j=0,10
                Do i=j,j+1
                    If(mobseqmtp)then
                        If(muobs.eq.mu_tp1)then
                            t1=j
                            t2=i
                            ! equation (52) in Yang & Wang (2012).
                            mu2p=-pp+two*(t1*p1+t2*p2)
                        else
                            t1=i
                            t2=j
                            mu2p=pp+two*(t1*p1+t2*p2)
                        endif
                    else
                        If(f12342.lt.zero)then
                            t1=i
                            t2=j
                            ! equation (52) in Yang & Wang (2012).
                            mu2p=pp+two*(t1*p1+t2*p2)
                        endif
                        If(f12342.gt.zero)then
                            t1=j
                            t2=i
                            ! equation (52) in Yang & Wang (2012).
                            mu2p=-pp+two*(t1*p1+t2*p2)
                        endif
                    endif
                    !write(*,*)p,mu2p,abs(p-mu2p),pp,p1,p2!
                    If(abs(p-mu2p).lt.1.D-4)goto 400
                enddo
            enddo
            !11111111*****************************************************************************************
            400 continue
            index_p4(1)=-1
            index_p4(2)=-2
            index_p4(3)=0
            index_p4(4)=-4
            !*****pp part***************************************
            If(lambda.ne.zero)then
                cases_int=2
                call weierstrass_int_J3(tinf,tmu,dd,del,-tplus,b4,index_p4,abs(pp),integ4,cases_int)
                call weierstrass_int_J3(tinf,tmu,dd,del,-tminus,b4,index_p4,abs(pp),integ14,cases_int)
                ! equation (72) in Yang & Wang (2012).
                pp_phi=lambda*(pp/(one-mu_tp1*mu_tp1)+integ4(2)*Wmup-integ14(2)*Wmum)
            else
                pp_phi=zero
            endif
            cases_int=4
            call weierstrass_int_J3(tinf,tmu,dd,del,h,b4,index_p4,abs(pp),integ,cases_int)
            ! equation (71) in Yang & Wang (2012).
            pp_t=a_spin**two*(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
            !*****p1 part***************************************
            If(t1.eq.0)then
                p1_phi=zero
                p1_t=zero
            else
                If(lambda.ne.zero)then
                    IF(PI1_phi .EQ. zero)THEN
                        cases_int=2
                        call weierstrass_int_J3(tinf,infinity,dd,del,-tplus,b4,index_p4,PI0,integ4,cases_int)
                        call weierstrass_int_J3(tinf,infinity,dd,del,-tminus,b4,index_p4,PI0,integ14,cases_int)
                        ! equation (72) in Yang & Wang (2012).
                        PI1_phi=lambda*(PI0/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum)
                    ENDIF
                    ! equation (51) in Yang & Wang (2012).
                    p1_phi=PI1_phi-pp_phi
                else
                    p1_phi=zero
                endif
                IF(PI1_time .EQ. zero)THEN
                    cases_int=4
                    call weierstrass_int_J3(tinf,infinity,dd,del,h,b4,index_p4,PI0,integ,cases_int)
                    ! equation (62) in Yang & Wang (2012).
                    PI1_time=a_spin**two*(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                ENDIF
                ! equation (51) in Yang & Wang (2012).
                p1_t=PI1_time-pp_t
            endif
            !*****p2 part***************************************
            If(t2.eq.0)then
                p2_phi=zero
                p2_t=zero
            else
                IF(lambda.ne.zero)then
                    IF(PI2_phi .EQ. zero)THEN
                        cases_int=2
                        call weierstrass_int_J3(tp2,tinf,dd,del,-tplus,b4,index_p4,PI2_p,integ4,cases_int)
                        call weierstrass_int_J3(tp2,tinf,dd,del,-tminus,b4,index_p4,PI2_p,integ14,cases_int)
                        ! equation (72) in Yang & Wang (2012).
                        PI2_phi=lambda*(PI2_p/(one-mu_tp1*mu_tp1)+integ4(2)*Wmup-integ14(2)*Wmum)
                    ENDIF
                    ! equation (51) in Yang & Wang (2012).
                    p2_phi=PI2_phi+pp_phi
                ELSE
                    p2_phi=zero
                ENDIF

                IF(PI2_time .EQ. zero)THEN
                    cases_int=4
                    call weierstrass_int_J3(tp2,tinf,dd,del,h,b4,index_p4,PI2_p,integ,cases_int)
                    ! equation (71) in Yang & Wang (2012).
                    PI2_time=a_spin**two*(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                ENDIF
                ! equation (51) in Yang & Wang (2012).
                p2_t=PI2_time+pp_t
            !write(*,*)'ynogk=',tp2,tinf,h,dd
            endif
            !**************************************************************
            ! equation (52) in Yang & Wang (2012).
            !write(*,*)'ynogk=',pp_t,p1_t,p2_t,t1,t2
            If(mobseqmtp)then
                If(muobs.eq.mu_tp1)then
                    phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                else
                    phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=pp_t+two*(t1*p1_t+t2*p2_t)
                endif
            else
                If(f12342.lt.zero)then
                    phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=pp_t+two*(t1*p1_t+t2*p2_t)
                endif
                If(f12342.gt.zero)then
                    phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                endif
            endif
            count_num=count_num+1
        ELSE
            If(f12343_1.eq.f12343.and.f12342_1.eq.f12342.and.lambda_1.eq.lambda.and.q_1.eq.q.and.sinobs_1.eq.sinobs&
                .and.muobs_1.eq.muobs.and.a_spin_1.eq.a_spin.and.scal_1.eq.scal)then
                !***************************************************************************
                t1=0
                t2=0
                If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
                    mucos = sign(one,muobs)
                    timet=zero      !this is because that mu==1 for ever
                    phyt=zero       !this is because that mu==1 for ever,this because that Theta_mu=-a^2(1-mu^2)
                    return          !so,mu must =+1 or -1 for ever.
                endif
                If(muobs.eq.zero.and.(abs(lambda).lt.abs(a_spin)).and.q.eq.zero)then
                    timet=zero
                    phyt=zero
                    mucos=zero
                    return
                endif
                If(mu_tp1.eq.zero)then
                    !photons are confined in the equatorial plane, so the integrations about \theta are valished.
                    timet=zero
                    phyt=zero
                    mucos=zero
                    return
                endif
                If(a_spin.eq.zero)then
                    timet=zero
                    CALL phyt_schwatz(p,f12343,f12342,lambda,q,sinobs,muobs,scal,phyt,mucos,t1,t2)
                    return
                endif

                tmu=weierstrassP(p+PI01,g2,g3,dd,del)
                mucos=mu_tp1+b0/(four*tmu-b1)
                !get numbers of turn points of t1 and t2.
                !111111111************************************************************
                index_p4(1)=0
                cases_int=1
                call weierstrass_int_J3(tinf,tmu,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                pp=integ4(1)
                p1=PI0-pp
                p2=p_mt1_mt2-p1
                !p1=zero
                !p2=zero
                Do j=0,10
                    Do i=j,j+1
                        If(mobseqmtp)then
                            If(muobs.eq.mu_tp1)then
                                t1=j
                                t2=i
                                ! equation (54) in Yang & Wang (2012).
                                mu2p=-pp+two*(t1*p1+t2*p2)
                            else
                                t1=i
                                t2=j
                                mu2p=pp+two*(t1*p1+t2*p2)
                            endif
                        else
                            If(f12342.lt.zero)then
                                t1=i
                                t2=j
                                ! equation (54) in Yang & Wang (2012).
                                mu2p=pp+two*(t1*p1+t2*p2)
                            endif
                            If(f12342.gt.zero)then
                                t1=j
                                t2=i
                                ! equation (54) in Yang & Wang (2012).
                                mu2p=-pp+two*(t1*p1+t2*p2)
                            endif
                        endif
                        !write(*,*)p,mu2p,t1,t2
                        If(abs(p-mu2p).lt.1.D-4)goto 410
                    enddo
                enddo
                !11111111*********************************************
                410 continue
                index_p4(1)=-1
                index_p4(2)=-2
                index_p4(3)=0
                index_p4(4)=-4
                !*****pp parts************************************
                If(lambda.ne.zero)then
                    cases_int=2
                    call weierstrass_int_J3(tinf,tmu,dd,del,-tplus,b4,index_p4,abs(pp),integ4,cases_int)
                    call weierstrass_int_J3(tinf,tmu,dd,del,-tminus,b4,index_p4,abs(pp),integ14,cases_int)
                    ! equation (72) in Yang & Wang (2012).
                    pp_phi=lambda*(pp/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum)
                else
                    pp_phi=zero
                endif
                cases_int=4
                call weierstrass_int_J3(tinf,tmu,dd,del,h,b4,index_p4,abs(pp),integ,cases_int)
                ! equation (71) in Yang & Wang (2012).
                pp_t=a_spin**two*(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                !*****p1 parts************************************
                If(t1.eq.0)then
                    p1_phi=zero
                    p1_t=zero
                else
                    If(lambda.ne.zero)then
                        IF(PI1_phi .EQ. zero)THEN
                            cases_int=2
                            call weierstrass_int_J3(tinf,infinity,dd,del,-tplus,b4,index_p4,PI0,integ4,cases_int)
                            call weierstrass_int_J3(tinf,infinity,dd,del,-tminus,b4,index_p4,PI0,integ14,cases_int)
                            ! equation (72) in Yang & Wang (2012).
                            PI1_phi=lambda*(PI0/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum)
                        ENDIF
                        ! equation (51) in Yang & Wang (2012).
                        p1_phi=PI1_phi-pp_phi
                    else
                        p1_phi=zero
                    endif
                    IF(PI1_time .EQ. zero)THEN
                        cases_int=4
                        call weierstrass_int_J3(tinf,infinity,dd,del,h,b4,index_p4,PI0,integ,cases_int)
                        ! equation (71) in Yang & Wang (2012).
                        PI1_time=a_spin**two*(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                    ENDIF
                    ! equation (51) in Yang & Wang (2012).
                    p1_t=PI1_time-pp_t
                endif
                !*****p2 parts************************************
                If(t2.eq.0)then
                    p2_phi=zero
                    p2_t=zero
                else
                    If(lambda.ne.zero)then
                        IF(PI2_phi .EQ. zero)THEN
                            cases_int=2
                            call weierstrass_int_J3(tp2,tinf,dd,del,-tplus,b4,index_p4,PI2_p,integ4,cases_int)
                            call weierstrass_int_J3(tp2,tinf,dd,del,-tminus,b4,index_p4,PI2_p,integ14,cases_int)
                            ! equation (72) in Yang & Wang (2012).
                            PI2_phi=lambda*(PI2_p/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum)
                        ENDIF
                        ! equation (51) in Yang & Wang (2012).
                        p2_phi=PI2_phi+pp_phi
                    else
                        p2_phi=zero
                    endif
                    IF(PI2_time .EQ. zero)THEN
                        cases_int=4
                        call weierstrass_int_J3(tp2,tinf,dd,del,h,b4,index_p4,PI2_p,integ,cases_int)
                        ! equation (71) in Yang & Wang (2012).
                        PI2_time=a_spin**two*(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                    ENDIF
                    ! equation (51) in Yang & Wang (2012).
                    p2_t=PI2_time+pp_t
                endif
                !**************************************************************
                ! equation (52) in Yang & Wang (2012).
                If(mobseqmtp)then
                    If(muobs.eq.mu_tp1)then
                        phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                    else
                        phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timet=pp_t+two*(t1*p1_t+t2*p2_t)
                    endif
                else
                    If(f12342.lt.zero)then
                        phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timet=pp_t+two*(t1*p1_t+t2*p2_t)
                    endif
                    If(f12342.gt.zero)then
                        phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                    endif
                endif
            else
                count_num=1
                goto 30
            endif
        ENDIF
        !write(*,*)'ff1=',phyt,timet,pp_phi,p1_phi,p2_phi,t1,t2
        RETURN
    END SUBROUTINE INTTPART

    !********************************************************************************************
    SUBROUTINE phyt_schwatz(p,f3,f2,lambda,q,sinobs,muobs,scal,phyc_schwatz,mucos,t1,t2)
        !********************************************************************************************
        !*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, expressed by equation (72)
        !*               in Yang & Wang (2012) with zero spin of black hole.
        !*     INPUTS:   p--------------independent variable, which must be nonnegative.
        !*               f2-------------p_\theta, which is the \theta component of four momentum of a photon
        !*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
        !*               f3-------------p_\phi, which is the \phi component of four momentum of a photon
        !*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  phyc_schwatz-----------value of integral \phi_\theta expressed by equation (71) in
        !*                              Yang & Wang (2012).
        !*               mucos----------value of function \mu(p) with zero spin.
        !*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
        !*                              respectively.
        !*     ROUTINES CALLED: schwatz_int
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  5 Jan 2012
        !*     REVISIONS: ******************************************
        USE constants
        implicit none
        Double precision phyc_schwatz,f3,f2,p,sinobs,muobs,pp,p1,p2,mu,&
        lambda,AA,BB,scal,q,mu_tp1,mu_tp2,&
        mu2p,mucos,Pt,PI1,PI1_phi,PI2_phi,f3_1,f2_1,lambda_1,q_1,&
        sinobs_1,muobs_1,scal_1,pp_phi,p1_phi,p2_phi
        !parameter(zero=0.D0,one=1.D0,two=2.D0)
        integer  :: t1,t2,i,j,count_num=1
        logical :: mobseqmtp
        save :: PI1,PI1_phi,PI2_phi,Pt,f3_1,f2_1,lambda_1,q_1,pp_phi,p1_phi,p2_phi,&
        sinobs_1,muobs_1,scal_1,mobseqmtp,AA,BB,&
        mu_tp1,mu_tp2

        60 continue
        IF(count_num .EQ. 1)THEN
            f3_1=f3
            f2_1=f2
            lambda_1=lambda
            q_1=q
            muobs_1=muobs
            sinobs_1=sinobs
            scal_1=scal
            t1=0
            t2=0
            mobseqmtp=.false.
            If(q.gt.zero)then
                AA=sqrt((lambda**two+q)/q)
                BB=sqrt(q)
                !*****************************************************
                If(f2.lt.zero)then
                    mu=sin(asin(muobs*AA)+p*BB*AA)/AA
                else
                    If(f2.eq.zero)then
                        mu=cos(p*AA*BB)*muobs
                    else
                        mu=sin(asin(muobs*AA)-p*AA*BB)/AA
                    endif
                endif
                mucos = mu
                !****************************************************
                If(f2.ne.zero)then
                    mu_tp1=sqrt(q/(lambda**two+q))
                    mu_tp2=-mu_tp1
                else
                    mu_tp1=abs(muobs)
                    mu_tp2=-mu_tp1
                    mobseqmtp=.true.
                endif
                If(abs(muobs).eq.one)mobseqmtp=.true.

                If(mu_tp1.eq.zero)then
                    !photons are confined in the equatorial plane,
                    !so the integrations about !\theta are valished.
                    phyc_schwatz=zero
                    return
                endif

                !***************************************************
                PI1=(PI/two-asin(muobs/mu_tp1))*mu_tp1/BB
                Pt=PI*mu_tp1/BB
                pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB
                p1=PI1-pp
                p2=Pt-p1
                PI1_phi=zero
                PI2_phi=zero
                Do j=0,100
                    Do i=j,j+1
                        If(mobseqmtp)then
                            If(muobs.eq.mu_tp1)then
                                t1=j
                                t2=i
                                mu2p=-pp+two*(t1*p1+t2*p2)
                            else
                                t1=i
                                t2=j
                                mu2p=pp+two*(t1*p1+t2*p2)
                            endif
                        else
                            If(f2.lt.zero)then
                                t1=i
                                t2=j
                                mu2p=pp+two*(t1*p1+t2*p2)
                            endif
                            If(f2.gt.zero)then
                                t1=j
                                t2=i
                                mu2p=-pp+two*(t1*p1+t2*p2)
                            endif
                        endif
                        If(abs(p-mu2p).lt.1.D-4)goto 300
                    enddo
                enddo
                !***************************************************************
                300 continue
                If(lambda.eq.zero)then
                    phyc_schwatz = zero
                    return
                endif
                pp_phi=lambda*schwatz_int(muobs,mu,AA)/BB
                If(t1.eq.0)then
                    p1_phi=zero
                else
                    IF(PI1_phi .eq. zero)THEN
                        PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB
                    ENDIF
                    p1_phi=PI1_phi-pp_phi
                endif
                If(t2.eq.0)then
                    p2_phi=zero
                else
                    IF(PI2_phi .EQ. zero)THEN
                        PI2_phi=lambda*schwatz_int(mu_tp2,muobs,AA)/BB
                    ENDIF
                    p2_phi=PI2_phi+pp_phi
                endif
                If(mobseqmtp)then
                    If(muobs.eq.mu_tp1)then
                        phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    else
                        phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    endif
                else
                    If(f2.lt.zero)then
                        phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    endif
                    If(f2.gt.zero)then
                        phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    endif
                endif
            else
                !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is a affending',&
                !                'value, the program should be',&
                !                'stoped! and q = ',q
                !stop
                mucos=muobs
                t1 = 0
                t2 = 0
                phyc_schwatz = zero
            endif
        ELSE
            IF(f3_1.eq.f3.and.f2_1.eq.f2.and.lambda_1.eq.lambda.and.q_1.eq.q.and.sinobs_1.eq.sinobs&
                .and.muobs_1.eq.muobs.and.scal_1.eq.scal)THEN
                If(q.gt.zero)then
                    !*****************************************************
                    If(f2.lt.zero)then
                        mu=sin(asin(muobs*AA)+p*BB*AA)/AA
                    else
                        If(f2.eq.zero)then
                            mu=cos(p*AA*BB)*muobs
                        else
                            mu=sin(asin(muobs*AA)-p*AA*BB)/AA
                        endif
                    endif
                    mucos = mu
                    !****************************************************
                    If(mu_tp1.eq.zero)then
                        !photons are confined in the equatorial plane,
                        !so the integrations about !\theta are valished.
                        phyc_schwatz=zero
                        return
                    endif

                    !***************************************************
                    pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB
                    p1=PI1-pp
                    p2=Pt-p1
                    Do j=0,100
                        Do i=j,j+1
                            If(mobseqmtp)then
                                If(muobs.eq.mu_tp1)then
                                    t1=j
                                    t2=i
                                    mu2p=-pp+two*(t1*p1+t2*p2)
                                else
                                    t1=i
                                    t2=j
                                    mu2p=pp+two*(t1*p1+t2*p2)
                                endif
                            else
                                If(f2.lt.zero)then
                                    t1=i
                                    t2=j
                                    mu2p=pp+two*(t1*p1+t2*p2)
                                endif
                                If(f2.gt.zero)then
                                    t1=j
                                    t2=i
                                    mu2p=-pp+two*(t1*p1+t2*p2)
                                endif
                            endif
                            If(abs(p-mu2p).lt.1.D-4)goto 310
                        enddo
                    enddo
                    !************************************************************
                    310 continue
                    If(lambda.eq.zero)then
                        phyc_schwatz = zero
                        return
                    endif
                    pp_phi=lambda*schwatz_int(muobs,mu,AA)/BB
                    If(t1.eq.0)then
                        p1_phi=zero
                    else
                        IF(PI1_phi .eq. zero)THEN
                            PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB
                        ENDIF
                        p1_phi=PI1_phi-pp_phi
                    endif
                    If(t2.eq.0)then
                        p2_phi=zero
                    else
                        IF(PI2_phi .EQ. zero)THEN
                            PI2_phi=lambda*schwatz_int(mu_tp2,muobs,AA)/BB
                        ENDIF
                        p2_phi=PI2_phi+pp_phi
                    endif
                    If(mobseqmtp)then
                        If(muobs.eq.mu_tp1)then
                            phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        else
                            phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        endif
                    else
                        If(f2.lt.zero)then
                            phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        endif
                        If(f2.gt.zero)then
                            phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        endif
                    endif
                else
                    !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is a affending',&
                    !                'value, the program should be',&
                    !                'stoped! and q = ',q
                    !stop
                    mucos=muobs
                    t1 = 0
                    t2 = 0
                    phyc_schwatz = zero
                endif
            ELSE
                count_num=1
                goto 60
            ENDIF
        ENDIF
        return
    End SUBROUTINE phyt_schwatz
    !*************************************************************************
    Function schwatz_int(y,x,AA)
        !*************************************************************************
        !*     PURPOSE:  Computes \int^x_y dt/(1-t^2)/sqrt(1-AA^2*t^2) and AA .gt. 1
        !*     INPUTS:   components of above integration.
        !*     OUTPUTS:  valve of integral.
        !*     ROUTINES CALLED: NONE.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  5 Jan 2012
        !*     REVISIONS: ******************************************
        USE constants
        implicit none
        Double precision y,x,yt,xt,AA,schwatz_int,ppx,ppy,A2

        xt=x
        yt=y
        If(yt.eq.xt)then
            schwatz_int=0.D0
            return
        endif
        If(abs(AA).ne.one)then
            A2=AA*AA
            ppx=atan(sqrt(A2-one)*xt/sqrt(abs(one-A2*xt*xt)))
            ppy=atan(sqrt(A2-one)*yt/sqrt(abs(one-A2*yt*yt)))
            schwatz_int=(ppx-ppy)/sqrt(A2-one)
        ELse
            If(abs(xt).eq.one)then
                schwatz_int=infinity
            Else
                If(abs(yt).eq.one)then
                    schwatz_int=-infinity
                Else
                    ppx=xt/sqrt(abs(one-xt*xt))
                    ppy=yt/sqrt(abs(one-yt*yt))
                    schwatz_int=ppx-ppy
                endif
            Endif
        Endif
        return
    End Function schwatz_int

    !********************************************************************************************
    !   SUBROUTINE INTRPART(p,f1234r,f1234t,lambda,q,sinobs,muobs,a_spin,&
    !                         robs,scal,phyr,timer,affr,r_coord,t1,t2)
    SUBROUTINE INTRPART(p,f1234r,f1234t,lambda,q,sinobs,muobs,a_spin,lambdaBar,&
        robs,scal,phyr,timer,affr,r_coord,t1,t2)
        !********************************************************************************************
        !*     PURPOSE:  Computes r part of integrals in coordinates \phi, t and affine parameter \sigma,
        !*               expressed by equations (62), (63), (65), (67), (68) and (69) in Yang & Wang (2012).
        !*     INPUTS:   p--------------independent variable, which must be nonnegative.
        !*               f1234r---------p_r, which is the r component of four momentum of a photon
        !*                              measured under the LNRF, see equation (83) in Yang & Wang (2012).
        !*               f1234t---------p_\theta, which is the \theta component of four momentum of a photon
        !*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of observer or the initial position of photon.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  phyr-----------value of integral \phi_r expressed by equation (65) or (69) in
        !*                              Yang & Wang (2012).
        !*               affr-----------value of integral \sigma_r expressed by equation (62) or (67) in
        !*                              Yang & Wang (2012).
        !*               timer----------value of integral t_r expressed by equation (63) or (68) in
        !*                              Yang & Wang (2012).
        !*               r_coord--------value of function r(p).
        !*               t1,t2----------number of times of photon meets turning points r_tp1 and r_tp2
        !*                              respectively.
        !*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  5 Jan 2012
        !*     REVISIONS: ******************************************
        USE constants
        IMPLICIT NONE
        DOUBLE PRECISION phyr,radius,p,sinobs,muobs,a_spin,lambdaBar,rhorizon,q,lambda,integ4(4),&
        cc,b0,b1,b2,b3,g2,g3,tobs,tp,pp,p1,p2,PI0,E_add,E_m,&
        u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,r_add,r_m,B_add,B_m,D_add,D_m,&
        y,x,f1,g1,h1,f2,h2,a5,b5,a4,b4,robs,&
        scal,tinf,integ04(4),integ14(4),integ5(5),integ15(5),&
        r_tp1,r_tp2,t_inf,tp2,f1234r,f1234t,p_temp,PI0_obs_inf,PI0_total,PI0_obs_hori,&
        PI0_obs_tp2,PI01,timer,affr,r_coord,cr,dr,rff_p,&
        Ap,Am,h,wp,wm,wbarp,wbarm,hm,hp,pp_time,pp_phi,pp_aff,p1_phi,p1_time,p1_aff,&
        p2_phi,p2_time,p2_aff,time_temp,sqt3,p_tp1_tp2,PI2_p,PI1_p,&
        PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_aff,PI2_aff
        DOUBLE PRECISION f1234r_1,f1234t_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,scal_1
        ! 在 SUBROUTINE INTRPART 内部的开头部分添加：
        DOUBLE PRECISION :: dt_corr, dp_corr, da_corr
        !PARAMETER(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0)
        COMPLEX*16 bb(1:4),dd(3)
        INTEGER ::  reals,i,j,t1,t2,index_p4(4),index_p5(5),del,cases_int,cases,count_num=1
        LOGICAL :: robs_eq_rtp,indrhorizon
        SAVE :: f1234r_1,f1234t_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,scal_1,rhorizon,r_add,r_m,a4,b4,B_add,&
        B_m,robs_eq_rtp,indrhorizon,r_tp1,r_tp2,reals,cases,bb,b0,b1,b2,b3,g2,g3,tobs,thorizon,&
        tp2,tinf,dd,E_add,E_m,D_add,D_m,PI0,PI0_obs_inf,PI0_total,PI0_obs_hori,PI0_obs_tp2,del,&
        u,v,w,L1,L2,m2,t_inf,pinf,f1,g1,h1,f2,h2,b5,Ap,Am,h,wp,wm,wbarp,wbarm,hm,hp,a5,cc,&
        PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_aff,PI2_aff,PI2_p,PI1_p,p_tp1_tp2,sqt3

        40 continue
        If(count_num.eq.1)then
            ! initialize wormhole correction accumulators for this call
            dt_corr = 0.d0
            dp_corr = 0.d0
            da_corr = 0.d0
            f1234r_1=f1234r
            f1234t_1=f1234t
            lambda_1=lambda
            q_1=q
            muobs_1=muobs
            sinobs_1=sinobs
            a_spin_1=a_spin
            robs_1=robs
            scal_1=scal
            !************************************************************************************
            ! Wormhole throat radius (replaces Kerr horizon for termination)
            rhorizon = (one + lambdaBar**two) + sqrt((one + lambdaBar**two)**two - a_spin**two)
            ! equation (64) in Yang & Wang (2012) - Kerr horizon roots for partial fractions
            r_add=one+sqrt(one-a_spin**two)
            r_m=one-sqrt(one-a_spin**two)
            ! equation (64) in Yang & Wang (2012).
            B_add=(two*r_add-a_spin*lambda)/(r_add-r_m)
            B_m=(two*r_m-a_spin*lambda)/(r_add-r_m)
            ! equation (64) in Yang & Wang (2012).
            Ap=(r_add*(four-a_spin*lambda)-two*a_spin**two)/sqrt(one-a_spin**two)
            Am=(r_m*(four-a_spin*lambda)-two*a_spin**two)/sqrt(one-a_spin**two)
            b4=one
            a4=zero
            cc=a_spin**2-lambda**2-q
            robs_eq_rtp=.false.
            indrhorizon=.false.
            ! call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
            !                   reals,robs_eq_rtp,indrhorizon,cases,bb)
            call radiustp(f1234r,a_spin,lambdaBar,robs,lambda,q,r_tp1,r_tp2,&
            reals,robs_eq_rtp,indrhorizon,cases,bb)

            ! equation (55) in Yang & Wang (2012).
            PI1_phi=zero
            PI2_phi=zero
            PI1_time=zero
            PI2_time=zero
            PI1_aff=zero
            PI2_aff=zero
            !** R(r)=0 has real roots and turning points exists in radial r.
            If(reals.ne.0)then
                ! equations (35)-(38) in Yang & Wang (2012).
                ! Use pure Kerr quartic coefficients (turning points are Kerr).
                b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
                b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
                b2=four/three*r_tp1
                b3=one
                g2=three/four*(b1**2-b0*b2)
                g3=one/16.D0*(three*b0*b1*b2-two*b1**three-b0**two*b3)
                ! equation (39) in Yang & Wang (2012).
                If(robs-r_tp1.ne.zero)then
                    tobs=b0/four/(robs-r_tp1)+b1/four
                else
                    tobs=infinity
                endif
                If(rhorizon-r_tp1.ne.zero)then
                    thorizon=b1/four+b0/four/(rhorizon-r_tp1)
                else
                    thorizon=infinity
                endif
                tp2=b0/four/(r_tp2-r_tp1)+b1/four
                tinf=b1/four
                h=-b1/four
                ! equation (64), (66) and (70) in Yang & Wang (2012).
                call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
                ! Use Delta_WH horizons: r_add, r_m are roots of Delta_WH
                E_add=b0/(four*(r_add-r_tp1))+b1/four
                E_m=b0/(four*(r_m-r_tp1))+b1/four
                D_add=b0/(four*(r_tp1-r_add)**2)
                D_m=b0/(four*(r_tp1-r_m)**2)

                wp=one/(r_tp1-r_add)
                wm=one/(r_tp1-r_m)
                wbarp=b0/four/(r_tp1-r_add)**two
                wbarm=b0/four/(r_tp1-r_m)**two
                hp=b0/four/(r_add-r_tp1)+b1/four
                hm=b0/four/(r_m-r_tp1)+b1/four

                index_p4(1)=0
                cases_int=1
                call weierstrass_int_J3(tobs,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                ! equation (42) in Yang & Wang (2012).
                PI0=integ04(1)
                select case(cases)
                CASE(1)
                    If(f1234r .ge. zero)then !**photon will goto infinity.
                        index_p4(1)=0
                        cases_int=1
                        call weierstrass_int_J3(tinf,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        PI0_obs_inf=integ04(1)
                        If(p.lt.PI0_obs_inf)then
                            ! equation (41) in Yang & Wang (2012).
                            tp=weierstrassP(p+PI0,g2,g3,dd,del)
                            r_coord = r_tp1+b0/(four*tp-b1)
                            pp=-p
                        else
                            tp=tinf! !Goto infinity, far away.
                            r_coord = infinity
                            pp=-PI0_obs_inf
                        endif
                        t1=0
                        t2=0
                    ELSE
                        If(.not.indrhorizon)then
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            PI0_total=PI0+integ14(1)
                            t2=0
                            If(p.le.PI0)then
                                t1=0
                                pp=p
                                ! equation (41) in Yang & Wang (2012).
                                tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                r_coord = r_tp1+b0/(four*tp-b1)
                            else
                                t1=1
                                PI1_p=PI0
                                If(p.lt.PI0_total)then
                                    ! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=two*PI0-p
                                    p1=abs(p-PI0)
                                else
                                    tp=tinf !Goto infinity, far away.
                                    r_coord = infinity
                                    pp=-PI0_total+two*PI0
                                    p1=pI0_total-PI0
                                endif
                            endif
                        ELSE     !f1234r<0, photon will fall into black hole unless something encountered.
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            PI0_obs_hori=integ04(1)
                            If(p.lt.PI0_obs_hori)then
                                ! equation (41) in Yang & Wang (2012).
                                tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                r_coord = r_tp1+b0/(four*tp-b1)
                                pp=p
                            else
                                tp=thorizon! !Fall into black hole.
                                r_coord = rhorizon
                                pp=PI0_obs_hori
                            endif
                            t1=0
                            t2=0
                        ENDIF
                    ENDIF
                CASE(2)
                    If(.not.indrhorizon)then
                        If(f1234r.lt.zero)then
                            PI01=-PI0
                        else
                            PI01=PI0
                        endif
                        ! equation (41) in Yang & Wang (2012).
                        tp=weierstrassP(p+PI01,g2,g3,dd,del)
                        r_coord = r_tp1+b0/(four*tp-b1)
                        index_p4(1)=0
                        cases_int=1
                        call weierstrass_int_J3(tobs,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                        pp=integ4(1)
                        ! equation (57) in Yang & Wang (2012).
                        p_tp1_tp2=integ14(1)
                        PI2_p=p_tp1_tp2-PI0
                        PI1_p=PI0
                        p1=PI0-pp
                        p2=p_tp1_tp2-p1
                        !p1=zero
                        !p2=zero
                        !*************************************************************************************
                        ! equation (58) in Yang & Wang (2012).
                        Do j=0,100
                            Do i=j,j+1
                                If(robs_eq_rtp)then
                                    If(robs.eq.r_tp1)then
                                        t1=j
                                        t2=i
                                        p_temp=-pp+two*(t1*p1+t2*p2)
                                    else
                                        t1=i
                                        t2=j
                                        p_temp=pp+two*(t1*p1+t2*p2)
                                    endif
                                else
                                    If(f1234r.gt.zero)then
                                        t1=j
                                        t2=i
                                        p_temp=-pp+two*(t1*p1+t2*p2)
                                    endif
                                    If(f1234r.lt.zero)then
                                        t1=i
                                        t2=j
                                        p_temp=pp+two*(t1*p1+t2*p2)
                                    endif
                                endif
                                If(abs(p-p_temp).lt.1.D-4)goto 200
                            Enddo
                        Enddo
                        !*************************************************************************************
                        200     continue
                    else  !photon has probability to fall into black hole.
                        If(f1234r.le.zero)then
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            PI0_obs_hori=integ04(1)
                            If(p.lt.PI0_obs_hori)then
                                ! equation (41) in Yang & Wang (2012).
                                tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                r_coord = r_tp1+b0/(four*tp-b1)
                                pp=p
                            else
                                tp=thorizon! !Fall into black hole.
                                r_coord = rhorizon
                                pp=PI0_obs_hori
                            endif
                            t1=0
                            t2=0
                        ELSE  !p_r>0, photon will meet the r_tp2 turning point and turn around then goto vevnt horizon.
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tp2,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            call weierstrass_int_j3(tp2,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            PI0_obs_tp2=integ04(1)
                            PI2_p=PI0_obs_tp2
                            PI0_total=integ14(1)+PI0_obs_tp2
                            If(p.le.PI0_obs_tp2)then
                                t1=0
                                t2=0
                                pp=-p
                                ! equation (41) in Yang & Wang (2012).
                                tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                r_coord = r_tp1+b0/(four*tp-b1)
                            else
                                t1=0
                                t2=1
                                If(p.lt.PI0_total)then
                                    ! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p-two*PI0_obs_tp2
                                    p2=p-PI0_obs_tp2
                                else
                                    tp=thorizon !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_total-two*PI0_obs_tp2
                                    p2=PI0_total-PI0_obs_tp2
                                endif
                            endif
                        ENDIF
                    ENDIF
                END SELECT
                !******************************************************************
                index_p4(1)=-1
                index_p4(2)=-2
                index_p4(3)=0
                index_p4(4)=-4
                !pp part ***************************************************
                cases_int=4
                call weierstrass_int_J3(tobs,tp,dd,del,h,b4,index_p4,abs(pp),integ4,cases_int)

                ! equation (62) in Yang & Wang (2012).
                pp_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pp*r_tp1**two
                ! equation (63) in Yang & Wang (2012).
                pp_time=integ4(2)*b0/two+pp*(two*r_tp1+four+Ap*wp)+pp_aff
                time_temp=pp*(-Am*wm)

                cases_int=2
                call weierstrass_int_J3(tobs,tp,dd,del,-E_add,b4,index_p4,abs(pp),integ4,cases_int)
                ! equation (63) in Yang & Wang (2012).
                pp_time=pp_time-Ap*wbarp*integ4(2)
                IF(a_spin.NE.zero)THEN
                    call weierstrass_int_J3(tobs,tp,dd,del,-E_m,b4,index_p4,abs(pp),integ14,cases_int)
                    ! equation (63) in Yang & Wang (2012).
                    pp_time=pp_time+Am*wbarm*integ14(2)+time_temp
                    ! equation (65) in Yang & Wang (2012).
                    pp_phi=pp*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                    -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                ELSE
                    pp_phi=zero
                ENDIF
                IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                    ! equation (18) in Yang & Wang (2012).
                    pp_phi=pp_phi+pp*lambda
                ENDIF
                !p1 part *******************************************************
                IF(t1 .EQ. 0)THEN
                    p1_phi=ZERO
                    p1_time=ZERO
                    p1_aff=ZERO
                ELSE
                    IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN
                        cases_int=4
                        call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p4,PI0,integ4,cases_int)
                        ! equation (62) in Yang & Wang (2012).
                        PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+PI0*r_tp1**two
                        ! equation (63) in Yang & Wang (2012).
                        PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff
                        time_temp=PI0*(-Am*wm)

                        cases_int=2
                        call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,b4,index_p4,PI0,integ4,cases_int)
                        ! equation (63) in Yang & Wang (2012).
                        PI1_time=PI1_time-Ap*wbarp*integ4(2)
                        IF(a_spin.NE.zero)THEN
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,b4,index_p4,PI0,integ14,cases_int)
                            ! equation (63) in Yang & Wang (2012).
                            PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp
                            ! equation (65) in Yang & Wang (2012).
                            PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                            -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                        ELSE
                            PI1_phi=zero
                        ENDIF
                        IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                            ! equation (18) in Yang & Wang (2012).
                            PI1_phi=PI1_phi+PI0*lambda
                        ENDIF
                    ENDIF
                    ! equation (55) in Yang & Wang (2012).
                    p1_aff=PI1_aff-pp_aff
                    p1_time=PI1_time-pp_time
                    P1_phi=PI1_phi-pp_phi
                ENDIF
                !p2 part *******************************************************
                IF(t2.EQ.ZERO)THEN
                    p2_phi=ZERO
                    p2_time=ZERO
                    p2_aff=ZERO
                ELSE
                    IF(PI2_aff .EQ. zero .AND. PI2_time .EQ. zero)THEN
                        cases_int=4
                        call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p4,PI2_p,integ4,cases_int)
                        ! equation (62) in Yang & Wang (2012).
                        PI2_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+PI2_p*r_tp1**two
                        ! equation (63) in Yang & Wang (2012).
                        PI2_time=integ4(2)*b0/two+PI2_p*(two*r_tp1+four+Ap*wp)+PI2_aff
                        time_temp=PI2_p*(-Am*wm)

                        cases_int=2
                        call weierstrass_int_J3(tp2,tobs,dd,del,-E_add,b4,index_p4,PI2_p,integ4,cases_int)
                        ! equation (63) in Yang & Wang (2012).
                        PI2_time=PI2_time-Ap*wbarp*integ4(2)!+Am*wbarm*integ14(2)
                        IF(a_spin.NE.zero)THEN
                            call weierstrass_int_J3(tp2,tobs,dd,del,-E_m,b4,index_p4,PI2_p,integ14,cases_int)
                            ! equation (63) in Yang & Wang (2012).
                            PI2_time=PI2_time+Am*wbarm*integ14(2)+time_temp
                            ! equation (65) in Yang & Wang (2012).
                            PI2_phi=PI2_p*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                            -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                        ELSE
                            PI2_phi=zero
                        ENDIF
                        IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                            ! equation (18) in Yang & Wang (2012).
                            PI2_phi=PI2_phi+PI2_p*lambda
                        ENDIF
                    ENDIF
                    ! equation (55) in Yang & Wang (2012).
                    p2_aff=PI2_aff+pp_aff
                    p2_time=PI2_time+pp_time
                    p2_phi=PI2_phi+pp_phi
                ENDIF
                !phi, aff,time part *******************************************************
                ! equation (56) in Yang & Wang (2012).
                If(f1234r.ne.zero)then
                    phyr=sign(one,-f1234r)*pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timer=sign(one,-f1234r)*pp_time+two*(t1*p1_time+t2*p2_time)
                    affr=sign(one,-f1234r)*pp_aff+two*(t1*p1_aff+t2*p2_aff)
                else
                    If(robs.eq.r_tp1)then
                        phyr=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timer=-pp_time+two*(t1*p1_time+t2*p2_time)
                        affr=-pp_aff+two*(t1*p1_aff+t2*p2_aff)
                    else
                        phyr=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timer=pp_time+two*(t1*p1_time+t2*p2_time)
                        affr=pp_aff+two*(t1*p1_aff+t2*p2_aff)
                    endif
                endif
                !************************************************************************************************
                If(a_spin.eq.zero)then
                    If(cc.eq.zero)then
                        If(f1234r.lt.zero)then
                            If(p.lt.one/rhorizon-one/robs)then
                                radius=robs/(robs*p+one)
                            else
                                radius=rhorizon
                            endif
                        else
                            If(p.lt.one/robs)then
                                radius=robs/(one-robs*p)
                            else
                                radius=infinity
                            endif
                        endif
                    endif
                    If(cc.eq.-27.D0)then
                        sqt3=sqrt(three)
                        If(f1234r.lt.zero)then
                            cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                            dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                            If(p.ne.zero)then
                                radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                            else
                                radius=robs!infinity
                            endif
                        else
                            cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                            dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                            PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                            -Log(one+two/sqt3)/three/sqt3
                            If(p.lt.PI0)then
                                radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                            else
                                radius=infinity
                            endif
                        endif
                    endif
                endif
            ELSE
                ! equation (44) in Yang & Wang (2012).   !equation R(r)=0 has no real roots. we use the Legendre elliptic
                u=real(bb(4))        !integrations and functions to compute the calculations.
                w=abs(aimag(bb(4)))
                v=abs(aimag(bb(2)))
                If(u.ne.zero)then
                    ! equation (45) in Yang & Wang (2012).
                    L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                    L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                    ! equation (46) in Yang & Wang (2012).
                    thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
                    ! equation (48) in Yang & Wang (2012).
                    m2=(L1-L2)/L1
                    tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                    t_inf=sqrt((L1-one)/(L1-L2))
                    ! equation (50) in Yang & Wang (2012).
                    pinf=EllipticF(tinf,m2)/(w*sqrt(L1))
                    call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                    f1=u**2+w**2
                    g1=-two*u
                    h1=one
                    f2=u**2+v**2
                    g2=-g1
                    h2=one
                    a5=zero
                    b5=one
                    index_p5(1)=-1
                    index_p5(2)=-2
                    index_p5(3)=2
                    index_p5(4)=-4
                    index_p5(5)=4
                    IF(f1234r.lt.zero)THEN
                        PI0=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))
                        if(p.lt.PI0)then
                            ! equation (49) in Yang & Wang (2012).
                            y=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                            r_coord = y
                            pp=p
                        else
                            y=rhorizon
                            r_coord = y
                            pp=PI0
                        endif
                        x=robs
                    ELSE
                        PI0=EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                        if(p.lt.PI0)then
                            ! equation (49) in Yang & Wang (2012).
                            x=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                            r_coord = x
                            pp=p
                        else
                            x=infinity
                            r_coord = x
                            pp=PI0
                        endif
                        y=robs
                    ENDIF
                    !affine parameter part integration **********************************************
                    cases_int=5
                    call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,a5,b5,index_p5,abs(pp),integ5,cases_int)
                    ! equation (67) in Yang & Wang (2012).
                    affr=integ5(5)
                    ! equation (68) in Yang & Wang (2012).
                    timer=two*integ5(3)+four*pp+affr
                    cases_int=2
                    call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_add,b5,index_p5,abs(pp),integ5,cases_int)
                    call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_m,b5,index_p5,abs(pp),integ15,cases_int)
                    !phy part**************************************************************************
                    ! equation (68) in Yang & Wang (2012).
                    timer=timer+Ap*integ5(2)-Am*integ15(2)
                    ! equation (69) in Yang & Wang (2012).
                    phyr=a_spin*(B_add*integ5(2)-B_m*integ15(2))
                    IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                        ! equation (18) in Yang & Wang (2012).
                        phyr=phyr+pp*lambda
                    ENDIF
                ELSE
                    If(f1234r.lt.zero)then
                        PI0=(atan(robs/w)-atan(rhorizon/w))/w
                        If(p.lt.PI0)then
                            radius=w*tan(atan(robs/w)-p*w)
                            r_coord = radius
                        else
                            radius=rhorizon
                            r_coord = radius
                        endif
                        !timer part ****************************************
                        y=radius
                        x=robs
                    ELSE
                        PI0=(PI/two-atan(robs/w))/w
                        If(p.lt.PI0)then
                            radius=w*tan(atan(robs/w)+p*w)
                            r_coord = radius
                        else
                            radius=infinity
                            r_coord = radius
                        endif
                        !timer part ************************************************************************8
                        y=robs
                        x=radius
                    ENDIF
                    pp_time=(x-y)+atan(x/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))-&
                    atan(y/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))
                    pp_time=pp_time+log(x**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two))-&
                    (log(y**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two)))
                    timer=pp_time+Ap*log(abs(x-r_add))/(w**two+r_add**two)-Am*log(abs(x-r_m))/(w**two+r_m**two)-&
                    (Ap*log(abs(y-r_add))/(w**two+r_add**two)-Am*log(abs(y-r_m))/(w**two+r_m**two))
                    !affine parameter part *****************************************************************
                    affr=(x-y)-w*atan(x/w)+w*atan(y/w)
                    !phy part ******************************************************************
                    IF(a_spin .NE. zero)THEN
                        phyr=(-B_add*r_add/w/(r_add**two+w**two)+B_m*r_m/w/(r_m**two+w**two))&
                        *(atan(x/w)-atan(y/w))+&
                        log(abs(x-r_add)/sqrt(x**two+w**two))*B_add/(r_add*two+w**two)-&
                        log(abs(y-r_add)/sqrt(y**two+w**two))*B_add/(r_add*two+w**two)-&
                        log(abs(x-r_m)/sqrt(x**two+w**two))*B_m/(r_m*two+w**two)+&
                        log(abs(y-r_m)/sqrt(y**two+w**two))*B_m/(r_m*two+w**two)
                        phyr=phyr*a_spin
                    ELSE
                        phyr=zero
                    ENDIF
                    If(muobs.eq.zero.and.f1234t.eq.zero)then
                        phyr=phyr+lambda*(atan(x/w)-atan(y/w))/w
                    ENDIF
                ENDIF
            ENDIF

            ! --- Apply wormhole g_rr correction (Simpson quadrature for sqrt(DeltaK/DeltaWH)-1 factor) ---
            IF (lambdaBar .GT. 1.d-10) THEN
                CALL INTRPART_WH_CORR(robs, r_coord, r_tp1, r_tp2, t1, t2, cases, reals, indrhorizon, &
                f1234r, lambda, q, a_spin, lambdaBar, dt_corr, dp_corr, da_corr)
                timer = timer + dt_corr
                phyr  = phyr  + dp_corr
                affr  = affr  + da_corr
            END IF
            count_num=count_num+1
        !*****************************************************************************************
        else
            !*****************************************************************************************
            dt_corr = 0.d0
            dp_corr = 0.d0
            da_corr = 0.d0

            If(f1234r.eq.f1234r_1.and.f1234t.eq.f1234t_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
                sinobs.eq.sinobs_1.and.muobs.eq.muobs_1.and.a_spin.eq.a_spin_1.and.robs.eq.robs_1.and.scal.eq.scal_1)then
                !***********************************************************************
                If(reals.ne.0)then  !** R(r)=0 has real roots and turning points exists in radial r.
                    !used in the geodesics I'm trying to understand March 6 2017
                    select case(cases)
                    CASE(1)
                        If(f1234r .ge. zero)then !**photon will goto infinity.
                            If(p.lt.PI0_obs_inf)then
                                ! equation (41) in Yang & Wang (2012).
                                tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                r_coord = r_tp1+b0/(four*tp-b1)
                                pp=-p
                            else
                                tp=tinf! !Goto infinity, far away.
                                r_coord = infinity
                                pp=-PI0_obs_inf
                            endif
                            t1=0
                            t2=0
                        ELSE
                            If(.not.indrhorizon)then
                                t2=0
                                If(p.le.PI0)then
                                    t1=0
                                    pp=p
                                    ! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                else
                                    t1=1
                                    If(p.lt.PI0_total)then
                                        ! equation (41) in Yang & Wang (2012).
                                        tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=two*PI0-p
                                        p1=abs(p-PI0)
                                    else
                                        tp=tinf !Goto infinity, far away.
                                        r_coord = infinity
                                        pp=-PI0_total+two*PI0
                                        p1=pI0_total-PI0
                                    endif
                                endif
                            ELSE     !f1234r<0, photon will fall into black hole unless something encountered.
                                If(p.lt.PI0_obs_hori)then
                                    ! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0
                            ENDIF
                        ENDIF
                    CASE(2)
                        If(.not.indrhorizon)then
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0
                            endif
                            ! equation (41) in Yang & Wang (2012).
                            tp=weierstrassP(p+PI01,g2,g3,dd,del)
                            r_coord = r_tp1+b0/(four*tp-b1)
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tobs,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            ! equation (57) in Yang & Wang (2012).
                            pp=integ4(1)
                            p1=PI0-pp
                            p2=p_tp1_tp2-p1
                            !*************************************************************************************
                            ! equation (58) in Yang & Wang (2012).
                            Do j=0,100
                                Do i=j,j+1
                                    If(robs_eq_rtp)then
                                        If(robs.eq.r_tp1)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        else
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    else
                                        If(f1234r.gt.zero)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        endif
                                        If(f1234r.lt.zero)then
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    endif
                                    If(abs(p-p_temp).lt.1.D-4)goto 210
                                Enddo
                            Enddo
                            !*************************************************************************************
                            210        continue
                        else  !photon has probability to fall into black hole.
                            If(f1234r.le.zero)then
                                If(p.lt.PI0_obs_hori)then
                                    ! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0
                            ELSE  !p_r>0, photon will meet the r_tp2 turning point and turn around then goto vevnt horizon.
                                If(p.le.PI0_obs_tp2)then
                                    t1=0
                                    t2=0
                                    pp=-p
                                    ! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                else
                                    t1=0
                                    t2=1
                                    If(p.lt.PI0_total)then
                                        ! equation (41) in Yang & Wang (2012).
                                        tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=p-two*PI0_obs_tp2
                                        p2=p-PI0_obs_tp2
                                    else
                                        tp=thorizon !Fall into black hole.
                                        r_coord = rhorizon
                                        pp=PI0_total-two*PI0_obs_tp2
                                        p2=PI0_total-PI0_obs_tp2
                                    endif
                                endif
                            ENDIF
                        ENDIF
                    END SELECT
                    !******************************************************************
                    index_p4(1)=-1
                    index_p4(2)=-2
                    index_p4(3)=0
                    index_p4(4)=-4

                    !pp part ***************************************************
                    cases_int=4
                    call weierstrass_int_J3(tobs,tp,dd,del,h,b4,index_p4,abs(pp),integ4,cases_int)
                    ! equation (62) in Yang & Wang (2012).
                    pp_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pp*r_tp1**two
                    ! equation (63) in Yang & Wang (2012).
                    pp_time=integ4(2)*b0/two+pp*(two*r_tp1+four+Ap*wp)+pp_aff
                    time_temp=pp*(-Am*wm)

                    cases_int=2
                    call weierstrass_int_J3(tobs,tp,dd,del,-E_add,b4,index_p4,abs(pp),integ4,cases_int)
                    ! equation (63) in Yang & Wang (2012).
                    pp_time=pp_time-Ap*wbarp*integ4(2)
                    IF(a_spin.NE.zero)THEN
                        call weierstrass_int_J3(tobs,tp,dd,del,-E_m,b4,index_p4,abs(pp),integ14,cases_int)
                        ! equation (63) in Yang & Wang (2012).
                        pp_time=pp_time+Am*wbarm*integ14(2)+time_temp
                        ! equation (65) in Yang & Wang (2012).
                        pp_phi=pp*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                        -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                    ELSE
                        pp_phi=zero
                    ENDIF
                    IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                        ! equation (18) in Yang & Wang (2012).
                        pp_phi=pp_phi+pp*lambda
                    ENDIF
                    !p1 part *******************************************************
                    IF(t1 .EQ. 0)THEN
                        p1_phi=ZERO
                        p1_time=ZERO
                        p1_aff=ZERO

                        !added by Pieter ************************************************************************
                        IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN !original
                            !here, PI1_time is calculated the first point after the turning point (t1>0)
                            cases_int=4
                            call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p4,PI0,integ4,cases_int)

                            ! equation (62) in Yang & Wang (2012).
                            !********************************************* old ******************************************
                            !PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+p1*r_tp1**two
                            !********************************************* old ******************************************
                            !********************************************* new ******************************************
                            PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pI0*r_tp1**two
                            !********************************************* new ******************************************
                            ! equation (63) in Yang & Wang (2012).
                            PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff
                            time_temp=PI0*(-Am*wm)

                            cases_int=2
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,b4,index_p4,PI0,integ4,cases_int)
                            ! equation (63) in Yang & Wang (2012).
                            PI1_time=PI1_time-Ap*wbarp*integ4(2)
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,b4,index_p4,PI0,integ14,cases_int)
                                ! equation (65) in Yang & Wang (2012).
                                PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp
                                PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI1_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                                ! equation (18) in Yang & Wang (2012).
                                PI1_phi=PI1_phi+PI0*lambda
                            ENDIF
                        ENDIF
                        p1_aff=PI1_aff-pp_aff
                        p1_time=PI1_time-pp_time
                        P1_phi=PI1_phi-pp_phi



                    !end added by Pieter ************************************************************************
                    ELSE
                        IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN !original
                            !here, PI1_time is calculated the first point after the turning point (t1>0)
                            cases_int=4
                            call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p4,PI0,integ4,cases_int)

                            ! equation (62) in Yang & Wang (2012).
                            !********************************************* old ******************************************
                            !PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+p1*r_tp1**two
                            !********************************************* old ******************************************
                            !********************************************* new ******************************************
                            PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pI0*r_tp1**two
                            !********************************************* new ******************************************
                            ! equation (63) in Yang & Wang (2012).
                            PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff
                            time_temp=PI0*(-Am*wm)

                            cases_int=2
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,b4,index_p4,PI0,integ4,cases_int)
                            ! equation (63) in Yang & Wang (2012).
                            PI1_time=PI1_time-Ap*wbarp*integ4(2)
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,b4,index_p4,PI0,integ14,cases_int)
                                ! equation (65) in Yang & Wang (2012).
                                PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp
                                PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI1_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                                ! equation (18) in Yang & Wang (2012).
                                PI1_phi=PI1_phi+PI0*lambda
                            ENDIF
                        !write(*,*)'hier:',PI1_time
                        ENDIF
                        ! equation (55) in Yang & Wang (2012).
                        p1_aff=PI1_aff-pp_aff
                        p1_time=PI1_time-pp_time
                        P1_phi=PI1_phi-pp_phi
                    ENDIF
                    !p2 part *******************************************************
                    IF(t2.EQ.ZERO)THEN
                        p2_phi=ZERO
                        p2_time=ZERO
                        p2_aff=ZERO
                    ELSE
                        IF(PI2_aff .EQ. zero .AND. PI2_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p4,PI2_p,integ4,cases_int)
                            ! equation (62) in Yang & Wang (2012).
                            PI2_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+PI2_p*r_tp1**two
                            ! equation (63) in Yang & Wang (2012).
                            PI2_time=integ4(2)*b0/two+PI2_p*(two*r_tp1+four+Ap*wp)+PI2_aff
                            time_temp=PI2_p*(-Am*wm)

                            cases_int=2
                            call weierstrass_int_J3(tp2,tobs,dd,del,-E_add,b4,index_p4,PI2_p,integ4,cases_int)
                            ! equation (63) in Yang & Wang (2012).
                            PI2_time=PI2_time-Ap*wbarp*integ4(2)!+Am*wbarm*integ14(2)
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tp2,tobs,dd,del,-E_m,b4,index_p4,PI2_p,integ14,cases_int)
                                ! equation (63) in Yang & Wang (2012).
                                PI2_time=PI2_time+Am*wbarm*integ14(2)+time_temp
                                ! equation (65) in Yang & Wang (2012).
                                PI2_phi=PI2_p*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI2_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                                ! equation (18) in Yang & Wang (2012).
                                PI2_phi=PI2_phi+PI2_p*lambda
                            ENDIF
                        ENDIF
                        ! equation (55) in Yang & Wang (2012).
                        p2_aff=PI2_aff+pp_aff
                        p2_time=PI2_time+pp_time
                        p2_phi=PI2_phi+pp_phi
                    ENDIF
                    !phi, aff,time part *******************************************************
                    ! equation (56) in Yang & Wang (2012).
                    If(f1234r.ne.zero)then
                        phyr=sign(one,-f1234r)*pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timer=sign(one,-f1234r)*pp_time+two*(t1*p1_time+t2*p2_time)
                        affr=sign(one,-f1234r)*pp_aff+two*(t1*p1_aff+t2*p2_aff)
                    else
                        If(robs.eq.r_tp1)then
                            phyr=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=-pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=-pp_aff+two*(t1*p1_aff+t2*p2_aff)
                        else
                            phyr=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=pp_aff+two*(t1*p1_aff+t2*p2_aff)
                        endif
                    endif
                    !************************************************************************************************
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                If(p.lt.one/rhorizon-one/robs)then
                                    radius=robs/(robs*p+one)
                                else
                                    radius=rhorizon
                                endif
                            else
                                If(p.lt.one/robs)then
                                    radius=robs/(one-robs*p)
                                else
                                    radius=infinity
                                endif
                            endif
                        endif
                        If(cc.eq.-27.D0)then
                            sqt3=sqrt(three)
                            If(f1234r.lt.zero)then
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                                If(p.ne.zero)then
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=robs!infinity
                                endif
                            else
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                                PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                -Log(one+two/sqt3)/three/sqt3
                                If(p.lt.PI0)then
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=infinity
                                endif
                            endif
                        endif
                    endif
                ELSE                        !equation R(r)=0 has no real roots. we use the Legendre elliptic
                    !integrations and functions to compute the calculations.
                    If(u.ne.zero)then
                        call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                        index_p5(1)=-1
                        index_p5(2)=-2
                        index_p5(3)=2
                        index_p5(4)=-4
                        index_p5(5)=4
                        IF(f1234r.lt.zero)THEN
                            if(p.lt.PI0)then
                                ! equation (49) in Yang & Wang (2012).
                                y=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = y
                                pp=p
                            else
                                y=rhorizon
                                r_coord = y
                                pp=PI0
                            endif
                            x=robs
                        ELSE
                            if(p.lt.PI0)then
                                ! equation (49) in Yang & Wang (2012).
                                x=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = x
                                pp=p
                            else
                                x=infinity
                                r_coord = x
                                pp=PI0
                            endif
                            y=robs
                        ENDIF
                        !affine parameter part integration **********************************************
                        cases_int=5
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,a5,b5,index_p5,abs(pp),integ5,cases_int)
                        ! equation (67) in Yang & Wang (2012).
                        affr=integ5(5)
                        ! equation (68) in Yang & Wang (2012).
                        timer=two*integ5(3)+four*pp+affr
                        cases_int=2
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_add,b5,index_p5,abs(pp),integ5,cases_int)
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_m,b5,index_p5,abs(pp),integ15,cases_int)
                        !phy part**************************************************************************
                        ! equation (68) in Yang & Wang (2012).
                        timer=timer+Ap*integ5(2)-Am*integ15(2)
                        ! equation (69) in Yang & Wang (2012).
                        phyr=a_spin*(B_add*integ5(2)-B_m*integ15(2))
                        IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
                            ! equation (18) in Yang & Wang (2012).
                            phyr=phyr+pp*lambda
                        ENDIF
                    ELSE
                        If(f1234r.lt.zero)then
                            If(p.lt.PI0)then
                                radius=w*tan(atan(robs/w)-p*w)
                                r_coord = radius
                            else
                                radius=rhorizon
                                r_coord = radius
                            endif
                            !timer part ******************************************************************
                            y=radius
                            x=robs
                        ELSE
                            If(p.lt.PI0)then
                                radius=w*tan(atan(robs/w)+p*w)
                                r_coord = radius
                            else
                                radius=infinity
                                r_coord = radius
                            endif
                            !timer part ************************************************************************8
                            y=robs
                            x=radius
                        ENDIF
                        pp_time=(x-y)+atan(x/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))-&
                        atan(y/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))
                        pp_time=pp_time+log(x**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two))-&
                        (log(y**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two)))
                        timer=pp_time+Ap*log(abs(x-r_add))/(w**two+r_add**two)-Am*log(abs(x-r_m))/(w**two+r_m**two)-&
                        (Ap*log(abs(y-r_add))/(w**two+r_add**two)-Am*log(abs(y-r_m))/(w**two+r_m**two))
                        !affine parameter part *****************************************************************
                        affr=(x-y)-w*atan(x/w)+w*atan(y/w)
                        !phy part ******************************************************************
                        IF(a_spin .NE. zero)THEN
                            phyr=(-B_add*r_add/w/(r_add**two+w**two)+B_m*r_m/w/(r_m**two+w**two))&
                            *(atan(x/w)-atan(y/w))+&
                            log(abs(x-r_add)/sqrt(x**two+w**two))*B_add/(r_add*two+w**two)-&
                            log(abs(y-r_add)/sqrt(y**two+w**two))*B_add/(r_add*two+w**two)-&
                            log(abs(x-r_m)/sqrt(x**two+w**two))*B_m/(r_m*two+w**two)+&
                            log(abs(y-r_m)/sqrt(y**two+w**two))*B_m/(r_m*two+w**two)
                            phyr=phyr*a_spin
                        ELSE
                            phyr=zero
                        ENDIF
                        If(muobs.eq.zero.and.f1234t.eq.zero)then
                            phyr=phyr+lambda*(atan(x/w)-atan(y/w))/w
                        ENDIF
                    ENDIF
                ENDIF

                ! --- Apply wormhole g_rr correction for repeated calls ---
                IF (lambdaBar .GT. 1.d-10) THEN
                    CALL INTRPART_WH_CORR(robs, r_coord, r_tp1, r_tp2, t1, t2, cases, reals, indrhorizon, &
                        f1234r, lambda, q, a_spin, lambdaBar, dt_corr, dp_corr, da_corr)
                    timer = timer + dt_corr
                    phyr  = phyr  + dp_corr
                    affr  = affr  + da_corr
                END IF

            ELSE
                count_num=1
                goto 40
            endif
        endif

        RETURN
    END SUBROUTINE INTRPART



    !======================================================================
    !  Wormhole correction helpers (g_rr only): p is treated as Kerr.
    !  r(p) mapping is unchanged; only (t_r, phi_r, sigma_r) radial integrals
    !  receive an additive correction:
    !     delta I = ∫ (sqrt(Delta/Delta_hat) - 1) * (Kerr integrand) dr
    !======================================================================
    SUBROUTINE INTRPART_WH_CORR(robs, r_coord, r_tp1, r_tp2, t1, t2, cases, reals, indrhorizon, &
        f1234r, lam_ph, q, a, l_wh, dt_corr, dp_corr, da_corr)
        USE constants
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN)  :: robs, r_coord, r_tp1, r_tp2, f1234r, lam_ph, q, a, l_wh
        INTEGER,          INTENT(IN)  :: t1, t2, cases, reals
        LOGICAL,          INTENT(IN)  :: indrhorizon
        DOUBLE PRECISION, INTENT(OUT) :: dt_corr, dp_corr, da_corr

        DOUBLE PRECISION :: r_end_eff, dt, dp, da, r_turn
        DOUBLE PRECISION, PARAMETER :: RINF = 1.d5
        INTEGER :: tot_turns

        dt_corr = 0.d0; dp_corr = 0.d0; da_corr = 0.d0

        IF (l_wh .LE. 1.d-10) RETURN
        IF (indrhorizon) RETURN 

        ! 处理无穷远
        IF (r_coord .GE. 0.5d0*infinity) THEN
            r_end_eff = RINF
        ELSE
            r_end_eff = r_coord
        END IF

        ! =========================================================
        ! 修正后的路径补偿逻辑 (针对铁线谱优化)
        ! =========================================================
        
        IF (reals .NE. 0) THEN
            tot_turns = t1 + t2
            
            ! -----------------------------------------------------
            ! 情况 A: 单调路径 (无转折)
            ! -----------------------------------------------------
            IF (tot_turns .EQ. 0) THEN
                CALL WH_PATH_COMPENSATION(robs, r_end_eff, lam_ph, q, a, l_wh, dt_corr, dp_corr, da_corr)

            ! -----------------------------------------------------
            ! 情况 B: 单次转折 (最常见的一级像)
            ! -----------------------------------------------------
            ELSE IF (tot_turns .EQ. 1) THEN
                ! 确定转折点是 r_tp1 还是 r_tp2
                IF (t1 .EQ. 1) THEN
                    r_turn = r_tp1
                ELSE
                    r_turn = r_tp2
                END IF

                ! 第一段：观测者 -> 转折点
                CALL WH_PATH_COMPENSATION(robs, r_turn, lam_ph, q, a, l_wh, dt, dp, da)
                dt_corr = dt_corr + dt; dp_corr = dp_corr + dp; da_corr = da_corr + da

                ! 第二段：转折点 -> 发射点(吸积盘)
                CALL WH_PATH_COMPENSATION(r_turn, r_end_eff, lam_ph, q, a, l_wh, dt, dp, da)
                dt_corr = dt_corr + dt; dp_corr = dp_corr + dp; da_corr = da_corr + da

            ! -----------------------------------------------------
            ! 情况 C: 多次转折 (高阶光子环)
            ! -----------------------------------------------------
            ELSE
                ! 对于铁线谱，这部分贡献极小。
                ! 为了避免代码崩溃，我们做一个近似：
                ! 假设光子主要在转折点附近徘徊，我们只计算最后一次从转折点飞向盘面的过程。
                ! 或者，如果你想更严谨一点但不想写循环，可以简单地由两段组成近似路径。
                ! 这里采用最安全的“忽略中间震荡，只算两端”的策略，误差在高阶像中可接受。
                
                IF (t1 .GE. 1) r_turn = r_tp1
                IF (t2 .GE. 1) r_turn = r_tp2 ! 取最靠近盘面的那个转折点
                
                ! 近似处理：分两段算，忽略中间的来回震荡 (Loop correction ignored)
                CALL WH_PATH_COMPENSATION(robs, r_turn, lam_ph, q, a, l_wh, dt, dp, da)
                dt_corr = dt_corr + dt; dp_corr = dp_corr + dp; da_corr = da_corr + da

                CALL WH_PATH_COMPENSATION(r_turn, r_end_eff, lam_ph, q, a, l_wh, dt, dp, da)
                dt_corr = dt_corr + dt; dp_corr = dp_corr + dp; da_corr = da_corr + da
            END IF

        ELSE
            ! 复根情况 (无转折)
            CALL WH_PATH_COMPENSATION(robs, r_end_eff, lam_ph, q, a, l_wh, dt_corr, dp_corr, da_corr)
        END IF

    END SUBROUTINE INTRPART_WH_CORR

    SUBROUTINE WH_PATH_COMPENSATION(r_start, r_end, lam_ph, q, a, l_wh, d_timer, d_phyr, d_affr)
        USE constants
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN)  :: r_start, r_end, lam_ph, q, a, l_wh
        DOUBLE PRECISION, INTENT(OUT) :: d_timer, d_phyr, d_affr

        INTEGER, PARAMETER :: N = 200   ! must be even for Simpson
        DOUBLE PRECISION :: r1, r2, h, r, wgt
        DOUBLE PRECISION :: Delta, Delta_h, P, Rk, sqRk, ratio, factor
        DOUBLE PRECISION :: sum_t, sum_p, sum_a
        DOUBLE PRECISION :: it, ip, ia
        DOUBLE PRECISION :: eps, epsR
        INTEGER :: i

        d_timer = 0.d0
        d_phyr  = 0.d0
        d_affr  = 0.d0

        IF (ABS(r_end - r_start) .LT. 1.d-14) RETURN

        r1 = r_start
        r2 = r_end
        h  = (r2 - r1) / DBLE(N)

        eps  = 1.d-7
        epsR = 1.d-30

        sum_t = 0.d0
        sum_p = 0.d0
        sum_a = 0.d0

        DO i = 0, N
            r = r1 + h*DBLE(i)

            ! Nudge endpoints away from exact singular surfaces
            IF (i .EQ. 0) r = r + SIGN(eps, h)
            IF (i .EQ. N) r = r - SIGN(eps, h)

            Delta   = r*r - 2.d0*r + a*a
            Delta_h = r*r - 2.d0*(1.d0 + l_wh*l_wh)*r + a*a

            ! Avoid invalid regions; if ratio <= 0, treat correction as zero (outside domain)
            ratio = Delta / Delta_h
            IF (ratio .LE. 0.d0) THEN
                factor = 0.d0
            ELSE
                factor = SQRT(ratio) - 1.d0
            END IF

            P  = (r*r + a*a) - a*lam_ph
            ! Use pure Kerr Delta for the radial potential R_K(r).
            ! The wormhole correction is already captured by the 'factor' multiplier.
            Rk = P*P - Delta * ( q + (lam_ph - a)*(lam_ph - a) )

            IF (Rk .LT. 0.d0) Rk = 0.d0
            sqRk = SQRT(Rk + epsR)

            ! If Delta is extremely small, skip this point (near horizon); avoid blow-up
            IF (ABS(Delta) .LT. 1.d-14) CYCLE

            it = factor * ( ((r*r + a*a)*P/Delta) - a*(a - lam_ph) ) / sqRk
            ip = factor * ( (a*P/Delta) - (a - lam_ph) ) / sqRk
            ia = factor * ( r*r ) / sqRk

            IF (i .EQ. 0 .OR. i .EQ. N) THEN
                wgt = 1.d0
            ELSEIF (MOD(i,2) .EQ. 1) THEN
                wgt = 4.d0
            ELSE
                wgt = 2.d0
            ENDIF

            sum_t = sum_t + wgt*it
            sum_p = sum_p + wgt*ip
            sum_a = sum_a + wgt*ia
        END DO

        d_timer = (h/3.d0) * sum_t
        d_phyr  = (h/3.d0) * sum_p
        d_affr  = (h/3.d0) * sum_a

    END SUBROUTINE WH_PATH_COMPENSATION


    !********************************************************************************************
    !   Function  Pemdisk(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mu,rout,rin)
    Function  Pemdisk(f1234,lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,scal,mu,rout,rin)
        !********************************************************************************************
        !*     PURPOSE:  Solves equation \mu(p)=mu, i.e. to search the value p_{em} of
        !*               parameter p, corresponding to the intersection point of geodesic with with
        !*               disk, i.e. a surface has a constants inclination angle with respect to
        !*               equatorial plane of black hole.
        !*
        !*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are defined by equation
        !*                              (82)-(85) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of observer or initial position of photon.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination
        !*                              angle of disk surface with respect to the equatorial plane of
        !*                              black hole.
        !*               rin, rout------inner and outer radius of disk.
        !*     OUTPUTS:  pemdisk--------value of root of equation \mu(p)= mu for p.
        !*                              pemdisk=-1.D0, if the photon goto infinity.
        !*                              pemdisk=-2.D0, if the photon fall into event horizon.
        !*     REMARKS:                 This routine just search the intersection points of geodesic with
        !*                              up surface of disk. Following routine Pemdisk_all will searches
        !*                              intersection points of geodesic with up and down surface of disk.
        !*     ROUTINES CALLED: mutp, mu2p.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  6 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision Pemdisk,f1234(4),f3,f2,sinobs,muobs,a_spin,lambdaBar,lambda,q,mu_tp,two,&
        mu,scal,zero,robs,&
        mu_tp2,four,one,pm,rout,rin,re
        parameter(zero=0.D0,two=2.D0,four=4.D0,one=1.D0)
        integer  t1,t2,reals
        logical :: mobseqmtp

        f3 = f1234(3)
        f2 = f1234(2)
        mobseqmtp=.false.
        call mutp(f2,f3,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp)

        If(reals.eq.2)then
            If(mobseqmtp)then
                t1=0
                t2=0
            else
                If(muobs.gt.zero)then
                    If(f2.lt.zero)then
                        t1=1
                        t2=0
                    endif
                    If(f2.gt.zero)then
                        t1=0
                        t2=0
                    endif
                else
                    If(muobs.eq.zero)then
                        If(f2.lt.zero)then
                            t1=1
                            t2=0
                        endif
                        If(f2.gt.zero)then
                            t1=0
                            t2=1
                        endif
                    else
                        If(f2.lt.zero)then
                            t1=0
                            t2=0
                        endif
                        If(f2.gt.zero)then
                            t1=0
                            t2=1
                        endif
                    endif
                endif
            endif
            !write(*,*)'mu=',mu,B,t1,t2
            pm=mu2p(f3,f2,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal)
            ! re = radius(pm,f1234(1),lambda,q,a_spin,robs,scal)
            re = radius(pm,f1234(1),lambda,q,a_spin,lambdaBar,robs,scal)
            IF(re .le. rout .and. re .ge. rin)THEN
                Pemdisk = pm
                return
            ENDIF
            IF(re .gt. rout)THEN
                Pemdisk = -one
                RETURN
            ENDIF
            IF(re .lt. rin)THEN
                Pemdisk = -two
                return
            ENDIF
        else
            Pemdisk=-one
        endif
        return
    End Function Pemdisk
    !*********************************************************************
    Function  Pemdisk_all(f1234,lambda,q,sinobs,muobs,a_spin,lambdaBar,robs,scal,mu,rout,rin)
        !*********************************************************************
        !*     PURPOSE:  Solves equation \mu(p)=mu, where \mu(p)=\mu(p), i.e. to search the value p_{em} of
        !*               parameter p, corresponding to the intersection point of geodesic with
        !*               disk, i.e. a surface has a constants inclination angle with respect to
        !*               equatorial plane of black hole.
        !*
        !*     INPUTS:   f1234(1:4)-----array of f_1, f_2, f_3, f_0, which are defined by equation
        !*                              (106)-(109) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of observer or initial position of photon.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination
        !*                              angle of disk surface with respect to the equatorial plane of
        !*                              black hole.
        !*               rin, rout------inner and outer radius of disk.
        !*     OUTPUTS:  pemdisk--------value of root of equation \mu(p)= mu for p.
        !*                              pemdisk=-1.D0, if the photon goto infinity.
        !*                              pemdisk=-2.D0, if the photon fall into event horizon.
        !*     REMARKS:                 This routine will searches intersection points of
        !*                              geodesic with double surfaces of disk.
        !*     ROUTINES CALLED: mutp, mu2p, radius.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  6 Jan 2012
        !*     REVISIONS: ******************************************
        implicit none
        Double precision Pemdisk_all,f3,f2,sinobs,muobs,a_spin,lambdaBar,lambda,q,mu_tp1,two,&
        mu,scal,robs,zero,&
        mu_tp2,four,one,pm,f1234(4),rout,rin,re,rhorizon
        parameter(zero=0.D0,two=2.D0,four=4.D0,one=1.D0)
        integer  t1,t2,reals,i,j
        logical :: mobseqmtp

        f3 = f1234(3)
        f2 = f1234(2)
        mobseqmtp=.false.
        ! rhorizon = one + dsqrt(one-a_spin*a_spin)
        rhorizon = (one + lambdaBar**two) + sqrt((one + lambdaBar**two)**two - a_spin**two)
        call mutp(f2,f3,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)

        If(reals.eq.2)then
            Do j = 0,10
                Do i=j,j+1
                    If(mobseqmtp)then
                        If(muobs.eq.mu_tp1)then
                            t1=j
                            t2=i
                        else
                            t1=i
                            t2=j
                        endif
                    else
                        If(muobs.gt.zero)then
                            If(f2.lt.zero)then
                                t1=i
                                t2=j
                            endif
                            If(f2.gt.zero)then
                                t1=j
                                t2=i
                            endif
                        else
                            If(muobs.eq.zero)then
                                If(f2.lt.zero)then
                                    t1=i
                                    t2=j
                                endif
                                If(f2.gt.zero)then
                                    t1=j
                                    t2=i
                                endif
                            else
                                If(f2.lt.zero)then
                                    t1=i
                                    t2=j
                                endif
                                If(f2.gt.zero)then
                                    t1=j
                                    t2=i
                                endif
                            endif
                        endif
                    endif
                    pm=mu2p(f3,f2,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal)
                    IF(pm .le. zero)cycle
                    ! re = radius(pm,f1234(1),lambda,q,a_spin,robs,scal)
                    re = radius(pm,f1234(1),lambda,q,a_spin,lambdaBar,robs,scal)
                    !write(*,*)'mu=',re,pm,t1,t2
                    IF(re .le. rout .and. re .ge. rin)THEN
                        Pemdisk_all = pm
                        return
                    ELSE
                        IF(re .ge. infinity)THEN
                            Pemdisk_all = -one
                            RETURN
                        ELSE
                            IF(re .le. rhorizon)THEN
                                Pemdisk_all = -two
                                return
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
        else
            pm=-two
        endif
        Pemdisk_all=pm
        return
    End Function Pemdisk_all
    !*****************************************************************************************************
    subroutine metricg(robs,sinobs,muobs,a_spin,lambdaBar,somiga,expnu,exppsi,expmu1,expmu2)
        !*****************************************************************************************************
        !*     PURPOSE:  Computes Kerr metric, exp^\nu, exp^\psi, exp^mu1, exp^\mu2, and omiga at position:
        !*               r_obs, \theta_{obs}.
        !*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*     OUTPUTS:  somiga,expnu,exppsi,expmu1,expmu2------------Kerr metrics under Boyer-Lindquist coordinates.
        !*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012).
        !*     DATE WRITTEN:  5 Jan 2012.
        !*     REVISIONS: ******************************************
        implicit none
        Double precision robs,a_spin,lambdaBar,DeltaK,DeltaWH,bigA,two,sinobs,muobs,one,sigma
        Double precision ,optional :: somiga,expnu,exppsi,expmu1,expmu2


        !
        ! Only-g_rr-modified model: g_tt, g_tphi, g_phiphi use Kerr DeltaK;
        ! only g_rr uses wormhole DeltaWH.
        two=2.D0
        one=1.D0
        ! Kerr Delta
        DeltaK = robs**two - two*robs + a_spin**two
        ! Wormhole Delta
        DeltaWH = robs**two - two*robs*(one+lambdaBar**two) + a_spin**two
        sigma=robs**two+(a_spin*muobs)**two
        ! bigA uses Kerr DeltaK
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*DeltaK
        ! Frame dragging: Kerr somiga = 2ar / bigA
        somiga=two*a_spin*robs/bigA
        ! exp(nu), exp(psi) use Kerr DeltaK; exp(mu1) uses wormhole DeltaWH
        expnu=sqrt(sigma*DeltaK/bigA)
        exppsi=sinobs*sqrt(bigA/sigma)
        expmu1=sqrt(sigma/DeltaWH)
        expmu2=sqrt(sigma)
        return
    End subroutine metricg


    !********************************************************************************************
    !   Subroutine lambdaq(alpha,beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q)
    Subroutine lambdaq(alpha,beta,robs,sinobs,muobs,a_spin,lambdaBar,scal,velocity,f1234,lambda,q)

        !********************************************************************************************
        !*     PURPOSE:  Computes constants of motion from impact parameters alpha and beta by using
        !*               formulae (86) and (87) in Yang & Wang (2012).
        !*     INPUTS:   alpha,beta-----Impact parameters.
        !*               robs-----------radial coordinate of observer or the initial position of photon.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
        !*                              LNRF.
        !*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of
        !*                              four momentum of a photon measured under the LNRF frame, and
        !*                              defined by equations (82)-(85) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*     ROUTINES CALLED: NONE.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012).
        !*     DATE WRITTEN:  5 Jan 2012.
        !*     REVISIONS: ******************************************
        implicit none
        Double precision f1234(4),robs,sinobs,muobs,a_spin,lambda,q,A1,&
        lambdaBar,DeltaK,DeltaWH,zero,one,two,four,Sigma,at,Bt,&
        scal,Vr,Vt,Vp,gama,expnu2,&
        eppsi2,epmu12,epmu22,bigA,three,alpha,beta,&
        velocity(3),somiga,prt,ptt,ppt
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)

        if(abs(beta).lt.1.D-7)beta=zero
        if(abs(alpha).lt.1.D-7)alpha=zero
        ! equations (94), (95) in Yang & Wang (2012).
        at=alpha/scal/robs
        Bt=beta/scal/robs
        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
        ! equation (90) in Yang & Wang (2012).
        gama=one/dsqrt(one-(Vr**two+Vt**two+Vp**two))

        ! equations (97), (98), (99) in Yang & Wang (2012).
        prt=-one/dsqrt(one+at**two+Bt**two)
        ptt=Bt*prt
        ppt=at*prt
        ! equations (89), (90) and (91) in Yang & Wang (2012).
        f1234(1)=( gama*Vr-prt*(one+gama*gama*Vr*Vr/(one+gama))-&
        ptt*gama*gama*Vr*Vt/(one+gama)-ppt*gama*gama*Vr*Vp/(one+gama) )*robs*scal
        f1234(2)=( gama*Vt-prt*gama*gama*Vt*Vr/(one+gama)-ptt*(one+&
        gama*gama*Vt*Vt/(one+gama))-ppt*gama*gama*Vt*Vp/(one+gama) )*robs*scal
        f1234(3)=( gama*Vp-prt*gama*gama*Vp*Vr/(one+gama)-ptt*gama*gama*Vp*Vt/(one+gama)-&
        ppt*(one+gama*gama*Vp*Vp/(one+gama)) )*robs*scal
        f1234(4)=gama*(one-prt*Vr-ptt*Vt-ppt*Vp)



        !in the above equations, f1234(1), f1234(2), f1234(3), f1234(4) denote pr, ptheta, pphi, pt
        !a transformation was done from the source frame (denoted p_accent in the paper) with velocities vr, vt, vp to the LNRF frame
        !in the paper, see 5th line below eq. 89
        ! Keep r component p_r of four momentum to be negative, so the photon will go
        ! to the central black hole.
        f1234(1)=-f1234(1)
        If(dabs(f1234(1)).lt.1.D-6)f1234(1)=zero
        If(dabs(f1234(2)).lt.1.D-6)f1234(2)=zero
        If(dabs(f1234(3)).lt.1.D-6)f1234(3)=zero
        If(dabs(f1234(4)).lt.1.D-6)f1234(4)=zero
        ! Only-g_rr-modified: Kerr DeltaK for bigA/somiga/expnu2; wormhole DeltaWH for radial terms.
        DeltaK = robs**two-two*robs + a_spin**two
        DeltaWH = robs**two-two*robs*(one+lambdaBar**two) + a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*DeltaK
        somiga=two*a_spin*robs/bigA
        expnu2=Sigma*DeltaK/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/DeltaWH
        epmu22=Sigma
        ! equations (86) and (87) in Yang & Wang (2012).  Use wormhole DeltaWH in radial terms.
        A1 = f1234(3)/(dsqrt(DeltaWH)*Sigma/bigA*f1234(4)*robs*scal+f1234(3)*somiga*sinobs)
        lambda = A1*sinobs
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/f1234(4)/robs/scal*&
        (one-lambda*somiga))**two*bigA/DeltaWH
        return
    End subroutine lambdaq

    !********************************************************************************************
    Subroutine initialdirection(pr,ptheta,pphi,sinobs,&
        muobs,a_spin,robs,velocity,lambda,q,f1234)
        !********************************************************************************************
        !*     PURPOSE:  Computes constants of motion from components of initial 4 momentum
        !*               of photon measured by emitter in its local rest frame, by using
        !*               formulae (86) and (87) in Yang & Wang (2012).
        !*     INPUTS:   pr,ptheta,pphi-----components of initial 4 momentum of photon measured by
        !*               emitter in its local rest frame.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of observer or the initial position of photon.
        !*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
        !*                              LNRF.
        !*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of
        !*                              four momentum of a photon measured under the LNRF frame, and
        !*                              defined by equations (82)-(85) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*     ROUTINES CALLED: NONE.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012).
        !*     DATE WRITTEN:  5 Jan 2012.
        !*     REVISIONS: ******************************************
        implicit none
        Double precision lambda,q,sinobs,muobs,a_spin,robs,zero,one,two,three,four,&
        velocity(3),Vr,Vt,Vp,&
        a1,gama,gama_tp,gama_p,f1234(4),&
        lambdaBar,Delta,DeltaK,DeltaWH,Sigma,bigA,eppsi2,epmu12,epmu22,somiga,&
        pr,ptheta,pphi
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)
        optional pr,ptheta,pphi

        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
        ! equation (92) in Yang & Wang (2012).
        gama=one/sqrt(one-(Vr**two+Vt**two+Vp**two))
        gama_tp=one/sqrt(one-(Vt**two+Vp**two))
        gama_p=one/sqrt(one-Vp**two)
        ! equation (106)-(109) in Yang & Wang (2012).
        f1234(1)=-(-gama*Vr+gama/gama_tp*pr)
        f1234(2)=-(-gama*Vt+gama*gama_tp*Vr*Vt*pr+gama_tp/gama_p*ptheta)
        f1234(3)=-(-gama*Vp+gama*gama_tp*Vr*Vp*pr+gama_tp*gama_p*Vt*Vp*ptheta+gama_p*pphi)
        f1234(4)=(gama-gama*gama_tp*Vr*pr-gama_tp*gama_p*Vt*ptheta-gama_p*Vp*pphi)


        !*************************** added by Pieter ********************************************** codewoord kreeft
        ! equations (89), (90) and (91) in Yang & Wang (2012).
        !        f1234(1)=( gama*Vr-pr*(one+gama*gama*Vr*Vr/(one+gama))-&
        !                 pt*gama*gama*Vr*Vt/(one+gama)-pp*gama*gama*Vr*Vp/(one+gama) )
        !        f1234(2)=( gama*Vt-pr*gama*gama*Vt*Vr/(one+gama)-pt*(one+&
        !                 gama*gama*Vt*Vt/(one+gama))-pp*gama*gama*Vt*Vp/(one+gama) )
        !        f1234(3)=( gama*Vp-pr*gama*gama*Vp*Vr/(one+gama)-pt*gama*gama*Vp*Vt/(one+gama)-&
        !                 pp*(one+gama*gama*Vp*Vp/(one+gama)) )
        !        f1234(4)=gama*(one-pr*Vr-pt*Vt-pp*Vp)
        !        write(*,*) f1234(1), f1234(2), f1234(3), f1234(4)
        !*************************** End of added by Pieter ***************************************
        If(abs(f1234(1)).lt.1.D-7)f1234(1)=zero
        If(abs(f1234(2)).lt.1.D-7)f1234(2)=zero
        If(abs(f1234(3)).lt.1.D-7)f1234(3)=zero
        If(abs(f1234(4)).lt.1.D-7)f1234(4)=zero

        ! equations (1), (2) in Yang & Wang (2012).
        ! Separate Kerr and wormhole Deltas.  bigA uses Kerr DeltaK, whereas
        ! the radial metric terms use wormhole DeltaWH.
        ! Full Delta_WH model: use DeltaWH consistently.
        DeltaWH = robs**two - two*robs*(one+lambdaBar**two) + a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*DeltaWH
        somiga=two*(one+lambdaBar**two)*a_spin*robs/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/DeltaWH
        epmu22=Sigma
        ! equations (86) and (87) in Yang & Wang (2012).
        A1 = f1234(3)/(dsqrt(DeltaWH)*Sigma/bigA*f1234(4)+f1234(3)*somiga*sinobs)
        lambda = A1*sinobs
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/f1234(4)*&
        (one-lambda*somiga))**two*bigA/DeltaWH
        return
    End subroutine initialdirection

    !********************************************************************************************
    Subroutine center_of_image(robs,scal,velocity,alphac,betac)
        !********************************************************************************************
        !*     PURPOSE:  Solves equations f_3(alphac,betac)=0, f_2(alphac,betac)=0, of (100)
        !*               and (101) in Yang & Wang (2012). alphac, betac are the coordinates of
        !*               center point of images on the screen of observer.
        !*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
        !*                              LNRF.
        !*     OUTPUTS:  alphac,betac---coordinates of center point of images on the screen of observer.
        !*     ROUTINES CALLED: NONE.
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012).
        !*     DATE WRITTEN:  9 Jan 2012.
        !*     REVISIONS: ******************************************
        implicit none
        Double precision robs,scal,velocity(3),alphac,betac,zero,one,two,four,&
        Vr,Vt,Vp,gama,a1,b1,c1,alphap,alpham,betap,betam
        parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0)

        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
        ! equation (90) in Yang & Wang (2012).
        gama=one/sqrt(one-(Vr*Vr+Vt*Vt+Vp*Vp))

        If(Vt.ne.zero)then
            If(Vp.ne.zero)then
                a1=(one+gama*gama*(Vt*Vt+Vp*Vp)/(one+gama))**two-gama*gama*(Vp*Vp+Vt*Vt)
                b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*(Vt*Vt+Vp*Vp))/(one+gama)**two
                c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt
                betap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
                betam=(-b1-sqrt(b1**two-four*a1*c1))/two/a1
                If(betap*Vp.lt.zero)then
                    betac=betap
                Else
                    betac=betam
                Endif
                alphac=Vp/Vt*betac
            Else
                alphac=zero
                a1=(one+gama*gama*Vt*Vt/(one+gama))**two-gama*gama*Vt*Vt
                b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*Vt*Vt)/(one+gama)**two
                c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt
                betap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
                betam=(-b1-sqrt(b1**two-four*a1*c1))/two/a1
                If(betap*Vt.lt.zero)then
                    betac=betap
                Else
                    betac=betam
                Endif
            Endif
        else
            betac=zero
            If(Vp.ne.zero)then
                a1=(one+gama*gama*Vp*Vp/(one+gama))**two-gama*gama*Vp*Vp
                b1=two*gama*gama*Vp*Vr*(one+gama+gama*gama*Vp*Vp)/(one+gama)**two
                c1=(gama*gama*Vp*Vr/(one+gama))**two-gama*gama*Vp*Vp
                alphap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
                alpham=(-b1-sqrt(b1**two-four*a1*c1))/two/a1
                If(alphap*Vp.lt.zero)then
                    alphac=alphap
                Else
                    alphac=alpham
                Endif
            Else
                alphac=zero
            Endif
        endif
        alphac=alphac*robs*scal
        betac=betac*robs*scal
        return
    End Subroutine center_of_image


    !********************************************************************************************
    FUNCTION p_total(f1234r,lambda,q,sinobs,muobs,a_spin,robs,scal)
        !********************************************************************************************
        !*     PURPOSE:  Computes the integral value of \int^r dr (R)^{-1/2}, from the starting position to
        !*               the termination----either the infinity or the event horizon.
        !*     INPUTS:   f1234r---------p_r, the r component of four momentum of the photon measured
        !*                              under the LNRF, see equation (83) in Yang & Wang (2012).
        !*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
        !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where
        !*                              \theta_{obs} is the inclination angle of the observer.
        !*               a_spin---------spin of black hole, on interval (-1,1).
        !*               robs-----------radial coordinate of observer or the initial position of photon.
        !*               scal-----------a dimentionless parameter to control the size of the images.
        !*                              Which is usually be set to 1.D0.
        !*     OUTPUTS:  p_total--------which is the value of integrals \int^r dr (R)^{-1/2}, along a
        !*                              whole geodesic, that is from the starting position to either go to
        !*                              infinity or fall in to black hole.
        !*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5
        !*     ACCURACY:   Machine.
        !*     AUTHOR:     Yang & Wang (2012)
        !*     DATE WRITTEN:  5 Jan 2012
        !*     REVISIONS: ******************************************
        IMPLICIT NONE
        DOUBLE PRECISION sinobs,muobs,a_spin,lambdaBar,rhorizon,q,lambda,&
        cc,b0,b1,b2,b3,g2,g3,tobs,pp,p1,p2,PI0,&
        u,v,w,L1,L2,thorizon,m2,pinf,r_add,r_m,&
        a4,b4,robs,&
        scal,tinf,integ04(4),integ14(4),&
        r_tp1,r_tp2,t_inf,tp2,f1234r,&
        PI01,rff_p,p_total,p_tp1_tp2,PI2_p,PI1_p,sqt3
        !PARAMETER(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0)
        COMPLEX*16 bb(1:4),dd(3)
        INTEGER ::  reals,t1,t2,index_p4(4),del,cases_int,cases
        LOGICAL :: robs_eq_rtp,indrhorizon

        ! rhorizon=one+sqrt(one-a_spin**two)
        rhorizon=one+sqrt(one-a_spin**two)
        ! equation (64) in Yang & Wang (2012).
        r_add=rhorizon
        ! r_m=one-sqrt(one-a_spin**two)
        r_m=one-sqrt(one-a_spin**two)

        ! equation (64) in Yang & Wang (2012).
        b4=one
        a4=zero
        cc=a_spin**2-lambda**2-q
        robs_eq_rtp=.false.
        indrhorizon=.false.
        ! call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
        !                   reals,robs_eq_rtp,indrhorizon,cases,bb)
        call radiustp(f1234r,a_spin,lambdaBar,robs,lambda,q,r_tp1,r_tp2,&
        reals,robs_eq_rtp,indrhorizon,cases,bb)
        !** R(r)=0 has real roots and turning points exists in radial r.
        If(reals.ne.0)then
            ! equations (35)-(38) in Yang & Wang (2012).
            ! Use pure Kerr quartic coefficients (turning points are Kerr).
            b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
            b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
            b2=four/three*r_tp1
            b3=one
            g2=three/four*(b1**2-b0*b2)
            g3=one/16.D0*(three*b0*b1*b2-two*b1**three-b0**two*b3)
            ! equation (39) in Yang & Wang (2012).
            If(robs-r_tp1.ne.zero)then
                tobs=b0/four/(robs-r_tp1)+b1/four
            else
                tobs=infinity
            endif
            If(rhorizon-r_tp1.ne.zero)then
                thorizon=b1/four+b0/four/(rhorizon-r_tp1)
            else
                thorizon=infinity
            endif
            tp2=b0/four/(r_tp2-r_tp1)+b1/four
            tinf=b1/four
            ! equation (64), (66) and (70) in Yang & Wang (2012).
            call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)

            index_p4(1)=0
            cases_int=1
            call weierstrass_int_J3(tobs,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
            ! equation (42) in Yang & Wang (2012).
            PI0=integ04(1)
            select case(cases)
            CASE(1)
                If(f1234r .ge. zero)then !**photon will goto infinity.
                    index_p4(1)=0
                    cases_int=1
                    call weierstrass_int_J3(tinf,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                    p_total = integ04(1)
                ELSE
                    If(.not.indrhorizon)then
                        index_p4(1)=0
                        cases_int=1
                        call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                        p_total = PI0+integ14(1)
                    ELSE     !f1234r<0, photon will fall into black hole unless something encountered.
                        index_p4(1)=0
                        cases_int=1
                        call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        p_total = integ04(1)
                    ENDIF
                ENDIF
            CASE(2)
                If(.not.indrhorizon)then
                    write(*,*)'we come here!'
                    If(f1234r.lt.zero)then
                        PI01=-PI0
                    else
                        PI01=PI0
                    endif
                    ! equation (41) in Yang & Wang (2012).
                    index_p4(1)=0
                    cases_int=1
                    call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                    pp=zero
                    ! equation (57) in Yang & Wang (2012).
                    p_tp1_tp2=integ14(1)
                    PI2_p=p_tp1_tp2-PI0
                    PI1_p=PI0
                    p1=PI0-pp
                    p2=p_tp1_tp2-p1
                    ! equation (58) in Yang & Wang (2012).
                    t1 = 2
                    t2 = 2
                    If(robs_eq_rtp)then
                        p_total=abs(pp)+two*(t1*p1+t2*p2)
                    else
                        If(f1234r.gt.zero)then
                            p_total=-pp+two*(t1*p1+t2*p2)
                        endif
                        If(f1234r.lt.zero)then
                            p_total=pp+two*(t1*p1+t2*p2)
                        endif
                    endif
                    !*************************************************************************************
                    200     continue
                else  !photon has probability to fall into black hole.
                    If(f1234r.le.zero)then
                        index_p4(1)=0
                        cases_int=1
                        call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        p_total = integ04(1)
                    ELSE  !p_r>0, photon will meet the r_tp2 turning point and turn around then goto vevnt horizon.
                        index_p4(1)=0
                        cases_int=1
                        call weierstrass_int_J3(tp2,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        call weierstrass_int_j3(tp2,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                        p_total = integ14(1)+integ04(1)
                    ENDIF
                ENDIF
            END SELECT
            !************************************************************************************************
            If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        p_total = one/rhorizon-one/robs
                    else
                        p_total = one/robs
                    endif
                endif
                If(cc.eq.-27.D0)then
                    sqt3=sqrt(three)
                    If(f1234r.lt.zero)then
                        p_total=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                        -Log(one+two/sqt3)/three/sqt3
                    else
                        p_total=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                        -Log(one+two/sqt3)/three/sqt3
                    endif
                endif
            endif
        ELSE
            ! equation (44) in Yang & Wang (2012).   !equation R(r)=0 has no real roots. we use the Legendre elliptic
            u=real(bb(4))        !integrations and functions to compute the calculations.
            w=abs(aimag(bb(4)))
            v=abs(aimag(bb(2)))
            If(u.ne.zero)then
                ! equation (45) in Yang & Wang (2012).
                L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                ! equation (46) in Yang & Wang (2012).
                thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
                ! equation (48) in Yang & Wang (2012).
                m2=(L1-L2)/L1
                tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                t_inf=sqrt((L1-one)/(L1-L2))
                ! equation (50) in Yang & Wang (2012).
                pinf=EllipticF(tinf,m2)/(w*sqrt(L1))
                IF(f1234r.lt.zero)THEN
                    p_total = pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))
                ELSE
                    p_total = EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                ENDIF
            ELSE
                If(f1234r.lt.zero)then
                    p_total = (atan(robs/w)-atan(rhorizon/w))/w
                ELSE
                    p_total = (PI/two-atan(robs/w))/w
                ENDIF
            ENDIF
        ENDIF
        RETURN
    END FUNCTION p_total

!********************************************************************************************
end module BLcoordinate
