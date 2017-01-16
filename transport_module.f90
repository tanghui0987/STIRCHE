    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This part is used to modify the fluid density.                                     !
        ! Based on Warner, 2006                                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine density(mbc,mx,my,c,rhom)

            use sediment_module, only: rhos,rho,gmax
            use Set_Precision, only: Prec

            implicit none

            ! inputs
            integer, intent(in) :: mbc, mx, my
            Real(kind=Prec), intent(in) :: c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            ! local
            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::  ctotal

            !output
            Real(kind=Prec),intent(out),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: rhom

            ctotal = 0.0
            rhom = 0.0

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax
                        ctotal(i,j) = ctotal(i,j) + c(i,j,k)
                    enddo
                enddo
            enddo

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    rhom(i,j) = rho + c(i,j)*(rhos-rho)
                enddo
            enddo

        end subroutine density

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This part is used to calculate bed roughness                                       !
        ! based on Wieberg 1989                                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine bedrough(mbc,mx,my,u,v,h,c,ub_cr,z0)

            use sediment_module, only: hcr,D,g,m0,rhos,rho,gmax,pbbed
            use Set_Precision, only: Prec

            implicit none

            ! argument
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc),&
            v(1-mbc:mx+mbc,1-mbc:my+mbc), c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax), ub_cr(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local

            integer :: i,j
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,hloc,Dmm,zon,a2,ustarc,ustarcrit,taub,taucrit,Tstar,delb,zos,rhom
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,ub_cr,us_cr1,us_cr2
            Real(kind=Prec) :: a1,gammaWs

            !output
            Real(kind=Prec),intent(out),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: z0

            gammaWs = 0.056

            vmag2 = u**2.0+v**2.0
            hloc = max(h,hcr)

            call density(mbc,mx,my,c,rhom)

            !call crtical_velocity1(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)
            !call crtical_velocity2(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)
            !call crtical_velocity3(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

            Dmm = meansize(D,pbbed(:,:,1,:),mx,my,mbc,gmax)
            zon = Dmm/30.0
            a1 = 0.68
            a2 = 0.0204*(log(Dmm*1000))**2.0+0.0220*log(Dmm*1000)+0.0709
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    ustarc(i,j) = sqrt(g*m0**2.0*vmag2(i,j)/((rhom(i,j)-rho)*Dmm(i,j)*hloc(i,j)**(1.0/3.0)))
                    ustarcrit(i,j) = sqrt(g*m0**2.0*ub_cr(i,j,int(gmax/2.0+1.0))**2.0/((rhom(i,j)-rho)*Dmm(i,j)*hloc(i,j)**(1.0/3.0)))
                enddo
            enddo
            taub = rho*ustarc**2.0
            taucrit=rho*ustarcrit**2.0
            do i = 1-mbc,mx+mbc
                do j =1-mbc, my+mbc
                    Tstar(i,j)=taub(i,j)/taucrit(i,j)
                    delb(i,j)=Dmm(i,j)*a1*Tstar(i,j)/(1.0+a2(i,j)*Tstar(i,j))
                enddo
            enddo
            zos = gammaWs*delb
            z0 = zon+zos
        end subroutine bedrough
        !*******************************************************************
        !
        !This function is used to caculate mean grain size
        !
        !*******************************************************************
        function meansize(D,fr,mx,my,mbc,gmax)

            implicit none

            integer :: mx,my,mbc,gmax

            Real(kind=Prec), Dimension(gmax) :: D
            Real(kind=Prec), Dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: fr
            Real(kind=Prec), Dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: meansize
            integer :: i,j,k

            meansize = 0.0

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax
                        meansize(i,j) = meansize(i,j)+D(k)*fr(i,j,k)
                    enddo
                enddo
            enddo
        end function meansize


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This part is used to estimate settling velocity.                                   !
        ! Ahrens (2000)                                                                      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine settling_velocity(mbc,mx,my,c,ws)

            use sediment_module, only: rhos,rho,gmax,g,D
            use Set_Precision, only: Prec

            !use params

            implicit none

            ! argument
            integer, intent(in) :: mbc,mx,my
            Real(kind=Prec), intent(in) :: c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local
            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::  ctotal,rhom,delta
            Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Sster,c1,c2,wster,w,R,alpha
            Real(kind=Prec) :: Te,vis

            !output
            Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ws

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax
                        ctotal(i,j) = ctotal(i,j) + c(i,j,k)
                    enddo
                enddo
            enddo

            call density(mbc,mx,my,c,rhom)

            delta  = (rhom-rho)/rho

            do i =1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax
                        Te    = 20.0
                        vis  = 0.1*1e-5 ! Van rijn, 1993
                        Sster(i,j,k) = D(k)/(4*vis)*sqrt(delta(i,j)*g*D(k))
                        c1(i,j,k)    = 1.06*tanh(0.064*Sster(i,j,k)*exp(-7.5/Sster(i,j,k)**2.0))
                        c2(i,j,k)    = 0.220*tanh(2.34*Sster(i,j,k)**(-1.180)*exp(-0.0064*Sster(i,j,k)**2.0))
                        wster(i,j,k) = c1(i,j,k)+c2(i,j,k)*Sster(i,j,k)
                        w(i,j,k) = wster(i,j,k)*sqrt(delta(i,j)*g*D(k))
                        !reduce settling velocity
                        R(i,j,k) = w(i,j,k)*D(k)/vis
                        alpha(i,j,k) = 2.35*(2.0+0.175*R(i,j,k)**(3.0/4.0))/(1.0+0.175*R(i,j,k)**(3.0/4.0))
                        ws(i,j,k) = (1-ctotal(i,j))**alpha(i,j,k)*w(i,j,k)
                        ! for gravel transport
                        ! dster(i,j,k)=(delta(i,j)*g/1.d-12)**(1.0/3.0)*D(k)
                    enddo
                enddo
            enddo

        end subroutine settling_velocity



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This Part is used to calculate critical velocity for each grain size classes!
        ! Van (1984)                                                                  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine crtical_velocity1(mbc,mx,my,h,u,v,c,ub_cr,us_cr1,us_cr2)

            use sediment_module, only: rhos,rho,gmax,g,D,hcr,k0,m0,pbbed
            use Set_Precision, only: Prec

            implicit none

            ! argument
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc) &
                            ,v(1-mbc:mx+mbc,1-mbc:my+mbc),c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local
            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: z0, rhom, delta, hloc
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ws,ub_crt, us_cr1t, us_cr2t, us_cr1s, us_cr2s

            !output
            Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ub_cr, us_cr1, us_cr2

            call settling_velocity(mbc,mx,my,c,ws)

            call density(mbc,mx,my,c,rhom)

            delta  = (rhom-rho)/rho

            hloc = max(h,hcr)

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax

                        ! only sand for this method
                        if (D(k) <= 0.0005) then
                            ub_crt(i,j,k) = 0.19*D(k)**0.1*dlog10(4*hloc(i,j)/D(k))
                        elseif (D(k) > 0.0005 .AND. D(k) <= 0.002) then
                            ub_crt(i,j,k) = 8.5*D(k)**0.6*dlog10(4.0*hloc(i,j)/D(k))
                        else
                            print *, "The Grain size is out of range"
                        endif

                        ! Use Rouse number

                        us_cr1s(i,j,k) = ws(i,j,k)/(2.5*k0)
                        us_cr2s(i,j,k) = ws(i,j,k)/(1.2*k0)

                        ! bed roughness

                        call bedrough(mbc,mx,my,u,v,h,ub_crt,z0) ! Rough estimation

                        ! Use the law of the wall

                        us_cr1t(i,j,k) = us_cr1s(i,j,k)/k0*(log(h(i,j)/k0-(1-z0(i,j)/h(i,j))))
                        us_cr2t(i,j,k) = us_cr2s(i,j,k)/k0*(log(h(i,j)/k0-(1-z0(i,j)/h(i,j))))

                        ! correction

                        ub_cr(i,j,k) = min(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr2(i,j,k) = max(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr1(i,j,k) =ub_crt(i,j,k)+us_cr1t(i,j,k)+us_cr2t(i,j,k)-us_cr2(i,j,k)-ub_cr(i,j,k)
                    enddo
                enddo
            enddo
        end subroutine crtical_velocity1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This Part is used to calculate critical velocity for each grain size classes!
        ! Van rijn, 1993                                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine crtical_velocity2(mbc,mx,my,h,u,v,c,ub_cr,us_cr1,us_cr2)

            use sediment_module, only: rhos,rho,gmax,g,D,hcr,Trep,k0,m0,pbbed
            use Set_Precision, only: Prec

            implicit none

            ! argument
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc) &
                                    ,v(1-mbc:mx+mbc,1-mbc:my+mbc),c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local
            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: z0, rhom, delta, hloc
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ws,ub_crt, us_cr1t, us_cr2t, us_cr1s, us_cr2s

            !output
            Real(kind=Prec),intent(inout) :: ub_cr(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
                                        ,us_cr1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
                                        ,us_cr2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            call settling_velocity(mbc,mx,my,c,ws)

            call density(mbc,mx,my,c,rhom)

            delta  = (rhom-rho)/rho

            hloc = max(h,hcr)

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax

                        if (D(k) <= 0.0005) then
                            ub_crt(i,j,k) = 0.19*D(k)**0.1*dlog10(4*hloc(i,j)/D(k))
                        else
                            ub_crt(i,j,k) = 1.3*sqrt(delta(i,j)*D(k))*(hloc(i,j)/D(k))**(1.0/6.0)
                        end if

                        ! use Rouse number
                        us_cr1s(i,j,k) = ws(i,j,k)/(2.5*k0)
                        us_cr2s(i,j,k) = ws(i,j,k)/(1.2*k0)

                        !bed roughness

                        call bedrough(mbc,mx,my,u,v,h,z0)

                        ! the law of the wall
                        us_cr1t(i,j,k) = us_cr1s(i,j,k)/k0*(log(h(i,j)/k0-(1-z0(i,j)/h(i,j))))
                        us_cr2t(i,j,k) = us_cr2s(i,j,k)/k0*(log(h(i,j)/k0-(1-z0(i,j)/h(i,j))))

                        !correction
                        ub_cr(i,j,k) = min(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr2(i,j,k) = max(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr1(i,j,k) =ub_crt(i,j,k)+us_cr1t(i,j,k)+us_cr2t(i,j,k)-us_cr2(i,j,k)-ub_cr(i,j,k)

                    end do
                end do
            end do
        end subroutine crtical_velocity2


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This Part is used to calculate critical velocity for each grain size classes!
        ! Shield
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine crtical_velocity3(mbc,mx,my,h,u,v,c,ub_cr,us_cr1,us_cr2)

            use sediment_module, only: rhos,rho,gmax,g,D,hcr,Trep,k0,m0,g,NL,pbbed
            use Set_Precision, only: Prec

            implicit none

            ! argument
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc) &
                                ,v(1-mbc:mx+mbc,1-mbc:my+mbc),c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local
            integer :: i,j,k,l
            Real(kind=Prec) :: vis, uscr, xr, yr
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: z0, rhom, gamma, Dmm
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ws, ub_crs, us_cr1s, us_cr2s &
                                                        ub_crt, us_cr1t, us_cr2t, ub_cr, us_cr2, us_cr1
            Real(kind=Prec),dimension(NL) :: z,K_b,K_s1,K_s2,spd_b,spd_s1,spd_s2
            Real(kind=Prec),dimension(NL-1) :: spd_b_profile, spd_s1_profile, spd_s2_profile

            !output
            Real(kind=Prec),intent(inout) :: ub_cr(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
                                            ,us_cr1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
                                            ,us_cr2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)


            call settling_velocity(mbc,mx,my,c,ws)

            call density(mbc,mx,my,c,rhom)

            gamma=g*(rhom-rho)

            vis  = 0.1*1e-5

            ! bed roughness
            Dmm = meansize(D,pbbed(:,:,1,:),mx,my,mbc,gmax)
            z0 = Dmm/30.0

            ! shield diagram
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax
                        do l = 1, 20
                            uscr=0.00001
                            xr=uscr*D(k)/vis
                            yr=0.15*xr**(-1.0)+0.05*exp(-8.0*xr**(-0.9))
                            uscr=sqrt(D(k)*gamma(i,j)*yr/rho)
                        end do
                        ub_crs(i,j,k) = uscr
                    enddo
                enddo
            enddo

            !
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax
                        ! use rouse number
                        us_cr1s(i,j,k) = ws(i,j,k)/(2.5*k0)
                        us_cr2s(i,j,k) = ws(i,j,k)/(1.2*k0)

                        ! use eddy viscosity
                        do l= 1, NL
                            z(l) = z0(i,j) + l*h(i,j)/NL
                            K_b(l)=k0*ub_crs(i,j,k)*z(l)*exp((-z(l)/h(i,j))-3.2*(z(l)/h(i,j))**2.0+2.13*(z(l)/h(i,j))**3.0)
                            K_s1(l)=k0*us_cr1s(i,j,k)*z(l)*exp((-z(l)/h(i,j))-3.2*(z(l)/h(i,j))**2+2.13*(z(l)/h(i,j))**3.0)
                            K_s2(l)=k0*us_cr2s(i,j,k)*z(l)*exp((-z(l)/h(i,j))-3.2*(z(l)/h(i,j))**2+2.13*(z(l)/h(i,j))**3.0)
                            spd_b(l)=h(i,j)/NL*ub_crs(i,j,k)**2.0/K_b(l)
                            spd_s1(l)=h(i,j)/NL*ub_cr1s(i,j,k)**2.0/K_s1(l)
                            spd_s2(l)=h(i,j)/NL*ub_cr2s(i,j,k)**2.0/K_s2(l)
                        enddo

                        do l=1, NL-1
                            spd_b_profile(l) = (spd_b(l)+spd_b(l+1))/2.0
                            spd_s1_profile(l) = (spd_s1(l)+spd_s1(l+1))/2.0
                            spd_s2_profile(l) = (spd_s2(l)+spd_s2(l+1))/2.0
                        enddo

                        ub_crt(i,j,k) = sum(spd_b_profile)/(NL-1.0)
                        us_cr1t(i,j,k) = sum(spd_s1_profile)/(NL-1.0)
                        us_cr2t(i,j,k) = sum(spd_s2_profile)/(NL-1.0)

                        ! correction
                        ub_cr(i,j,k) = min(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr2(i,j,k) = max(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr1(i,j,k) =ub_crt(i,j,k)+us_cr1t(i,j,k)+us_cr2t(i,j,k)-us_cr2(i,j,k)-ub_cr(i,j,k)

                    end do
                end do
            end do
        end subroutine critical_velocity3

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Soulsby-VanRijn Method                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine sb_vr(mbc,mx,my,u,v,h,c,Tsg,ceqbg,ceqsg)

            use sediment_module, only: eps,rhos,rho,gmax,g,D,hcr,tsfac,Tsmin,cmax,sws
            use Set_Precision, only: Prec
            use Set_variable, only:pbbed

            implicit none

            ! Arguments
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc),c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local
            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: wet,rhom,delta,urms,vmg,urms2,hloc,z0,Asb,Ass,term1,Cd
            Real(kind=Prec) :: vis,perc
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: dster,ws,ub_cr,us_cr1,us_cr2,Ts,term2,ceqb,ceqs,ceq

            !output
            Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Tsg,ceqbg,ceqsg


            wet = 0.0

            vis = 1e-6

            ! density correction

            call density(mbc,mx,my,c,rhom)

            call settling_velocity(mbc,mx,my,c,ws) !settling velocity

            ! calculate threshold velocity

            call crtical_velocity1(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

            !call crtical_velocity3(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

            delta = rhom - rho

            hloc = max(h,hcr)

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    if (h(i,j)>eps) then
                        wet(i,j) = 1.0
                    endif
                    do k=1,gmax

                        ! Dimensionless grain size
                        dster(i,j,k)=(delta(i,j)*g/vis**2)**(1.0/3.0)*D(k)

                        !response time
                        Ts(i,j,k) = tsfac*hloc(i,j)/ws(i,j,k)
                        Tsg(i,j,k) = max(Ts(i,j,k),Tsmin)
                    enddo
                enddo
            enddo


            urms = 0.0 !assuming the root mean sqaure velocity is zero

            vmg  = dsqrt(u**2+v**2) !velocity magnitude

            urms2  = urms**2.0


            ! bed roughness

            call bedrough(mbc,mx,my,u,v,h,c,ub_cr,z0)

            do k = 1, gmax
                ! drag coefficient
                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        Cd(i,j)=(0.40/(log(max(hloc(i,j),10.0*z0(i,j))/z0(i,j))-1.0))**2.0
                    enddo
                enddo

                ! transport parameters
                Asb=0.005*hloc*(D(k)/hloc/(delta*g))**1.20/D(k)**0.75  ! bed load coefficent
                Ass=0.012*D(k)*dster(:,:,k)**(-0.60)/(delta*g*D(k))**1.20  ! suspended load coeffient

                term1 = (vmg**2.0+0.018/Cd*sws*urms2)
                term1 = sqrt(term1)
                term2 = 0.0

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc
                        if(term1(i,j)>ub_cr(i,j,k) .and. term1(i,j)<us_cr1(i,j,k) .and. hloc(i,j)>eps) then
                            term2(i,j,k)=(term1(i,j)-ub_cr(i,j,k))**2.40
                            ceqb(i,j,k) = min(Asb(i,j)*term2(i,j)/hloc(i,j),cmax/gmax)
                            ceqs(i,j,k) = 0.0
                            ceq(i,j,k) = ceqb(i,j,k)
                        endif
                        if(term1(i,j)>us_cr2(i,j,k) .and. hloc(i,j)>eps) then
                            term2(i,j,k)=(term1(i,j)-us_cr1(i,j,k))**2.40
                            ceqs(i,j,k) = min(Ass(i,j)*term2(i,j)/hloc(i,j),cmax/gmax)
                            ceqb(i,j,k) = 0.0
                            ceq(i,j,k) = ceqs(i,j,k)
                        endif
                        if (term1(i,j)>ub_cr1(i,j,k) .and. term1(i,j)<us_cr2(i,j,k) .and. hloc(i,j)>eps) then
                            term2(i,j,k)=(term1(i,j)-us_cr1(i,j,k))**2.40
                            ceq(i,j,k) = (Ass(i,j)+Asb(i,j))*term2(i,j)/hloc(i,j)
                            perc = term1(i,j)/Us_cr2(i,j,k)
                            ceqb(i,j,k) = perc*min(ceq(i,j,k),cmax/gmax/2.0)
                            ceqs(i,j,k) = (1-perc)*min(ceq(i,j,k),cmax/gmax/2.0)
                            ! or just directly calculate
                            !ceqb(i,j,k) = min(Asb(i,j)*term2(i,j)/hloc(i,j),cmax/gmax/2.0)
                            !ceqs(i,j,k) = min(Ass(i,j)*term2(i,j)/hloc(i,j),cmax/gmax/2.0)
                            !ceq(i,j,k) = ceqb(i,j,k)+ceqs(i,j,k)
                        endif
                        ceqbg(i,j,k) = ceqb(i,j,k)*wet(i,j) ! make sure only in water
                        ceqsg(i,j,k) = ceqs(i,j,k)*wet(i,j)
                    enddo
                enddo
            enddo

        end subroutine sb_vr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Van Thiel-Van Rijn Method                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine vt_vr(mbc,mx,my,u,v,h,c,Tsg,ceqbg,ceqsg)

            use sediment_module, only: eps,rhos,rho,gmax,g,D,hcr,tsfac,Tsmin,cmax,sws
            use Set_Precision, only: Prec
            use Set_variable, only:pbbed

            implicit none

            ! Arguments
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc),c(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local
            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: wet,rhom,delta,urms,vmg,urms2,hloc,Asb,Ass,term1,term2,term3
            Real(kind=Prec) :: vis
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ws,ub_cr,us_cr1,us_cr2,dster,ceq,ceqs,ceqb,Ts

            !output
            Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Tsg,ceqsg,ceqbg

            wet = 0.0

            vis = 1e-6

            ! correct density

            call density(mbc,mx,my,c,rhom)

            call settling_velocity(mbc,mx,my,c,ws) !settling velocity

            ! calculate threshold velocity

            call crtical_velocity2(mbc,mx,my,h,u,v,c,ub_cr,us_cr1,us_cr2)

            !call crtical_velocity3(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

            delta = rhom - rho

            hloc = max(h,hcr)

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc

                    ! mark the water cover area
                    if (h(i,j)>eps) then
                        wet(i,j) = 1.0
                    endif

                    do k=1,gmax
                        dster(i,j,k)=(delta(i,j)*g/vis**2.0)**(1.0/3.0)*D(k)
                        Ts(i,j,k) = tsfac*hloc(i,j)/ws(i,j,k)
                        Tsg(i,j,k) = max(Ts(i,j,k),Tsmin)
                    enddo
                enddo
            enddo

            urms = 0.0

            vmg  = dsqrt(u**2+v**2)

            urms2  = urms**2.0

            do k = 1,gmax

                ! transport parameters

                Asb=0.015*hloc*(D(k)/hloc)**1.20/(delta*g*D(k))**0.75        !bed load coefficent
                Ass=0.012*D(k)*dster(:,:,k)**(-0.60)/(delta*g*D(k))**1.20        !suspended load coeffient

                term1=vmg**2+0.640*sws*urms2

                term1=sqrt(term1)

                !initialize ceqs

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc
                        if(term1(i,j)>ub_cr(i,j,k) .and. h(i,j)>eps .and. term1(i,j)< us_cr1(i,j,k)) then
                            term2(i,j)=(term1(i,j)-Ub_cr(i,j,k))**1.50
                            term3(i,j)=0.0
                            ceq(i,j,k) = (Asb(i,j)*term2(i,j)+Ass(i,j)*term3(i,j))/hloc(i,j)
                            ceqb(i,j,k) = min(ceq(i,j,k),cmax/gmax)
                            ceqs(i,j,k) = 0.0
                        endif
                        if(term1(i,j)>us_cr2(i,j,k) .and. h(i,j)>eps) then
                            term2(i,j)=0.0
                            term3(i,j)=(term1(i,j)-us_cr1(i,j,k))**2.4
                            ceq(i,j,k) = (Asb(i,j)*term2(i,j)+Ass(i,j)*term3(i,j))/hloc(i,j)
                            ceqb(i,j,k) = 0.0
                            ceqs(i,j,k) = min(ceq(i,j,k),cmax/gmax)
                        end if
                        if(term1(i,j)>us_cr1(i,j,k) .and. h(i,j)>eps .and. term1(i,j)< us_cr2(i,j,k)) then
                            term2(i,j)=(term1(i,j)-Ub_cr(i,j,k))**1.50
                            term3(i,j)=(term1(i,j)-Us_cr1(i,j,k))**2.40
                            ceqb(i,j,k) = min(Asb(i,j)*term2(i,j)/hloc(i,j),cmax/gmax/2.0)
                            ceqs(i,j,k) = min(Ass(i,j)*term3(i,j)/hloc(i,j),cmax/gmax/2.0)
                        end if
                        ! make sure only in water
                        ceqbg(i,j,k) = ceqb(i,j,k)*wet(i,j)
                        ceqsg(i,j,k) = ceqs(i,j,k)*wet(i,j)
                    enddo
                enddo
            enddo
        end subroutine vt_vr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Part is used to calculate Source term for finite volume method                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine transus(mbc,mx,my,dx,dy,time,u,v,h,dt,cut,cubt,cvt,cvbt,ccgt,ccbgt,Susgt,Subgt,Svbgt,Svsgt,ero1,ero2,depo_ex1,depo_ex2)

            use flux, only: Flux_vector
            use sediment_module, only: trim,gmax,morfac,por,D,thetanum,cmax,lmax,eps,facDC,nuh,nuhfac,rho,thick,method
            use Set_Precision, only: Prec
            use Set_variable,only: zb

            implicit none

            ! Arguments
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: dx,dy,dt,time
            real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)
            !local

            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,hold,DR,Dc,wet,dzbdx, dzbdy, &
                                                                        pbbedu,pbbedv,dzbdy!hold?
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,fac,ero1,ero2,depo_ex1,depo_ex2, &
                                                                        cc,ccb,dcsdy,dcbdy,dcsdx,dcbdx，Tsg
            Real(kind=Prec) :: exp_ero,facsl

            !out
            Real(kind=Prec),intent(inout) :: cut(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cubt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: cvt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cvbt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: ccgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),ccbgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: Susgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Subgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: Svbgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Svsgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: ero1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),ero2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: depo_ex1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),depo_ex2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            vmag2     = u**2+v**2
            cut = 0.0
            cubt = 0.0
            cvt = 0.0
            cvbt = 0.0

            dcsdx = 0.0
            dcbdx = 0.0
            dcsdy = 0.0
            dcbdy = 0.0

            DR = 0.0 !turn off roller model
            wet = 0.0
            exp_ero = 0.0
            facsl = 1.6

            ! calculate equibrium sediment concentration
            if (trim=='soulsby_vanrijn') then           ! Soulsby van Rijn
                call sb_vr(mbc,mx,my,u,v,h,c,Tsg,ceqbg,ceqsg)
            elseif (trim=='vanthiel_vanrijn') then       ! Van Thiel de Vries & Reniers 2008
                call vt_vr(mbc,mx,my,u,v,h,c,Tsg,ceqbg,ceqsg)
            end if

            ! mark water cover area
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    if (h(i,j)>eps) then
                        wet(i,j) = 1.0
                    endif
                enddo
            enddo

            ! x-slope
            do i = 1-mbc, mx+mbc-1
                do j = 1-mbc, my+mbc
                    dzbdx(i,j) = (zb(i+1,j)-zb(i,j))/dx
                enddo
            enddo

            dzbdx(mx+mbc,:) = dzbdx(mx+mbc-1,:)

            ! y-slope
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc-1
                    dzbdy(i,j) = (zb(i,j+1)-zb(i,j))/dy
                enddo
            enddo

            dzbdy(:,my+mbc) = dzbdy(:,my+mbc-1)

            ! compute reduction factor for sediment sources due to presence of hard layers
            frc = pbbed(:,:,1,:)

            do k = 1,gmax
                do j= 1-mbc,my+mbc
                    do i= 1-mbc,mx+mbc
                        exp_ero = morfac*dt/(1.0-por)*h(i,j)*(ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k) &
                                + ceqbg(i,j,k)*frc(i,j,k)/dt) !ceqsg, ceqbg from vt_vr or sb_vr
                        fac(i,j,k) =min(1.0,thick*frc(i,j,k)/max(tiny(0.0),exp_ero))! limit to first layer
                    enddo
                enddo
            enddo

            ! compute diffusion coefficient, usually just based on horizontal viscosity
            Dc = facDc*(nuh+nuhfac*h*(DR/rho)**(1.0/3.0))
            cc = ccgt
            ccb = ccbgt
            !Sus = 0.0
            !Sub = 0.0
            do k = 1,gmax
                if (D(k)>0.002) then
                    print *, "WARNING: Grain size is larger than 2 mm"
                    ! set ceqsg to zero for gravel; in the future, we will add this module
                endif
                ! x-direction
                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc-1
                        if(u(i,j)>0.0) then
                            cut(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(min(i+1,mx+mbc-1),j,k)
                            cubt(i,j,k)=thetanum*frc(i,j,k)*ceqbg(i,j,k)+(1.0-thetanum)&
                                *frc(min(i+1,mx+mbc-1),j,k)*ceqbg(min(i+1,mx+mbc-1),j,k)
                        elseif(u(i,j)<0.0) then
                            cut(i,j,k)=thetanum*cc(i+1,j,k)+(1.0-thetanum)*cc(max(i,2-mbc),j,k)
                            cubt(i,j,k)=thetanum*frc(i+1,j,k)*ceqbg(i+1,j,k)+(1.0-thetanum)&
                                *frc(max(i,2-mbc),j,k)*ceqbg(max(i,2-mbc),j,k)
                            !cubt(i,j,k)=thetanum*ccb(i+1,j,k)+(1.0-thetanum)*ccb(max(i,2-mbc),j,k)
                        else
                            cut(i,j,k)=0.50*(cc(i,j,k)+cc(i+1,j,k))
                            cubt(i,j,k)=0.50*(frc(i,j,k)*ceqbg(i,j,k)+frc(i+1,j,k)*ceqbg(i+1,j,k))
                            !cubt(i,j,k)=0.50*(ccb(i,j,k)+ccb(i+1,j,k))
                        endif
                        dcsdx(i,j,:)=(cc(i+1,j,:)-cc(i,j,:))/dx
                        dcbdx(i,j,:)=(ccb(i+1,j,:)-ccb(i,j,:))/dx
                    enddo
                enddo
                cut(mx+mbc,:,:) = cc(mx+mbc,:,:)
                cubt(mx+mbc,:,:) = ccb(mx+mbc,:,:)

                Sus = 0.0
                Sub = 0.0


                ! maybe move to the end
                if (method=='Van_leer') then           ! Van leer
                    call Flux_vector(mbc,mx,my,u,v,h,cc,ccb,Sus,Sub,Svs,Svb)
                elseif (method=='normal') then       ! Van Thiel de Vries & Reniers 2008
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            ! suspended load
                            Sus(i,j,k)=(cut(i,j,k)*u(i,j)*h(i,j)-Dc(i,j)*h(i,j)*dcsdx(i,j,k) &
                                    -facsl*cut(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdx(i,j,k))*wet(i,j)   !No bed slope term in suspended transport?
                            ! bed load
                            Sub(i,j,k)=(cubt(i,j,k)*u(i,j)*h(i,j)-facsl*cubt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdx(i,j,k))*wet(i,j)
                        enddo
                    enddo
                end if

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc-1
                        if(Sub(i,j,k)>0.0) then
                            pbbedu(i,j) = pbbed(i,j,1,k)
                        elseif(Sub(i,j,k)<0.0) then
                            pbbedu(i,j)= pbbed(i+1,j,1,k)
                        else
                            pbbedu(i,j)=0.50*(pbbed(i,j,1,k)+pbbed(i+1,j,1,k))
                        endif
                    enddo
                enddo

                pbbedu(mx+mbc,:) = pbbedu(mx+mbc-1,:)

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc
                        Sub(i,j,k) = pbbedu(i,j)*Sub(i,j,k) ! bed load are limited to the first sediment layer
                    enddo
                enddo

                !y-direction
                if (my>0) then
                    do j=1-mbc,my+mbc-1
                        do i=1-mbc,mx+mbc
                            if(v(i,j)>0) then
                                cvt(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(i,min(j+1,my+mbc),k)
                                cvbt(i,j,k)=thetanum*frc(i,j,k)*ceqbg(i,j,k)+(1.0-thetanum)&
                                            *frc(i,min(j+1,my+mbc),k)*ceqbg(i,min(j+1,my+mbc),k)
                                !cvbt(i,j,k)=thetanum*ccb(i,j,k)+(1.0-thetanum)*ccb(i,min(j+1,my+mbc),k)
                            elseif(v(i,j)<0) then
                                cvt(i,j,k)=thetanum*cc(i,j+1,k)+(1.0-thetanum)*cc(i,max(j,2-mbc),k)
                                cvbt(i,j,k)=thetanum*frc(i,j+1,k)*ceqbg(i,j+1,k)+(1.0-thetanum)&
                                            *frc(i,max(j,2),k)*ceqbg(i,max(j,2),k)
                                !cvbt(i,j,k)=thetanum*ccb(i,j+1,k)+(1.0-thetanum)*ccb(i,max(j,2-mbc),k)
                            else
                                cvt(i,j,k)=0.50*(cc(i,j,k)+cc(i,j+1,k)) !Jaap: cc instead of cv
                                cvbt(i,j,k)=0.50*(frc(i,j,k)*ceqbg(i,j,k)+frc(i,j+1,k)*ceqbg(i,j+1,k))
                                !cvbt(i,j,k)=0.50*(ccb(i,j,k)+ccb(i,j+1,k))
                            end if
                            dcsdy(i,j,:)=(cc(i,j+1,:)-cc(i,j,:))/dy !Jaap
                            dcbdy(i,j,:)=(ccb(i,j+1,:)-ccb(i,j,:))/dy
                        end do
                    end do
                    cvt(:,my+mbc,:) = cc(:,my+mbc,:)
                    cvbt(:,my+mbc,:) = ccb(:,my+mbc,:)
                else
                    cvt = cc
                    cvbt = ceqbg
                endif ! my>0

                ! Compute sedimnent transport in v-direction
                Svs = 0.0
                Svb = 0.0

                if (method=='Van_leer') then           ! Van leer
                    call Flux_vector(mbc,mx,my,u,v,h,cc,ccb,Sus,Sub,Svs,Svb)
                elseif (method=='normal') then       ! Van Thiel de Vries & Reniers 2008
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            ! suspended load
                            Svs(i,j,k)=(cvt(i,j,k)*v(i,j)*h(i,j)-Dc(i,j)*h(i,j)*dcsdy(i,j,k) &
                                    -facsl*cvt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdy(i,j))*wet(i,j)   !No bed slope term in suspended transport?
                            !bed load
                            Svb(i,j,k)=(cvbt(i,j,k)*v(i,j)*h(i,j)-facsl*cvbt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdy(i,j))*wet(i,j)
                        enddo
                    enddo
                end if


                do j=1-mbc,my+mbc-1
                    do i=1-mbc,mx+mbc
                        if(Svb(i,j,k)>0.0) then
                            pbbedv(i,j) = pbbed(i,j,1,k)
                        elseif(Svb(i,j,k)<0.0) then
                            pbbedv(i,j)= pbbed(i,j+1,1,k)
                        else
                            pbbedv(i,j)=0.50*(pbbed(i,j,1,k)+pbbed(i,j+1,1,k))
                        endif
                    enddo
                enddo

                pbbedv(:,my+mbc) = pbbedv(:,my+mbc-1)

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc
                        Svb(i,j,k) = pbbedv(i,j)*Svb(i,j,k)
                    enddo
                enddo

                ! solve concentration update

                if (my>0) then
                    do j=2-mbc,my+mbc-1
                        do i=2-mbc,mx+mbc-1
                            !suspended sediment transport
                            ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                            cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                    (h(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dy-Sus(i-1,j,k)*dy+&
                                    Svs(i,j,k)*dx-Svs(i,j-1,k)*dx)/(dx*dy)-ero1(i,j,k))) !please check

                            cc(i,j,k)=max(cc(i,j,k),0.00)
                            cc(i,j,k)=min(cc(i,j,k),cmax/2.0*h(i,j))
                            depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)

                            !bed sediment tranpsort
                            ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                            ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                    (h(i,j)*cc(i,j,k)/dt -((Sub(i,j,k)*dy-Sub(i-1,j,k)*dy+&
                                    Svb(i,j,k)*dx-Svb(i,j-1,k)*dx)/(dx*dy)-ero2(i,j,k))) !please check
                            ccb(i,j,k)=max(ccb(i,j,k),0.00)
                            ccb(i,j,k)=min(ccb(i,j,k),cmax/2.0*h(i,j))
                            depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                        enddo
                    enddo
                else
                    j=1
                    do i=2-mbc,mx+mbc-1
                        !suspended sediment transport
                        ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                        cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                    (h(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dy-Sus(i-1,j,k)*dy)/(dx*dy)-&
                                    ero1(i,j,k)))
                        cc(i,j,k)=max(cc(i,j,k),0.00)
                        cc(i,j,k)=min(cc(i,j,k),cmax/2.0*h(i,j))
                        depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)
                        !bed sediment tranpsort
                        ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                        ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                    (h(i,j)*ccb(i,j,k)/dt -((Sub(i,j,k)*dy-Sub(i-1,j,k)*dy)/(dx*dy)-&
                                    ero2(i,j,k)))
                        ccb(i,j,k)=max(ccb(i,j,k),0.00)
                        ccb(i,j,k)=min(ccb(i,j,k),cmax/2.0*h(i,j))
                        depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                    enddo
                endif

                do j= 1-mbc,my+mbc
                    do i= 1-mbc,mx+mbc
                        if (h(i,j) .lt. eps) then
                            cc(i,j,k) = 0.0
                            ccb(i,j,k) = 0.0
                        else
                            cc(i,j,k) = cc(i,j,k)/h(i,j)
                            ccb(i,j,k) = ccb(i,j,k)/h(i,j)
                        end if
                    end do
                end do
            end do
            ccbgt = ccb
            ccgt = cc
            Svsgt = Svs
            Susgt = Sus
            Svbgt = Svb
            Subgt = Sub
        end subroutine transus
    end module transport_module
