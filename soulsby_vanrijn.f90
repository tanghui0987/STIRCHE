    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

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
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: wet,rhom,delta,urms,vmg,urms2,hloc,z0,Asb,Ass,term1
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

            call bedrough(mbc,mx,my,u,v,h,ub_cr,z0)

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
                            ceqb(i,j,k) = min(Asb(i,j)*term2(i,j)/hloc(i,j),cmax/gmax/2.0)
                            ceq(i,j,k) = ceqb(i,j,k)
                            ceqs(i,j,k) = 0.0
                        endif
                        if(term1(i,j)>us_cr2(i,j,k) .and. hloc(i,j)>eps) then
                            term2(i,j,k)=(term1(i,j)-us_cr1(i,j,k))**2.40
                            ceqs(i,j,k) = min(Ass(i,j)*term2(i,j)/hloc(i,j),cmax/gmax/2.0)
                            ceq(i,j,k) = ceqs(i,j,k)
                            ceqb(i,j,k) = 0.0
                        endif
                        if (term1(i,j)>ub_cr1(i,j,k) .and. term1(i,j)<us_cr2(i,j,k) .and. hloc(i,j)>eps) then
                            term2(i,j,k)=(term1(i,j)-us_cr1(i,j,k))**2.40
                            ceq(i,j,k) = (Ass(i,j)+Asb(i,j))*term2(i,j)/hloc(i,j)
                            perc = term1(i,j)/Us_cr2(i,j,k)
                            ceqb(i,j,k) = (1-perc)*min(ceq(i,j,k),cmax/gmax/2.0)
                            ceqs(i,j,k) = perc*min(ceq(i,j,k),cmax/gmax/2.0)
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
    end module transport_module
