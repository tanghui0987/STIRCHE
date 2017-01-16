    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

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

            call crtical_velocity2(mbc,mx,my,h,u,v,ub_cr,us_cr1,us_cr2)

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
                            ceqb(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                            ceqs(i,j,k) = 0.0
                        endif
                        if(term1(i,j)>us_cr2(i,j,k) .and. h(i,j)>eps) then
                            term2(i,j)=0.0
                            term3(i,j)=(term1(i,j)-Us_cr1(i,j,k))**2.4
                            ceq(i,j,k) = (Asb(i,j)*term2(i,j)+Ass(i,j)*term3(i,j))/hloc(i,j)
                            ceqb(i,j,k) = 0.0
                            ceqs(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                        end if
                        if(term1(i,j)>us_cr1(i,j,k) .and. h(i,j)>eps .and. term1(i,j)< us_cr2(i,j,k)) then
                            term2(i,j)=(term1(i,j)-Ub_cr(i,j,k))**1.50
                            term3(i,j)=(term1(i,j)-Us_cr1(i,j,k))**2.40
                            ceqb(i,j,k) = Asb(i,j)*term2(i,j)/hloc(i,j)
                            ceqs(i,j,k) = Ass(i,j)*term3(i,j)/hloc(i,j)
                        end if
                        ! make sure only in water
                        ceqbg(i,j,k) = ceqb(i,j,k)*wet(i,j)
                        ceqsg(i,j,k) = ceqs(i,j,k)*wet(i,j)
                    enddo
                enddo
            enddo
        end subroutine vt_vr
    end module transport_module
