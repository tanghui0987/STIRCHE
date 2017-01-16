    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This Part is used to calculate critical velocity for each grain size classes!
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
            , us_cr1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
            , us_cr2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            call settling_velocity(mbc,mx,my,c,ws)

            call density(mbc,mx,my,c,rhom)

            delta  = (rhom-rho)/rho

            hloc = max(h,hcr)

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do k = 1, gmax

                        if (D(k) <= 0.0005) then
                            ub_crt(i,j,k) = 0.24*(delta*g)**(2.0/3.0)*(D(k)*Trep)**(1.0/3.0)
                        else
                            ub_crt(i,j,k) = 0.95*(delta*g)**(0.57)*D(k)**0.43*Trep**0.14
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


    end module transport_module
