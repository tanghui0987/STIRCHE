    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

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

                        ! only sand
                        if (D(k) <= 0.0005) then
                            ub_crt(i,j,k) = 0.19*D(k)**0.1*dlog10(4*hloc(i,j)/D(k))
                        elseif (D(k) > 0.0005 .AND. D(k) <= 0.002) then
                            ub_crt(i,j,k) = 8.5*D(k)**0.6*dlog10(4*hloc(i,j)/D(k))
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

    end module transport_module
