    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This Part is used to calculate critical velocity for each grain size classes!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine crtical_velocity3(mbc,mx,my,h,u,v,c,ub_cr,us_cr1,us_cr2)

            use sediment_module, only: rhos,rho,gmax,g,D,hcr,Trep,k0,m0,g,NL
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
            , us_cr1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) &
            , us_cr2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)


            call settling_velocity(mbc,mx,my,c,ws)

            !call bedrough(mbc,mx,my,u,v,h,c,z0)

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

                        ub_crt(i,j,k) = sum(spd_b_profile)/h(i,j)
                        us_cr1t(i,j,k) = sum(spd_s1_profile)/h(i,j)
                        us_cr2t(i,j,k) = sum(spd_s2_profile)/h(i,j)

                        ! correction 
                        ub_cr(i,j,k) = min(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr2(i,j,k) = max(ub_crt(i,j,k),us_cr1t(i,j,k),us_cr2t(i,j,k))
                        us_cr1(i,j,k) =ub_crt(i,j,k)+us_cr1t(i,j,k)+us_cr2t(i,j,k)-us_cr2(i,j,k)-ub_cr(i,j,k)

                    end do
                end do
            end do
        end subroutine critical_velocity3
    end module transport_module
