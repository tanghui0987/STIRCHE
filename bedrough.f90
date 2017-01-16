    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This part is used to calculate bed roughness                                       !
        ! based on Wieberg 1989                                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine bedrough(mbc,mx,my,u,v,h,c,ub_cr,z0)

            use sediment_module, only: hcr,D,g,m0,rhos,rho,gmax
            use Set_Precision, only: Prec

            implicit none

            ! argument
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc),&
            v(1-mbc:mx+mbc,1-mbc:my+mbc), ub_cr(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)

            !local

            integer :: i,j
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,hloc,Dmm,zon,a2,ustarc,ustarcrit, &
                taub,taucrit,Tstar,delb,zos,rhom
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
                    ustarcrit(i,j) = sqrt(g*m0**2.0*ub_cr(i,j,int(gmax/2.0+1.0))**2.0/((rhom(i,j)-rho)*Dmm(i,j)&
                        *hloc(i,j)**(1.0/3.0)))
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
!********************************************************************
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
    end module transport_module


