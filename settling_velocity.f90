    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

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

