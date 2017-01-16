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
            integer :: i,j
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

    end module transport_module
