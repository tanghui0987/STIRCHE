    module fraction_module

        use Set_precision, only: Prec
        use set_variable
        !use sed

        IMPLICIT NONE

        contains


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used for remapping and recalculate grain-size distributin                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine update_fractions(is,js,dzs,pbs,dzbs,dpbs,dt,pbt)

            use sediment_module, only: lmax,gmax,por,nd_var,split,merge,toler,thick
            use Set_precision, only: Prec

            IMPLICIT NONE

            !Argument
            integer, intent(in) :: is, js, NLD !is: i; js: j; NL: numer of layers
            real(Kind= Prec), intent(in)                              :: dzbs,dzs ! dab and total thickness
            !real(Kind= Prec), dimension(gmax),intent(in)          :: edgs !
            real(Kind= Prec) , dimension(lmax),intent(in)    :: dzs !totalthickness
            real(Kind= Prec) , dimension(lmax,gmax),intent(in) :: pbs ! grain-size distribution
            real(Kind= Prec) , dimension(gmax),intent(in) :: dpbs ! change part of grain-size distribution

            !local
            integer                                         :: nl,nlf,nlb,i,j
            real(Kind= Prec)                                :: dzst,deltath_o,deltath_d,deltath_f
            real(Kind= Prec) , dimension(lmax,gmax)         :: pbt

            real(Kind= Prec) , dimension(lmax,gmax),intent(inout) :: pbt
            dzst = dzs+dzbs

            ! figure out number of layer and thickness of first layer
            if (dzs>thick) then
                nl = nint(dzs/thick)
                if (dzs>nl*thick) then
                    nl=nl+1
                endif
                deltath_o = dzs - nl*thick
            else
                nl = 1
                deltath_o = dzs
            endif

            if (dzst>thick) then
                nlf = nint(dzst/thick)
                if (dzst>nlf*thick) then
                    nlf=nlf+1
                endif
                deltath_f = dzst - nlf*thick
            else
                nlf = 1
                deltath_f = dzst
            endif

            if (abs(dzbs)>thick) then
                nlb = nint(abs(dzbs)/thick)
                if (abs(dzbs)>nlb*thick) then
                    nlb=nlb+1
                endif
                if (dzbs<0) then
                    deltath_d = abs(dzbs) - nlb*thick
                else
                    deltath_d = dzbs - nlb*thick
            else
                nlb = 1
                if (dzbs<0) then
                    deltath_d = abs(dzbs)
                else
                    deltath_d = dzbs
            endif


            if (dzbs<0.0) then ! if erosion, from bottom to top

                if (deltath_o > deltath_d)
                    if (nlf == 1) then
                        pbt(nlf,:) = pbs(nlf,:)
                    else
                        pbt(nlf,:) = pbs(nlf,:)
                        do j =1, gmax
                            pbt(1,j) = (pbs(1,j)*(thick-deltath_d)+pbs(2,j)*deltath_d)/thick
                        end do
                        if (nlf > 2) then
                            do i = 2, nlf-1
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i,j)*(thick-deltath_d)+pbs(i+1,j)*deltath_d)/thick
                                end do
                            end do
                        endif
                    endif
                else
                    if (nlf == 1) then
                        do j =1, gmax
                            pbt(1,j) = (pbs(1,j)*(thick-deltath_d)+pbs(2,j)*deltath_d)/thick
                        end do
                    else
                        do j =1, gmax
                            pbt(nlf,j) = (pbs(nlf-1,j)*(thick-deltath_d)+pbs(nlf,j)*deltath_d)/thick
                        end do
                        do i = 1, nlf-1
                            do j =1, gmax
                                pbt(i,j) = (pbs(i,j)*(thick-deltath_d)+pbs(i+1,j)*deltath_d)/thick
                            end do
                        end do
            else ! if deposition, from top to bottom
                if (dzbs>thick) then
                    do i = 1, nlb-1
                        pbt(i,:) = dpbs
                    end do
                    if (deltath_d>thick-deltath_o) then
                        pbt(nlf,:) = pbs(nl,:)
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(nlb,j) = (dpbs(j)*deltath_d+pbs(1,j)*(thick-deltath_d))/thick
                            end do
                        else
                            do i = nlb+1, nlf
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i-nlb,j)*deltath_d+pbs(i+1-nlb,j)*(thick-deltath_d)/thick
                                end do
                            end do
                        end if
                    else
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(nlf,j) = (dpbs(j)*deltath_d+pbs(nlo,j)*deltath_o)/thick
                            end do
                        else
                            do j =1, gmax
                                pbt(nlf,j) = (pbs(nlo-1,j)*deltath_d+pbs(nlo,j)*deltath_o)/thick
                            end do
                            do i = nlb+1, nlf
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i-nlb,j)*deltath_d+pbs(i+1-nlb,j)*(thick-deltath_d)/thick
                                end do
                            end do
                        end if
                    end if
                else
                    if (deltath_d>thick-deltath_o) then
                        pbt(nlf,:) = pbs(nl,:)
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(1,j) = (dpbs(j)*deltath_d+pbs(1,j)*(thick-deltath_d))/thick
                            end do
                        else
                            do i = 1, nlf
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i,j)*deltath_d+pbs(i+1,j)*(thick-deltath_d)/thick
                                end do
                            end do
                        end if
                    else
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(1,j) = (dpbs(j)*deltath_d+pbs(1,j)*deltath_o)/thick
                            end do
                        else
                            do j =1, gmax
                                pbt(nlf,j) = (pbs(nlo-1,j)*deltath_d+pbs(nlo,j)*deltath_o)/thick
                            end do
                        do i = 1, nlf
                            do j =1, gmax
                                pbt(i,j) = (pbs(i,j)*deltath_d+pbs(i+1,j)*(thick-deltath_d)/thick
                            end do
                        end do
                    end if
                end if
            end if
            do i = 1, nlf
                do j = 1, gmax
                    pbt(i,j) = pbt(i,j)/sum(pbt(i,:))
                end do
            end do
        end subroutine update_fractions
    end module fraction_module






