    module Flux

    !use params
        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used to put flux limiter
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine flux_limitor(mbc,mx,my,V_new,C_new_x,C_new_y,psix1,psix2,psiy1,psiy2)

            use Set_Precision, only: Prec
            use sediment_module, only: toler,limit_method,gmax,vareps
            implicit none

            ! Arguments
            integer, intent(in) :: mbc,mx,my

            real(kind=Prec), intent(in) :: V_new(1-mbc:mx+mbc,1-mbc:my+mbc,3),C_new_x(1-mbc:mx+mbc,1-mbc:my+mbc,gmax,2),&
                C_new_y(1-mbc:mx+mbc,1-mbc:my+mbc,gmax,2)

            !output
            real(kind=Prec), intent(inout) :: psix1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3),psix2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3), &
                psiy1(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3),psiy2(1-mbc:mx+mbc,1-mbc:my+mbc,gmax+3)

            !local

            integer :: i,j,k,ii
            Real(kind=Prec),Dimension(1-mbc:mx+mbc-1,1-mbc:my+mbc-1,gmax+3) :: Rx1,Rx2 !
            Real(kind=Prec),Dimension(1-mbc:mx+mbc-1,1-mbc:my+mbc-1,gmax+3) :: Ry1,Ry2  !
            Real(kind=Prec),Dimension(1-mbc:mx+mbc-1,1-mbc:my+mbc-1,gmax+3) :: denx, deny
            Real(kind=Prec),Dimension(1-mbc:mx+mbc,1-mbc:my+mbc-1,gmax) :: C_new_x_t,C_new_y_t


            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    do ii = 1, gmax
                        C_new_x_t(i,j,ii) = sum(C_new_x(i,j,ii,:))
                        C_new_y_t(i,j,ii) = sum(C_new_y(i,j,ii,:))
                    enddo
                enddo
            enddo

            do i = 1-mbc, mx+mbc-1
                do j = 1-mbc, my+mbc-1
                    do ii = 1, 3
                        denx(i,j,ii) = V_new(i+1,j,ii)-V_new(i,j,ii)
                        denx(i,j,ii) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))
                        deny(i,j,ii) = V_new(i,j+1,ii)-V_new(i,j,ii)
                        deny(i,j,ii) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))

                        Rx1(i,j,ii) = (V_new(i+2,j,ii)-V_new(i+1,j,ii))/denx(i,j,ii)
                        Rx2(i,j,ii) = (V_new(i,j,ii)-V_new(i-1,j,ii))/denx(i,j,ii)
                        Ry1(i,j,ii) = (V_new(i,j+2,ii)-V_new(i,j+1,ii))/deny(i,j,ii)
                        Ry2(i,j,ii) = (V_new(i,j,ii)-V_new(i,j-1,ii))/deny(i,j,ii)

                    enddo
                    do ii = 1, gmax
                        denx(i,j,ii+3) = C_new_x_t(i+1,j,ii)-C_new_x_t(i,j,ii)
                        denx(i,j,ii+3) = sign(max(abs(denx(i,j,ii)),toler),denx(i,j,ii))
                        deny(i,j,ii+3) = C_new_y_t(i,j+1,ii)-C_new_y_t(i,j,ii)
                        deny(i,j,ii+3) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))

                        Rx1(i,j,ii+3) = (C_new_x_t(i+2,j,ii)-C_new_x_t(i+1,j,ii))/denx(i,j,ii)
                        Rx2(i,j,ii+3) = (C_new_x_t(i,j,ii)-C_new_x_t(i-1,j,ii))/denx(i,j,ii)
                        Ry1(i,j,ii+3) = (C_new_y_t(i,j+2,ii)-C_new_y_t(i,j+1,ii))/deny(i,j,ii)
                        Ry2(i,j,ii+3) = (C_new_y_t(i,j,ii)-C_new_y_t(i,j-1,ii))/deny(i,j,ii)
                    enddo

                    if (limit_method == 'Vanleer') then
                        psix1(i,j,:) = (Rx1(i,j,:)+abs(Rx1(i,j,:)))/(1.0+Rx1(i,j,:))
                        psix2(i,j,:) = (Rx2(i,j,:)+abs(Rx2(i,j,:)))/(1.0+Rx2(i,j,:))
                        psiy1(i,j,:) = (Ry1(i,j,:)+abs(Ry1(i,j,:)))/(1.0+Ry1(i,j,:))
                        psiy2(i,j,:) = (Ry2(i,j,:)+abs(Ry2(i,j,:)))/(1.0+Ry2(i,j,:))
                    elseif (limit_method == 'VanAlbada') then
                        psix1(i,j,:) = (Rx1(i,j,:)+(Rx1(i,j,:))**2.0)/(1.0+(Rx1(i,j,:))**2.0)
                        psix2(i,j,:) = (Rx2(i,j,:)+(Rx2(i,j,:))**2.0)/(1.0+(Rx2(i,j,:))**2.0)
                        psiy1(i,j,:) = (Ry1(i,j,:)+(Ry1(i,j,:))**2.0)/(1.0+(Ry1(i,j,:))**2.0)
                        psiy2(i,j,:) = (Ry2(i,j,:)+(Ry2(i,j,:))**2.0)/(1.0+(Ry2(i,j,:))**2.0)
                    elseif (limit_method == 'minmod') then
                        do ii = 1, gmax+3
                            if (Rx1(i,j,ii)>0) then
                                psix1(i,j,ii) = min(1.0, Rx1(i,j,ii))
                            else
                                psix1(i,j,ii) = 0.0
                            endif
                            if (Rx2(i,j,ii)>0) then
                                psix2(i,j,ii) = min(1.0, Rx2(i,j,ii))
                            else
                                psix2(i,j,ii) = 0.0
                            endif
                            if (Ry1(i,j,ii)>0) then
                                psiy1(i,j,ii) = min(1.0, Ry1(i,j,ii))
                            else
                                psiy1(i,j,ii) = 0.0
                            endif
                            if (Ry2(i,j,ii)>0) then
                                psiy2(i,j,ii) = min(1.0, Ry2(i,j,ii))
                            else
                                psiy2(i,j,ii) = 0.0
                            endif
                        enddo
                    endif
                enddo
            enddo
            psix1(2-mbc,:,:)=psix1(3-mbc,:,:)
            psix1(1-mbc,:,:)=psix1(2-mbc,:,:)
            psix1(mx,:,:)=psix1(mx-1,:,:)
            psix1(mx+mbc-1,:,:)=psix1(mx,:,:)
            psix1(mx+mbc,:,:)=psix1(mx+mbc-1,:,:)

            psix2(2-mbc,:,:)=psix2(3-mbc,:,:)
            psix2(1-mbc,:,:)=psix2(2-mbc,:,:)
            psix2(mx,:,:)=psix2(mx-1,:,:)
            psix2(mx+mbc-1,:,:)=psix2(mx,:,:)
            psix2(mx+mbc,:,:)=psix1(mx+mbc-1,:,:)

            psix1(:,2-mbc,:)=psix1(:,3-mbc,:)
            psix1(:,1-mbc,:)=psix1(:,2-mbc,:)
            psix1(:,my,:)=psix1(:,my-1,:)
            psix1(:,my+mbc-1,:)=psix1(:,my,:)
            psix1(:,my+mbc,:)=psix1(:,my+mbc-1,:)

            psix2(:,2-mbc,:)=psix2(:,3-mbc,:)
            psix2(:,1-mbc,:)=psix2(:,2-mbc,:)
            psix2(:,my,:)=psix2(:,my-1,:)
            psix2(:,my+mbc-1,:)=psix2(:,my,:)
            psix2(:,my+mbc,:)=psix2(:,my+mbc-1,:)

            psiy1(2-mbc,:,:)=psiy1(3-mbc,:,:)
            psiy1(1-mbc,:,:)=psiy1(2-mbc,:,:)
            psiy1(mx,:,:)=psiy1(mx-1,:,:)
            psiy1(mx+mbc-1,:,:)=psiy1(mx,:,:)
            psiy1(mx+mbc,:,:)=psiy1(mx+mbc-1,:,:)

            psiy2(2-mbc,:,:)=psiy2(3-mbc,:,:)
            psiy2(1-mbc,:,:)=psiy2(2-mbc,:,:)
            psiy2(mx,:,:)=psiy2(mx-1,:,:)
            psiy2(mx+mbc-1,:,:)=psiy2(mx,:,:)
            psiy2(mx+mbc,:,:)=psiy2(mx+mbc-1,:,:)

            psiy1(:,2-mbc,:)=psiy1(:,3-mbc,:)
            psiy1(:,1-mbc,:)=psiy1(:,2-mbc,:)
            psiy1(:,my,:)=psiy1(:,my-1,:)
            psiy1(:,my+mbc-1,:)=psiy1(:,my,:)
            psiy1(:,my+mbc,:)=psiy1(:,my+mbc-1,:)

            psiy2(:,2-mbc,:)=psiy2(:,3-mbc,:)
            psiy2(:,1-mbc,:)=psiy2(:,2-mbc,:)
            psiy2(:,my,:)=psiy2(:,my-1,:)
            psiy2(:,my+mbc-1,:)=psiy2(:,my,:)
            psiy2(:,my+mbc,:)=psiy2(:,my+mbc,:)

            if (vareps ==0)  then
                psix1 = 0.0
                psix2 = 0.0
                psiy1 = 0.0
                psiy2 = 0.0
            end if
        end subroutine flux_limitor
    end module flux
