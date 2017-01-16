!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (C) 2015 Virginia Tech, Sediment transport Processes Group    !
! Hui Tang and Robert Weiss                                               !
!                                                                         !
! tanghui@vt.edu                                                          !
!                                                                         !
! This library is free software; you can redistribute it and/or           !
! modify it under the terms of the GNU Lesser General Public              !
! License as published by the Free Software Foundation; either            !
! version 2.1 of the License, or (at your option) any later version.      !
!                                                                         !
! This library is distributed in the hope that it will be useful,         !
! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
! Lesser General Public License for more details.                         !
!                                                                         !
! You should have received a copy of the GNU Lesser General Public        !
! License along with this library; if not, write to the Free Software     !
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
! USA                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module bedupdate_module

        use Set_precision, only: Prec
        use transport_module, only: vt_vr,sb_vr,transus
        use update
        use Set_variable
        use sediment_module, only : morstart,morfac,struct,gmax,sourcesink,por, & !dzbdt add to sediment_module
            thick,toler,nd_var,trim,eps
        !use Flux, only: Flux_vector


        IMPLICIT NONE

        contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This part is designed for avanlanching scheme
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine avalanch(mbc,mx,my,dx,dy,dt,naux,aux,h,t)

            use sediment_module, only : avalanching,morfac,gmax,hswitch,eps,wetslp,dryslp,&
                por,nd_var,aval,trim,totalthick!dzbed,sedcal

            IMPLICIT NONE
            ! argument
            integer, intent(in) :: mbc,mx,my,naux
            real(kind=Prec), intent(in) :: t,dt,dx,dy
            real(kind=Prec), intent(inout) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)
            real(kind=Prec), intent(inout) :: aux(naux,1-mbc:mx+mbc,1-mbc:my+mbc)
            !local
            Integer         :: ie,id,jdz,je,jd,i,j,k,ndz,ii,indx
            Real(kind=Prec) :: sdz,dzb,one=1.00,dzmax,dzleft,dzavt,dzt
            Real(kind=Prec),dimension(lmax,gmax) :: pb
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: dzbdx,dzbdy,zb,dh,sedero,dzbdt,wet
            !
            zb = aux(1,:,:)

            wet = 0.0

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    if (h(i,j)>eps) then
                        wet(i,j) = 1.0
                    endif
                enddo
            enddo

            if (avalanching==1) then
                do ii=1,morfac
                    aval=.false.
                    dzbdx=0.0
                    dzbdy=0.0
                    ! x slope
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc-1
                            dzbdx(i,j)=(zb(i+1,j)-zb(i,j))/dx
                        enddo
                    enddo
                    ! y slope
                    do j=1-mbc,my+mbc-1
                        do i=1-mbc,mx+mbc
                            dzbdy(i,j)=(zb(i,j+1)-zb(i,j))/dy
                        enddo
                    enddo

                    indx = mx+1
                    do i=1-mbc,mx+mbc-1
                        do j=1-mbc,my+mbc
                            !decide Maximum bedlevel change due to avalanching
                            if(max(h(i,j),h(i+1,j))>hswitch+eps) then ! Jaap instead of hh
                                dzmax=wetslp
                                if (wet(i,j)==1) then ! tricks: seaward of indx (transition from sand to structure) wetslope is set to 0.03; TODO
                                    dzmax = max(wetslp,abs(dzbdx(i,j))*0.990)
                                endif
                            else
                                dzmax=dryslp
                            endif
                            if(abs(dzbdx(i,j))>dzmax .and. h(i+nint(max(0.0,sign(one,dzbdx(i,j)))),j)>eps) then
                                aval=.true.
                                ! the total mass need to be moved
                                dzb=min(sign(1.0,dzbdx(i,j))*(abs(dzbdx(i,j))-dzmax)*dx,thick) !limit to the first layer
                                if (dzb >= 0.0) then
                                    ie = i+1                                        ! index erosion point
                                    id = i                                          ! index deposition point
                                    dzb=min(dzb,dzmax*dt/dx)                        ! make sure dzb is not in conflict with maximum erosion rate dzmax
                                    dzb=min(dzb,totalthick(i+1,j))                 ! make sure dzb is not larger than sediment layer thickness
                                else
                                    ie = i                                          ! index erosion point
                                    id = i+1                                        ! index deposition point
                                    dzb=max(dzb,-dzmax*dt/dx)
                                    dzb=max(dzb,-totalthick(i,j))
                                endif
                                dzleft = abs(dzb)
                                zb(id,j) = zb(id,j)+dzleft
                                zb(ie,j) = zb(ie,j)-dzleft
                                dzbdt(id,j) = dzbdt(id,j)+dzleft
                                dzbdt(ie,j) = dzbdt(ie,j)-dzleft
                                sedero(id,j) = sedero(id,j)+dzleft
                                sedero(ie,j) = sedero(ie,j)-dzleft
                                totalthick(id,j) = max(0.0,totalthick(id,j)+dzleft)
                                totalthick(ie,j) = max(0.0,totalthick(ie,j)-dzleft)
                                h(id,j)  = h(id,j)+dzleft
                                h(ie,j)  = h(ie,j)-dzleft
                                call update_fractions(ie,j,totalthick(i,j),pbbed(i,j,:,:),-dzleft,pbbed(i,j,1,:),dt,pb)
                                pbbed(ie,j,:,:) = pb
                                call update_fractions(id,j,totalthick(i,j),pbbed(i,j,:,:),dzleft,pbbed(i,j,1,:),dt,pb)
                                pbbed(id,j,:,:) = pb
                            end if !dzmax
                        end do !jmax
                    end do !imax
                    !JJ: update y slopes after avalanching in X-direction seems more appropriat

                    do j=1-mbc,my+mbc-1 !
                        do i=1-mbc,mx+mbc
                            if(max(h(i,j),h(i,j+1))>hswitch+eps) then
                                dzmax=wetslp
                            else
                                dzmax=dryslp
                            end if

                            if(abs(dzbdy(i,j))>dzmax .and. h(i,j+nint(max(0.0,sign(one,dzbdy(i,j)))))>eps) then
                                aval=.true.
                                dzb=(sign(one,dzbdy(i,j))*(abs(dzbdy(i,j))-dzmax)*dy,thick) ! limit to the first layer
                                if (dzb >= 0.0) then
                                    je = j+1                                        ! index erosion point
                                    jd = j                                          ! index deposition point
                                    dzb=min(dzb,dzmax*dt/dy)
                                    dzb=min(dzb,totalthick(i,j+1))
                                else
                                    je = j                                          ! index erosion point
                                    jd = j+1                                        ! index deposition point
                                    dzb=max(dzb,-dzmax*dt/dy)
                                    dzb=max(dzb,-totalthick(i,j))
                                endif
                                dzleft = abs(dzb)
                                zb(i,jd) = zb(i,jd)+dzleft
                                zb(i,je) = zb(i,je)-dzleft
                                dzbdt(i,jd) = dzbdt(i,jd)+dzleft
                                dzbdt(i,je) = dzbdt(i,je)-dzleft
                                sedero(i,jd) = sedero(i,jd)+dzleft
                                sedero(i,je) = sedero(i,je)-dzleft
                                totalthick(i,jd) = max(0.0,totalthick(i,jd)+dzleft)
                                totalthick(i,je) = max(0.0,totalthick(i,je)-dzleft)
                                h(i,jd)  = h(i,jd)+dzleft
                                h(i,je)  = h(i,je)-dzleft
                                call update_fractions(i,jd,totalthick(i,j),pbbed(i,j,:,:),dzleft,pbbed(i,j,1,:),dt,pb)
                                pbbed(i,jd,:,:) = pb
                                call update_fractions(i,je,totalthick(i,j),pbbed(i,j,:,:),-dzleft,pbbed(i,j,1,:),dt,pb)
                                pbbed(i,je,:,:) = pb
                            end if !dzmax
                        end do !imax
                    end do !jmax
                if (.not.aval) exit
                end do !morfac
            end if ! avalanching
        end subroutine avalanch
    end module bedupdate_module
