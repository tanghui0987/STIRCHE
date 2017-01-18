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
        !use update
        use Set_variable
        use sediment_module, only : morstart,morfac,struct,lmax,gmax,sourcesink,por, & !dzbdt add to sediment_module
            thick,toler,nd_var,trim,eps,totalthick
        !use Flux, only: Flux_vector


        IMPLICIT NONE

        contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used for containting the bed updating part                                    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine bed_update(mbc,mx,my,t,dt,dx,dy,u,v,h,aux,naux)

            !use set_variable

            IMPLICIT NONE

            ! argument
            integer, intent(in) :: mbc,mx,my,naux
            real(kind=Prec), intent(in) :: t,dt,dx,dy
            real(kind=Prec), intent(in) :: u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc)
            real(kind=Prec), intent(inout) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)
            real(kind=Prec), intent(inout) :: aux(naux,1-mbc:mx+mbc,1-mbc:my+mbc)

            !local
            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: dzbdt,wet,sederold
            integer,dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: indSus,indSub,indSvs,indSvb
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Sout,fac,dzg,edg
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: cu,cub,cv,cvb,ccg,ccbg,Susg, &
                    Subg,Svbg,Svsg,ero1,ero2,depo_ex1,depo_ex2
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ftem
            Real(kind=Prec) :: Savailable
            Real(kind=Prec),dimension(lmax,gmax) :: pb

            wet = 0.0

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    if (h(i,j)>eps) then
                        wet(i,j) = 1.0
                    endif
                enddo
            enddo
            zb = aux(1,:,:)

            call transus(mbc,mx,my,dx,dy,t,u,v,h,dt,cu,cub,cv,cvb,ccg,ccbg,Susg,Subg,Svbg,Svsg,ero1,ero2,depo_ex1,depo_ex2)


            dzbdt  = 0.0

            if (t>=morstart .and. morfac > .9990) then

                ! reduce sediment transports when hard layer comes to surface
                if (struct == 1) then
                    do k = 1,gmax
                        indSus = 0
                        indSub = 0
                        indSvs = 0
                        indSvb = 0
                        Sout   = 0.0
                        do j=2,my+mbc
                            do i=2,mx+mbc
                                ! fluxes at i,j
                                if (Subg(i,j,k) > 0.0) then      ! bed load x-direction
                                    indSub(i,j,k) = 1
                                    Sout(i,j,k) = Sout(i,j,k) + Subg(i,j,k)*dy
                                endif
                                ! fluxes at i-1,j
                                if (Subg(i-1,j,k) < 0.0 ) then   ! bed load x-direction
                                    Sout(i,j,k) = Sout(i,j,k) - Subg(i-1,j,k)*dy
                                endif
                                if (sourcesink==0) then
                                ! fluxes at i,j
                                    if (Susg(i,j,k) > 0.0 ) then     ! suspended load x-direction
                                        indSus(i,j,k) = 1
                                        Sout(i,j,k) = Sout(i,j,k) + Susg(i,j,k)*dy
                                    endif
                                    ! fluxes at i-1,j
                                    if (Susg(i-1,j,k) < 0.0 ) then   ! suspended load x-direction
                                        Sout(i,j,k) = Sout(i,j,k) - Susg(i-1,j,k)*dy
                                    endif
                                endif !sourcesink==0
                            enddo! imax
                        enddo! jmax
                        if (my>0) then
                            do j=2,my+mbc
                                do i=2,mx+mbc
                                    if (Svbg(i,j,k) > 0.0 ) then     ! bed load y-direction
                                        indSvb(i,j,k) = 1
                                        Sout(i,j,k) = Sout(i,j,k) + Svbg(i,j,k)*dx
                                    endif
                                    ! fluxes at i,j-1
                                    if (Svbg(i,j-1,k) < 0.0 ) then   ! bed load y-direction
                                        Sout(i,j,k) = Sout(i,j,k) - Svbg(i,j-1,k)*dx
                                    endif
                                    if (sourcesink==0) then
                                        if (Svsg(i,j,k) > 0.0 ) then     ! suspended load y-direction
                                            indSvs(i,j,k) = 1
                                            Sout(i,j,k) = Sout(i,j,k) + Svsg(i,j,k)*dx
                                        endif
                                    ! fluxes at i,j-1
                                        if (Svsg(i,j-1,k) < 0.0 ) then   ! suspended load y-direction
                                            Sout(i,j,k) = Sout(i,j,k) - Svsg(i,j-1,k)*dx
                                        endif
                                    endif ! sourcesink = 0
                                enddo !imax
                            enddo !jmax
                        endif !jmax>0

                        ! reduce the sediment flux/ transport due to hard structure layers
                        do j=2,my+mbc
                            do i=2,mx+mbc
                                ! reduction factor for cell outgoing sediment transports
                                Savailable = thick*pbbed(i,j,1,k)/morfac/dt*(1.0-por)*dx*dy
                                fac(i,j,k)  = min(1.0,Savailable/max(Sout(i,j,k),tiny(0.0)) )
                                ! fix sediment transports for the presence of a hard layer; remind indSus etc are 1 in cases of cell outgoing transports
                                ! updated S         all outgoing transports                  cell incoming transports
                                if (fac(i,j,k) < 1.0)then
                                    Subg(i,j,k)   = fac(i,j,k)*indSub(i,j,k)*Subg(i,j,k) &
                                                    + (1-indSub(i,j,k))*Subg(i,j,k)
                                    Subg(i-1,j,k) = fac(i-1,j,k)*(1-indSub(i-1,j,k))*Subg(i-1,j,k) &
                                                    + indSub(i-1,j,k)*Subg(i-1,j,k)
                                    if (my>0) then
                                        Svbg(i,j,k)   = fac(i,j,k)*indSvb(i,j,k)*Svbg(i,j,k) &
                                                        + (1-indSvb(i,j,k))*Svbg(i,j,k)
                                        Svbg(i,j-1,k) = fac(i,j-1,k)*(1-indSvb(i,j-1,k))*Svbg(i,j-1,k) &
                                                        + indSvb(i,j-1,k)*Svbg(i,j-1,k)
                                    endif
                                    if (sourcesink==0) then
                                        Susg(i,j,k)   = fac(i,j,k)*indSus(i,j,k)*Susg(i,j,k) &
                                                    + (1-indSus(i,j,k))*Susg(i,j,k)
                                        Susg(i-1,j,k) = fac(i-1,j,k)*(1-indSus(i-1,j,k))*Susg(i-1,j,k) &
                                                        + indSus(i-1,j,k)*Susg(i-1,j,k)
                                        if (my>0) then
                                            Svsg(i,j,k)   = fac(i,j,k)*indSvs(i,j,k)*Svsg(i,j,k) &
                                                        + (1-indSvs(i,j,k))*Svsg(i,j,k)
                                            Svsg(i,j-1,k) = fac(i,j,k)*(1-indSvs(i,j-1,k))*Svsg(i,j-1,k) &
                                                        + indSvs(i,j-1,k)*Svsg(i,j-1,k)
                                        endif !jmax > 0
                                    endif ! sourcesink = 0
                                endif !fac<1.0
                            enddo ! imax
                        enddo !jmax
                    enddo !gmax
                endif !struct == 1
                !if (trim=='soulsby_vanrijn') then           ! Soulsby van Rijn
                !   call sb_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)
                !elseif (trim=='vanthiel_vanrijn') then       ! Van Thiel de Vries & Reniers 2008
                !    call vt_vr(mbc,mx,my,u,v,h,Tsg,ceqbg,ceqsg)
                !elseif (trim=='flux') then
                !    call Flux_vector(mbc,mx,my,u,v,h,cc,ccb)
                !end if
                if (my>0) then
                    do j=2,my
                        do i=2,mx
                            ! bed level changes per fraction in this morphological time step in meters sand including pores
                            ! positive in case of erosion

                            if (sourcesink==0) then
                                do k = 1, gmax
                                    dzg(i,j,k)=wet(i,j)*morfac*dt/(1.0-por)*( &
                                    ! dz from sus transport gradients
                                            (Susg(i,j,k)*dy-Susg(i-1,j,k)*dy +&
                                            Svsg(i,j,k)*dx-Svsg(i,j-1,k)*dx +&
                                    ! dz from bed load transport gradients
                                            Subg(i,j,k)*dy-Subg(i-1,j,k)*dy+&
                                            Svbg(i,j,k)*dx-Svbg(i,j-1,k)*dx))/(dx*dy)! +&
                                end do
                                !source term
                            elseif (sourcesink==1) then
                                dzg(i,j,:)=morfac*dt/(1.0-por)* &
                                ! dz from sus transport gradients
                                    (ero1(i,j,:)-depo_ex1(i,j,:) +&
                                    (Subg(i,j,:)*dy-Subg(i-1,j,:)*dy +&
                                    Svbg(i,j,:)*dx-Svbg(i,j-1,:)*dx)/(dx*dy))
                                !source term
                                !(cu(i,j,:)+cv(i,j,:)-ceqsg(i,j,:))*h(i,j)/Tsg(i,j,:))
                            endif
                            if (gmax==1) then ! Simple bed update in case one fraction
                                zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                                dzbdt(i,j) = dzbdt(i,j)-sum(dzg(i,j,:))
                                sedero(i,j) = sedero(i,j)-sum(dzg(i,j,:))
                                totalthick(i,j) = max(0.0,totalthick(i,j)-sum(dzg(i,j,:)))
                            else
                                edg(i,j,:) = dzg(i,j,:)*(1.0-por)/dt
                                !call update_fractions(i,j,totalthick(i,j),pbbed(i,j,:,:),sum(dzg(i,j,:)),&
                                !       edg(i,j,:)/sum(edg(i,j,:)),dt,pb)
                                !pbbed(i,j,:,:)=pb
                                !edg(i,j,:) = dzg(i,j,:)*(1.0-por)/dt
                                !if (totalthick(i,j)>thick) then
                                !    totalnum(i,j) = nint(totalthick(i,j)/thick)
                                !    if (totalthick(i,j)>totalnum(i,j)*thick) then
                                !        totalnum(i,j)=totalnum(i,j)+1
                                !    endif
                                !else
                                !    totalnum(i,j) = 1
                                !endif
                                !if(dzbed(i,j,totalnum(i,j))<toler*thick) then
                                !    totalnum(i,j) = totalnum(i,j) - 1
                                !endif
                                !nd_var = totalnum(i,j)
                                !zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                                !sederold(i,j)=sedero(i,j)
                                !sedero(i,j) = sedero(i,j)+sum(dzg(i,j,:))
                                !if (abs(sedero(i,j)) > 1e-10) then
                                !    do k=1, gmax
                                !        pbbed(i,j,1,k) = abs((pbbed(i,j,1,k)*sederold(i,j)+dzg(i,j,k))/sedero(i,j))
                                !    end do
                                !end if
                                do k=1,gmax
                                    pbbed(i,j,1,k) = pbbed(i,j,1,k)/sum(pbbed(i,j,1,:))
                                end do
                                !call update_fractions(i,j,dzbed(i,j,:),pbbed(i,j,:,:),edg(i,j,:),sum(dzg(i,j,:)),dt)
                            endif
                        enddo ! imax+1
                        zb(2-mbc,j) = zb(1-mbc,j)
                    enddo ! jmax+1
                    zb(:,2-mbc) = zb(:,1-mbc)
                else
                    j=1
                    do i=1-mbc,mx+mbc
                        ! bed level changes per fraction in this morphological time step in meters sand including pores
                        ! positive in case of erosion
                        if (sourcesink==0) then
                            dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                            ! dz from sus transport gradients
                                (Susg(i,j,:)*dy-Susg(i-1,j,:)*dy +&
                            ! dz from bed load transport gradients
                                Subg(i,j,:)*dy-Subg(i-1,j,:)*dy)/(dx*dy))

                        elseif (sourcesink==1) then
                            dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                ! dz from sus transport gradients
                                ero1(i,j,:)-depo_ex1(i,j,:) + &
                                (Subg(i,j,:)*dy-Subg(i-1,j,:)*dy)/(dx*dy))
                        endif
                        if (gmax==1) then ! Simple bed update in case one fraction
                            zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                            dzbdt(i,j) = dzbdt(i,j)-sum(dzg(i,j,:))
                            sedero(i,j) = sedero(i,j)-sum(dzg(i,j,:))
                            totalthick(i,j) = max(0.0,totalthick(i,j)-sum(dzg(i,j,:)))
                        else ! multiple fractions...
                            edg(i,j,:) = dzg(i,j,:)*(1.0-por)/dt
                            !call update_fractions(i,j,totalthick(i,j),pbbed(i,j,:,:),sum(dzg(i,j,:)), &
                            !        edg(i,j,:)/sum(edg(i,j,:)),dt,pb)
                            !pbbed(i,j,:,:)=pb
                            !if (totalthick(i,j)>thick) then
                            !    totalnum(i,j) = nint(totalthick(i,j)/thick)
                            !    if (totalthick(i,j)>totalnum(i,j)*thick) then
                            !        totalnum(i,j)=totalnum(i,j)+1
                            !    endif
                            !else
                            !    totalnum(i,j) = 1
                            !endif
                            !if(dzbed(i,j,totalnum(i,j))<toler*thick) then
                            !    totalnum(i,j) = totalnum(i,j) - 1
                            !endif
                            !nd_var = totalnum(i,j)
                            !call update_fractions(i,j,dzbed(i,j,:),pbbed(i,j,:,:),edg(i,j,:),sum(dzg(i,j,:)),dt)
                        endif
                    enddo ! imax
                endif !jmax = 1
            endif !t
            !call avalanch(mbc,mx,my,dx,dy,dt,naux,aux,h,t)
        end subroutine bed_update

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This part is designed for avanlanching scheme
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine avalanch(mbc,mx,my,dx,dy,dt,naux,aux,h,t)

            use sediment_module, only : avalanching,morfac,lmax,gmax,hswitch,eps,wetslp,dryslp,&
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
                                dzb=sign(1.0d0,dzbdx(i,j))*(min(abs(dzbdx(i,j))-dzmax,thick))*dx !limit to the first layer
                                if (dzb >= 0.0) then
                                    ie = i+1                                        ! index erosion point
                                    id = i                                          ! index deposition point
                                    dzb=min(dzb,dzmax*dt/dx)                        ! make sure dzb is not in conflict with maximum erosion rate dzmax
                                    !dzb=min(dzb,totalthick(i+1,j))                 ! make sure dzb is not larger than sediment layer thickness
                                else
                                    ie = i                                          ! index erosion point
                                    id = i+1                                        ! index deposition point
                                    dzb=max(dzb,-dzmax*dt/dx)
                                    !dzb=max(dzb,-totalthick(i,j))
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
                                dzb=sign(1.0d0,dzbdy(i,j))*(min(abs(dzbdx(i,j))-dzmax,thick))*dy ! limit to the first layer
                                if (dzb >= 0.0) then
                                    je = j+1                                        ! index erosion point
                                    jd = j                                          ! index deposition point
                                    dzb=min(dzb,dzmax*dt/dy)
                                    !dzb=min(dzb,totalthick(i,j+1))
                                else
                                    je = j                                          ! index erosion point
                                    jd = j+1                                        ! index deposition point
                                    dzb=max(dzb,-dzmax*dt/dy)
                                    !dzb=max(dzb,-totalthick(i,j))
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
                                !call update_fractions(i,jd,totalthick(i,j),pbbed(i,j,:,:),dzleft,pbbed(i,j,1,:),dt,pb)
                                !pbbed(i,jd,:,:) = pb
                                !call update_fractions(i,je,totalthick(i,j),pbbed(i,j,:,:),-dzleft,pbbed(i,j,1,:),dt,pb)
                                !pbbed(i,je,:,:) = pb
                            end if !dzmax
                        end do !imax
                    end do !jmax
                    if (.not.aval) exit
                end do !morfac
            end if ! avalanching
        end subroutine avalanch

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used for remapping and recalculate grain-size distributin                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine update_fractions(is,js,dzs,pbs,dzbs,dpbs,dt,pbt)

            use sediment_module, only: lmax,gmax,por,nd_var,split,merge,toler,thick
            use Set_precision, only: Prec

            IMPLICIT NONE

            !Argument
            integer, intent(in) :: is, js !is: i; js: j; NL: numer of layers
            real(Kind= Prec), intent(in)                              :: dzbs,dzs,dt ! dab and total thickness
            !real(Kind= Prec), dimension(gmax),intent(in)          :: edgs !
            !real(Kind= Prec) , dimension(lmax),intent(in)    :: dzs !totalthickness
            real(Kind= Prec) , dimension(lmax,gmax),intent(in) :: pbs ! grain-size distribution
            real(Kind= Prec) , dimension(gmax),intent(in) :: dpbs ! change part of grain-size distribution

            !local
            integer                                         :: nl,nlf,nlb,nlo,i,j
            real(Kind= Prec)                                :: dzst,deltath_o,deltath_d,deltath_f
            !real(Kind= Prec) , dimension(lmax,gmax)         :: pbt

            real(Kind= Prec) , dimension(lmax,gmax),intent(inout) :: pbt

            dzst = dzs+dzbs

            ! figure out original number of layer and thickness of last layer
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

            ! figure out final number of layer and thickness of last layer
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

            ! figure out the number of layer and thickness of last layer for change part
            if (abs(dzbs)>thick) then
                nlb = nint(abs(dzbs)/thick)
                if (abs(dzbs)>nlb*thick) then
                    nlb=nlb+1
                endif
                if (dzbs<0) then
                    deltath_d = abs(dzbs) - nlb*thick
                else
                    deltath_d = dzbs - nlb*thick
                end if
            else
                nlb = 1
                if (dzbs<0) then
                    deltath_d = abs(dzbs)
                else
                    deltath_d = dzbs
                endif
            endif



            if (dzbs<0.0) then ! if erosion, from bottom to top need to change

                if (deltath_o > deltath_d) then
                    if (nlf == 1) then
                        pbt(nlf,:) = pbs(nlf,:)
                    else
                        pbt(nlf,:) = pbs(nlf,:)
                        do j =1, gmax
                            !pbt(1,j) = (pbs(1,j)*(thick-deltath_d)+pbs(2,j)*deltath_d)/thick ! problem?
                            pbt(1,j) = (pbs(i,j)*thick - dpbs(j)*deltath_d+pbs(2,j)*deltath_d)/thick
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
                            !pbt(1,j) = (pbs(1,j)*(thick-deltath_d)+pbs(2,j)*deltath_d)/thick ! problem?
                            pbt(1,j) = (pbs(i,j)*thick - dpbs(j)*deltath_d+pbs(2,j)*deltath_d)/thick
                        end do
                    else ! problem for first layer
                        do j =1, gmax
                            pbt(nlf,j) = (pbs(nlf-1,j)*(thick-deltath_d)+pbs(nlf,j)*deltath_d)/thick
                        end do
                        do i = 1, nlf-1
                            do j =1, gmax
                                pbt(i,j) = (pbs(i,j)*(thick-deltath_d)+pbs(i+1,j)*deltath_d)/thick
                            end do
                        end do
                    end if
                end if
            else ! if deposition, from top to bottom
                if (dzbs>thick) then
                    do i = 1, nlb-1
                        pbt(i,:) = dpbs
                    end do
                    if (deltath_d>thick-deltath_o) then
                        pbt(nlf,:) = pbs(nlo,:)
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(1,j) = (dpbs(j)*deltath_d+pbs(1,j)*(thick-deltath_d))/thick
                            end do
                        else
                            do i = nlb+1, nlf
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i-nlb,j)*deltath_d+pbs(i+1-nlb,j)*(thick-deltath_d))/thick
                                end do
                            end do
                        end if
                    else
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(nlf,j) = (dpbs(j)*deltath_d+pbs(nlo,j)*deltath_o)/(deltath_d+deltath_o)
                            end do
                        else
                            do j =1, gmax
                                pbt(nlf,j) = (pbs(nlo-1,j)*deltath_d+pbs(nlo,j)*deltath_o)/(deltath_d+deltath_o)
                            end do
                            do i = nlb+1, nlf
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i-nlb,j)*deltath_d+pbs(i+1-nlb,j)*(thick-deltath_d))/thick
                                end do
                            end do
                        end if
                    end if
                else
                    if (deltath_d>thick-deltath_o) then
                        pbt(nlf,:) = pbs(nlo,:)
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(1,j) = (dpbs(j)*deltath_d+pbs(1,j)*(thick-deltath_d))/thick
                            end do
                        else
                            do i = 1, nlf-1
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i,j)*deltath_d+pbs(i+1,j)*(thick-deltath_d))/thick
                                end do
                            end do
                        end if
                    else
                        if (nlo == 1) then
                            do j =1, gmax
                                pbt(1,j) = (dpbs(j)*deltath_d+pbs(1,j)*deltath_o)/(deltath_d+deltath_o)
                            end do
                        else
                            do j =1, gmax
                                pbt(nlf,j) = (pbs(nlo-1,j)*deltath_d+pbs(nlo,j)*deltath_o)/(deltath_d+deltath_o)
                            end do
                            do i = 1, nlf-1
                                do j =1, gmax
                                    pbt(i,j) = (pbs(i,j)*deltath_d+pbs(i+1,j)*(thick-deltath_d))/thick
                                end do
                            end do
                        end if
                    end if
                end if
            endif
            do i = 1, nlf
                do j = 1, gmax
                    pbt(i,j) = pbt(i,j)/sum(pbt(i,:))
                end do
            end do
        end subroutine update_fractions
    end module bedupdate_module
