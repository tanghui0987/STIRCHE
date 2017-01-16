    module transport_module

        use Set_Precision, only: Prec
        use Set_variable

        implicit none

        contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This Part is used to calculate Source term for finite volume method                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine transus(mbc,mx,my,dx,dy,time,u,v,h,dt,cut,cubt,cvt,cvbt,ccgt,ccbgt,Susgt,Subgt,Svbgt,Svsgt)

            use flux, only: Flux_vector
            use sediment_module, only: trim,gmax,morfac,por,D,thetanum,cmax,lmax,eps,facDC,nuh,nuhfac,rho,thick,method
            use Set_Precision, only: Prec
            use Set_variable,only: zb

            implicit none

            ! Arguments
            integer, intent(in) :: mbc,mx,my
            real(kind=Prec), intent(in) :: dx,dy,dt,time
            real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)
            !local

            integer :: i,j,k
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,hold,DR,Dc,wet,dzbdx, dzbdy, &
                pbbedu,pbbedv,dzbdy!hold?
            Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,fac,ero1,ero2,depo_ex1,depo_ex2, &
                cc,ccb,dcsdy,dcbdy,dcsdx,dcbdx
            Real(kind=Prec) :: exp_ero,facsl

            !out
            Real(kind=Prec),intent(inout) :: cut(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cubt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: cvt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),cvbt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: ccgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),ccbgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: Susgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Subgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)
            Real(kind=Prec),intent(inout) :: Svbgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax),Svsgt(1-mbc:mx+mbc,1-mbc:my+mbc,gmax)


            vmag2     = u**2+v**2
            cut = 0.0
            cubt = 0.0
            cvt = 0.0
            cvbt = 0.0
            dcsdx = 0.0
            dcbdx = 0.0
            dcsdy = 0.0
            dcbdy = 0.0
            DR = 0.0 !turn off roller model
            wet = 0.0
            exp_ero = 0.0
            facsl = 1.6

            ! calculate equibrium sediment concentration
            if (trim=='soulsby_vanrijn') then           ! Soulsby van Rijn
                call sb_vr(mbc,mx,my,u,v,h,c,Tsg,ceqbg,ceqsg)
            elseif (trim=='vanthiel_vanrijn') then       ! Van Thiel de Vries & Reniers 2008
                call vt_vr(mbc,mx,my,u,v,h,c,Tsg,ceqbg,ceqsg)
            end if

            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    if (h(i,j)>eps) then
                        wet(i,j) = 1.0
                    endif
                enddo
            enddo

            ! x-slope
            do i = 1-mbc, mx+mbc-1
                do j = 1-mbc, my+mbc
                    dzbdx(i,j) = (zb(i+1,j)-zb(i,j))/dx
                enddo
            enddo
            dzbdx(mx+mbc,:) = dzbdx(mx+mbc-1,:)
            ! y-slope
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc-1
                    dzbdy(i,j) = (zb(i,j+1)-zb(i,j))/dy
                enddo
            enddo
            dzbdy(:,my+mbc) = dzbdy(:,my+mbc-1)

            ! compute reduction factor for sediment sources due to presence of hard layers
            frc = pbbed(:,:,1,:)
            do k = 1,gmax
                do j= 1-mbc,my+mbc
                    do i= 1-mbc,mx+mbc
                        exp_ero = morfac*dt/(1.0-por)*h(i,j)*(ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k) &
                                + ceqbg(i,j,k)*frc(i,j,k)/dt) !ceqsg, ceqbg from vt_vr or sb_vr
                        fac(i,j,k) =min(1.0,thick*frc(i,j,k)/max(tiny(0.0),exp_ero))! limit to first layer
                    enddo
                enddo
            enddo

            ! compute diffusion coefficient, usually just based on horizontal viscosity
            Dc = facDc*(nuh+nuhfac*h*(DR/rho)**(1.0/3.0))
            cc = ccgt
            ccb = ccbgt
            !Sus = 0.0
            !Sub = 0.0
            do k = 1,gmax
                if (D(k)>0.002) then
                    print *, "WARNING: Grain size is larger than 2 mm"
                    ! set ceqsg to zero for gravel; in the future, we will add this module
                endif
                ! x-direction
                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc-1
                        if(u(i,j)>0.0) then
                            cut(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(min(i+1,mx+mbc-1),j,k)
                            cubt(i,j,k)=thetanum*frc(i,j,k)*ceqbg(i,j,k)+(1.0-thetanum)&
                                    *frc(min(i+1,mx+mbc-1),j,k)*ceqbg(min(i+1,mx+mbc-1),j,k)
                        elseif(u(i,j)<0.0) then
                            cut(i,j,k)=thetanum*cc(i+1,j,k)+(1.0-thetanum)*cc(max(i,2-mbc),j,k)
                            cubt(i,j,k)=thetanum*frc(i+1,j,k)*ceqbg(i+1,j,k)+(1.0-thetanum)&
                                    *frc(max(i,2-mbc),j,k)*ceqbg(max(i,2-mbc),j,k)
                            !cubt(i,j,k)=thetanum*ccb(i+1,j,k)+(1.0-thetanum)*ccb(max(i,2-mbc),j,k)
                        else
                            cut(i,j,k)=0.50*(cc(i,j,k)+cc(i+1,j,k))
                            cubt(i,j,k)=0.50*(frc(i,j,k)*ceqbg(i,j,k)+frc(i+1,j,k)*ceqbg(i+1,j,k))
                            !cubt(i,j,k)=0.50*(ccb(i,j,k)+ccb(i+1,j,k))
                        endif
                        dcsdx(i,j,:)=(cc(i+1,j,:)-cc(i,j,:))/dx
                        dcbdx(i,j,:)=(ccb(i+1,j,:)-ccb(i,j,:))/dx
                    enddo
                enddo
                cut(mx+mbc,:,:) = cc(mx+mbc,:,:)
                cubt(mx+mbc,:,:) = ccb(mx+mbc,:,:)

                Sus = 0.0
                Sub = 0.0

                if (method=='Van_leer') then           ! Van leer
                    call Flux_vector(mbc,mx,my,u,v,h,cc,ccb)
                elseif (method=='normal') then       ! Van Thiel de Vries & Reniers 2008
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                        ! suspended load
                            Sus(i,j,k)=(cut(i,j,k)*u(i,j)*h(i,j)-Dc(i,j)*h(i,j)*dcsdx(i,j,k) &
                            -facsl*cut(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdx(i,j,k))*wet(i,j)   !No bed slope term in suspended transport?
                        ! bed load
                            Sub(i,j,k)=(cubt(i,j,k)*u(i,j)*h(i,j)-facsl*cubt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdx(i,j,k))*wet(i,j)
                        enddo
                    enddo

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc-1
                        if(Sub(i,j,k)>0.0) then
                            pbbedu(i,j) = pbbed(i,j,1,k)
                        elseif(Sub(i,j,k)<0.0) then
                            pbbedu(i,j)= pbbed(i+1,j,1,k)
                        else
                            pbbedu(i,j)=0.50*(pbbed(i,j,1,k)+pbbed(i+1,j,1,k))
                        endif
                    enddo
                enddo

                pbbedu(mx+mbc,:) = pbbedu(mx+mbc-1,:)

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc
                        Sub(i,j,k) = pbbedu(i,j)*Sub(i,j,k)
                    enddo
                enddo

                !y-direction
                if (my>0) then
                    do j=1-mbc,my+mbc-1
                        do i=1-mbc,mx+mbc
                            if(v(i,j)>0) then
                                cvt(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(i,min(j+1,my+mbc),k)
                                cvbt(i,j,k)=thetanum*frc(i,j,k)*ceqbg(i,j,k)+(1.0-thetanum)&
                                        *frc(i,min(j+1,my+mbc),k)*ceqbg(i,min(j+1,my+mbc),k)
                                !cvbt(i,j,k)=thetanum*ccb(i,j,k)+(1.0-thetanum)*ccb(i,min(j+1,my+mbc),k)
                            elseif(v(i,j)<0) then
                                cvt(i,j,k)=thetanum*cc(i,j+1,k)+(1.0-thetanum)*cc(i,max(j,2-mbc),k)
                                cvbt(i,j,k)=thetanum*frc(i,j+1,k)*ceqbg(i,j+1,k)+(1.0-thetanum)&
                                            *frc(i,max(j,2),k)*ceqbg(i,max(j,2),k)
                                !cvbt(i,j,k)=thetanum*ccb(i,j+1,k)+(1.0-thetanum)*ccb(i,max(j,2-mbc),k)
                            else
                                cvt(i,j,k)=0.50*(cc(i,j,k)+cc(i,j+1,k)) !Jaap: cc instead of cv
                                cvbt(i,j,k)=0.50*(frc(i,j,k)*ceqbg(i,j,k)+frc(i,j+1,k)*ceqbg(i,j+1,k))
                                !cvbt(i,j,k)=0.50*(ccb(i,j,k)+ccb(i,j+1,k))
                            end if
                            dcsdy(i,j,:)=(cc(i,j+1,:)-cc(i,j,:))/dy !Jaap
                            dcbdy(i,j,:)=(ccb(i,j+1,:)-ccb(i,j,:))/dy
                        end do
                    end do
                    cvt(:,my+mbc,:) = cc(:,my+mbc,:)
                    cvbt(:,my+mbc,:) = ccb(:,my+mbc,:)
                else
                    cvt = cc
                    cvbt = ceqbg
                endif ! my>0

                ! Compute sedimnent transport in v-direction
                Svs = 0.0
                Svb = 0.0

                if (method=='Van_leer') then           ! Van leer
                    call Flux_vector(mbc,mx,my,u,v,h,cc,ccb)
                elseif (method=='normal') then       ! Van Thiel de Vries & Reniers 2008
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            ! suspended load
                            Svs(i,j,k)=(cvt(i,j,k)*v(i,j)*h(i,j)-Dc(i,j)*h(i,j)*dcsdy(i,j,k) &
                            -facsl*cvt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdy(i,j))*wet(i,j)   !No bed slope term in suspended transport?
                            !bed load
                            Svb(i,j,k)=(cvbt(i,j,k)*v(i,j)*h(i,j)-facsl*cvbt(i,j,k)*sqrt(vmag2(i,j))*h(i,j)*dzbdy(i,j))*wet(i,j)
                        enddo
                    enddo
                end if
                do j=1-mbc,my+mbc-1
                    do i=1-mbc,mx+mbc
                        if(Svb(i,j,k)>0.0) then
                            pbbedv(i,j) = pbbed(i,j,1,k)
                        elseif(Svb(i,j,k)<0.0) then
                            pbbedv(i,j)= pbbed(i,j+1,1,k)
                        else
                            pbbedv(i,j)=0.50*(pbbed(i,j,1,k)+pbbed(i,j+1,1,k))
                        endif
                    enddo
                enddo

                pbbedv(:,my+mbc) = pbbedv(:,my+mbc-1)

                do j=1-mbc,my+mbc
                    do i=1-mbc,mx+mbc
                        Svb(i,j,k) = pbbedv(i,j)*Svb(i,j,k)
                    enddo
                enddo

                ! solve concentration update

                if (my>0) then
                    do j=2-mbc,my+mbc-1
                        do i=2-mbc,mx+mbc-1
                            !suspended sediment transport
                            ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                            cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (h(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dy-Sus(i-1,j,k)*dy+&
                                        Svs(i,j,k)*dx-Svs(i,j-1,k)*dx)/(dx*dy)-ero1(i,j,k))) !please check

                            cc(i,j,k)=max(cc(i,j,k),0.00)
                            cc(i,j,k)=min(cc(i,j,k),cmax/2.0*h(i,j))
                            depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)

                            !bed sediment tranpsort
                            ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                            ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (h(i,j)*cc(i,j,k)/dt -((Sub(i,j,k)*dy-Sub(i-1,j,k)*dy+&
                                        Svb(i,j,k)*dx-Svb(i,j-1,k)*dx)/(dx*dy)-ero2(i,j,k))) !please check
                            ccb(i,j,k)=max(ccb(i,j,k),0.00)
                            ccb(i,j,k)=min(ccb(i,j,k),cmax/2.0*h(i,j))
                            depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                        enddo
                    enddo
                else
                    j=1
                    do i=2-mbc,mx+mbc-1
                        !suspended sediment transport
                        ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                        cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                (h(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dy-Sus(i-1,j,k)*dy)/(dx*dy)-&
                                ero1(i,j,k)))
                        cc(i,j,k)=max(cc(i,j,k),0.00)
                        cc(i,j,k)=min(cc(i,j,k),cmax/2.0*h(i,j))
                        depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)
                        !bed sediment tranpsort
                        ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*frc(i,j,k)/Tsg(i,j,k)
                        ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                    (h(i,j)*ccb(i,j,k)/dt -((Sub(i,j,k)*dy-Sub(i-1,j,k)*dy)/(dx*dy)-&
                                    ero2(i,j,k)))
                        ccb(i,j,k)=max(ccb(i,j,k),0.00)
                        ccb(i,j,k)=min(ccb(i,j,k),cmax/2.0*h(i,j))
                        depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                    enddo
                endif

                do j= 1-mbc,my+mbc
                    do i= 1-mbc,mx+mbc
                        if (h(i,j) .lt. eps) then
                            cc(i,j,k) = 0.0
                            ccb(i,j,k) = 0.0
                        else
                            cc(i,j,k) = cc(i,j,k)/h(i,j)
                            ccb(i,j,k) = ccb(i,j,k)/h(i,j)
                        end if
                    end do
                end do
            end do
            ccbgt = ccb
            ccgt = cc
            Svsgt = Svs
            Susgt = Sus
            Svbgt = Svb
            Subgt = Sub
        end subroutine transus
    end module transport_module
