      subroutine vibpot(rij,v,n,cart)
      implicit real*8 (a-h,o-z)
      dimension v(n),cart(n,3,3),cartin(3,3)
        dimension icoeopt(400),iopttmp(3)

       common/pcorcm/ ipcor
       common/ipcorcm/ icoeopt

      data ifirst/0/
      save
      if(ifirst.eq.0)then
       ifirst=1
       write(6,*)('cs2 --- tq5z+core')
       call setupcs2longrange()
       vzero=0d0
      end if

        b2a=0.529177249d0
	pi=dacos(-1.d0)

      do i=1,n
       cartin(1,:)=cart(i,:,1)
       cartin(2,:)=cart(i,:,2)
       cartin(3,:)=cart(i,:,3)
       rcs1=dsqrt(sum((cartin(1,:)-cartin(2,:))**2))
       rcs2=dsqrt(sum((cartin(1,:)-cartin(3,:))**2))
       rss=dsqrt(sum((cartin(2,:)-cartin(3,:))**2))

       ang=(rcs1**2+rcs2**2-rss**2)/(2.d0*rcs1*rcs2)
       ang=dacos(min(1.d0,max(-1.d0,ang)))*180.d0/pi
       call cs2potlongrange(rcs1*b2a,rcs2*b2a,ang,v(i))
       v(i)=v(i)+vzero

      end do
      return
      end

      subroutine cs2potlongrange(rcs1,rcs2,ang,v)
      implicit real*8 (a-h,o-z)

        dimension coe(1000),icoe(1000,3)
        dimension icoeopt(400),iopttmp(3)

        common/cs2pes/coe,emin
        common/ics2pes/icoe,Ncoe

        common/DeList/De1,De_1,De2,De_2,De3,De_3
        common/edamp/edp1,edp2,edp3,edp4,edp5,edp6
        common/pesref/rssref,rcsref,alpha1,alpha1b,alpha2,alpha2b

        common/pcorcm/ ipcor
        common/ipcorcm/ icoeopt

	pi=dacos(-1.d0)

! now we use scaling factor to calc sumstr
        rrefa=rcsref; rrefb=rcsref
        rcsd1=rcs1-rrefa; rcsd2=rcs2-rrefb
        rcsd1=rcsd1*rcsd1; rcsd2=rcsd2*rcsd2
        sumstr2=rcsd1+rcsd2
        sumstr4=rcsd1*rcsd1+rcsd2*rcsd2

        r1=rcs1-rcsref
        r2=rcs2-rcsref

	a3=1.d0+dcos(ang*dacos(-1.d0)/180.d0)

        v0=0.d0
        do i=1,Ncoe
          v0=v0+coe(i)*r1**icoe(i,1)*r2**icoe(i,2)*a3**(icoe(i,3))
	  if(icoe(i,1).ne.icoe(i,2))then
	    v0=v0+coe(i)*r2**icoe(i,1)*r1**icoe(i,2)*a3**(icoe(i,3))
	  end if
        end do

        ang2=dacos(-1.d0)-ang*dacos(-1.d0)/180.d0
        angref=150.d0*pi/180.d0-pi
        angx=min(-dabs(ang2)-angref,0.d0)
        bdamp2=angx**2;  bdamp4=angx**4
        etmp2=dexp(edp1*sumstr2+edp2*sumstr4+edp3*bdamp2+edp4*bdamp4)

	sumx0=V0
        V0=V0*etmp2

          enetmp1=De1*(1-dexp(-alpha1*(rcs1-rcsref)))**2 + De_1*(1-dexp(-alpha1*(rcs1-rcsref)))**4
          enetmp2=De2*(1-dexp(-alpha2*(rcs2-rcsref)))**2 + De_2*(1-dexp(-alpha2*(rcs2-rcsref)))**4

          acs2=ang; ax=dacos(-1.d0)*(180.d0-acs2)/360.d0
          enetmp3=De3*dsin(ax)**2 + De_3*dsin(ax)**4  !bending simulation

          edamp1=dexp(edp5*sumstr2+1.d0*edp6*sumstr4)
          enetmp3=edamp1*enetmp3

          etmp1=enetmp1+enetmp2+enetmp3 !+ enetmp2B

        V0=V0+etmp1
	
        V=V0/219474.63067d0-emin

      return
      end

      subroutine setupcs2longrange
      implicit real*8 (a-h,o-z)
        dimension coe(1000),icoe(1000,3)

        common/cs2pes/coe,emin
        common/ics2pes/icoe,Ncoe

        common/DeList/De1,De_1,De2,De_2,De3,De_3
        common/edamp/edp1,edp2,edp3,edp4,edp5,edp6
        common/pesref/rssref,rcsref,alpha1,alpha1b,alpha2,alpha2b

        common/ics2positive/icoe22,Ncoe22

        open(20,file='./cs2peslongrange.coeff.dat',status='old')
	read(20,*); read(20,*)emin
        read(20,*); read(20,*)rssref,rcsref,alpha1,alpha2
        read(20,*); read(20,*)De1,De_1,De2,De_2,De3,De_3
        read(20,*); read(20,*)edp1,edp2,edp3,edp4,edp5,edp6
        read(20,*)Ncoe
        do i=1,Ncoe
          read(20,*)icoe(i,1:3),coe(i)
        end do
        close(20)

        return
        end

