! testing subroutine for CS2 Ames-1 PES (internal version E3b)
! --- prepared by xhuang@seti.org --- 09/22/2023
! gfortran test output: 
! V=     3.27037 cm-1
   
        program testcs2pes
        implicit double precision (a-h,o-z)
        implicit integer (i-n)

        dimension cart(3,3)

! initilize to read in pes data        
        call setupcs2longrange()

        b2a=0.529177249d0
        pi=dacos(-1.d0)

! cart(3,3) contains cartesian coordinates of C S1 S2 atoms, in bohr,
! C atom is at center
        cart=0.d0
        cart(2,3)= 1.55d0/b2a
        cart(3,3)=-1.55d0/b2a

       rcs1=dsqrt(sum((cart(1,:)-cart(2,:))**2))
       rcs2=dsqrt(sum((cart(1,:)-cart(3,:))**2))
       rss=sum((cart(2,:)-cart(3,:))**2)

       ang=(rcs1**2+rcs2**2-rss**2)/(2.d0*rcs1*rcs2)
       ang=dacos(min(1.d0,max(-1.d0,ang)))*180.d0/pi

! use angstrom and degrees as input r1, r2, and bond angle
! V is output in hartree       
       call cs2potlongrange(rcs1*b2a,rcs2*b2a,ang,V)       

        write(*,'(A,F12.5,A)')'V=',V*219474.63067d0,' cm-1'

       end

       
