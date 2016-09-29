      program qsz
      implicit none
      real*8 r1,r2,cs,rho,q,s,z,x,gi,dr,up,upd,pwc
      real*8 ep(39:200)
      integer*4 i,loc1,loc2
      parameter(dr=0.04d0,pwc=0.001d0)

      open(3,file='qsz7.txt')
      open(4,file='pdm7.txt',action='read',status='old')
      open(5,file='ep7.txt',action='read',status='old')

      do i = 45,200
      read(5,*) r1,x
      ep(i)=x
      enddo

      do i = 1,601596
      read(4,*) r1,r2,cs,rho
      if(r1 .ge. 1.8d0 .and. r2 .ge. 1.8d0 .and. rho .ge. pwc) then
      q=0.5d0*(r1+r2)
      z=r1-r2
      s=sqrt(r1*r1+r2*r2-2.d0*cs*r1*r2)
!      if ( z .le. 1.d0 .and. s .le. 1.d0 ) then
      up=-log(gi(s)*rho)
      loc1=nint(r1/dr)
      loc2=nint(r2/dr)
      upd=up- 0.5d0*( ep(loc1) + ep(loc2) )
      write(3,*) q,s,z,upd
!      endif
      endif
      enddo


      end program


      real*8 function gi(x)
      implicit none
      real*8 x,tau,pi,lamb
      parameter(pi=acos(-1.d0))
      parameter(tau=0.025d0,lamb=12.1192393d0)
      
      gi=((4.d0*pi*lamb*tau)**(1.5d0))*exp((x*x)/(4.d0*lamb*tau))

      end function
     




