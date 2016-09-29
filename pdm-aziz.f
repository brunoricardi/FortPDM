      program sqpw
      implicit none
      real*8 dr,r1,r2,r12,tau_max,lamb,pi,besvar,alpha1
      real*8 alpha2,y(1),asbes,arg_exp,dm,pwc,da,darg,dpn(1)
      real*8 ipn(1),s,qp,quad,ilj,epa,ljp,dmf,reff
      real*8 iaznum,azp,it
      real*8, allocatable :: grid(:),pw(:,:,:)
      real*8, allocatable :: pws(:,:,:),pw2(:,:,:)
      real*8, allocatable :: fc(:),gc(:),tau(:),twl(:)
      integer*4 i,ng,rn,j,kode,n,nz,awi,awj,nw,k,nsq
      integer*4 na,a,h,mu1,mu2,mode,isig,ierror,nad,nau
      integer*4 ncu,ncd,l,nas,nal
      parameter(ng=2000,rn=10,nw=40)
      parameter(dr=0.004d0,da=0.005d0)
      parameter(tau_max=1.d0/10240.d0)
      parameter(lamb=12.1192393d0)
      parameter(pi=acos(-1.d0))
      parameter(pwc=10.d0**(-5.d0))
      parameter(nsq=8,na=40)
!-----------------------------------------------------
!constants
      allocate(tau(0:nsq))
      allocate(fc(0:nsq))
      allocate(gc(0:nsq))
      allocate(twl(0:nsq))
      do i=0,nsq
      tau(i)=tau_max*(2.d0**(real(i)))
      fc(i)=(4.d0*pi*lamb*tau(i))**(-1.5d0)
      gc(i)=1.d0/(4.d0*lamb*tau(i))
      twl(i)=2.d0*sqrt(2.d0*lamb*tau(i))
      enddo
      quad=0.d0
!--------------------------------------------------------
!parameters to the bessel routine
      alpha1=0.5d0
      alpha2=alpha1-0.5d0
      kode=1
      n=1
!legendre routine
      mu1=0
      mu2=0
      mode=1
!-------------------------------------------------------------
!making up linear grid
      allocate(grid(0:ng))
      do i = 1,ng
      grid(i)=i*dr
      enddo
!angular control
      nad=int(1.d0/da) - na
      nau=int(1.d0/da) - 1
      allocate(pws(0:ng,0:ng,nad:nau+1))
      do i = 0,ng
      do j = 0,ng
      do a = nad,nau+1
      pws(i,j,a)=0.d0
      enddo
      enddo
      enddo
!----------------------------------------------------------
!allocating memory for the waves
      allocate(pw(0:ng,0:ng,0:nw))
      allocate(pw2(0:ng,0:ng,0:nw))
      do i=0,ng
      do j=0,ng
      do k=0,nw
      pw(i,j,k)=0.d0
      pw2(i,j,k)=0.d0
      enddo
      enddo
      enddo
!------------------------------------------------------------
!calculating the initial partial wave
      do k=0,nw
      alpha1=real(k)+0.5d0
      alpha2=alpha1-0.5d0
      do i = 1,ng
      do j = 1,ng
      r1=grid(i)
      r2=grid(j)
      r12=r1*r2
      dm=4.d0*pi*r12*fc(0)
      besvar=2.d0*gc(0)*r12
      if (besvar .lt. 50.d0) then
      call DBESI (besvar,alpha1,kode,n,y,nz)
      y(1)=y(1)*sqrt(pi/(2.d0*besvar))
      dm=dm*y(1)
      arg_exp=-(r1*r1+r2*r2)*gc(0)
      dm=dm*exp(arg_exp)
      else
      arg_exp=besvar-(alpha2*(alpha2+1.d0)/(2.d0*besvar))
     + -(r1*r1+r2*r2)*gc(0)
      dm=dm*exp(arg_exp)
      dm=dm*asbes(besvar,alpha2)
      endif
!interaction term
!      if (i==j) then
!      dm=dm*exp(-epa(r1,tau(0)))
!      else
!      dm=dm*exp(ilj(r1,r2,tau(0)))
!      endif

!      dm = dm*exp(-iaznum(r1,r2,tau(0)))

      pw(i,j,k)=dm
      enddo
      enddo
      enddo


      do i = 1, ng
      do j = 1, ng
      r1=grid(i)
      r2=grid(j)
      it = exp(-iaznum(r1,r2,tau(0)))
      do k = 0,nw
      pw(i,j,k) = pw(i,j,k)*it
      enddo
      enddo
      enddo





!-------------------------------------------------------------------------------
!squaring procedure
      do l=1,nsq      

      do h=0,nw
      do i = 1,ng
      do j = 1,ng
      do k = 1,ng
      qp=pw(i,k,h)*pw(k,j,h)*dr
      quad=quad+qp
      enddo
      pw2(i,j,h)=quad
      quad=0.d0
      enddo
      enddo


!correction near truncage point
      alpha1=real(h)+0.5d0
      alpha2=alpha1-0.5d0
      ncu=ng - int(twl(l)/dr)
      do i = ncu,ng
      do j = ncu,ng
      r1=grid(i)
      r2=grid(j)
      r12=r1*r2
      dm=4.d0*pi*r12*fc(l)
      besvar=2.d0*gc(l)*r12
      if (besvar .lt. 50.d0) then
      call DBESI (besvar,alpha1,kode,n,y,nz)
      y(1)=y(1)*sqrt(pi/(2.d0*besvar))
      dm=dm*y(1)
      arg_exp=-(r1*r1+r2*r2)*gc(l)
      dm=dm*exp(arg_exp)
      else
      arg_exp=besvar-(alpha2*(alpha2+1.d0)/(2.d0*besvar))
     + -(r1*r1+r2*r2)*gc(l)
      dm=dm*exp(arg_exp)
      dm=dm*asbes(besvar,alpha2)
      endif
!interaction term
      if (r1==r2) then
      dm=dm*exp(-tau(l)*azp(r1))
      else
      dm=dm*exp(-tau(l-1)*(azp(r1)+azp(r2)))
      endif
      pw2(i,j,h)=dm
      enddo
      enddo



!correction for small/large arguments
      nas=int(twl(l)/dr) + 1
      nal=ncu
      do i = 1,nas
      do j = nal,ng
      r1=grid(i)
      r2=grid(j)
      r12=r1*r2
      dm=4.d0*pi*r12*fc(l)
      besvar=2.d0*gc(l)*r12
      if (besvar .lt. 50.d0) then
      call DBESI (besvar,alpha1,kode,n,y,nz)
      y(1)=y(1)*sqrt(pi/(2.d0*besvar))
      dm=dm*y(1)
      arg_exp=-(r1*r1+r2*r2)*gc(l)
      dm=dm*exp(arg_exp)
      else
      arg_exp=besvar-(alpha2*(alpha2+1.d0)/(2.d0*besvar))
     + -(r1*r1+r2*r2)*gc(l)
      dm=dm*exp(arg_exp)
      dm=dm*asbes(besvar,alpha2)
      endif
      dmf=4.d0*pi*r1*r1*fc(l-1)
      besvar=2.d0*gc(l-1)*r1*r1     
      if (besvar .lt. 50.d0) then
      call DBESI (besvar,alpha1,kode,n,y,nz)
      y(1)=y(1)*sqrt(pi/(2.d0*besvar))
      dmf=dmf*y(1)
      arg_exp=-(2.d0*r1*r1)*gc(l-1)
      dmf=dmf*exp(arg_exp)
      else
      arg_exp=besvar-(alpha2*(alpha2+1.d0)/(2.d0*besvar))
     + -(2.d0*r1*r1)*gc(l-1)
      dmf=dmf*exp(arg_exp)
      dmf=dmf*asbes(besvar,alpha2)
      endif
!interaction term
      reff=pw(i,i,h)/dmf 
      dm=dm*reff*exp(-tau(l-1)*azp(r2))
      pw2(i,j,h)=dm
      enddo
      enddo
      

      do i = nal,ng
      do j = 1,nas
      r1=grid(i)
      r2=grid(j)
      r12=r1*r2
      dm=4.d0*pi*r12*fc(l)
      besvar=2.d0*gc(l)*r12
      if (besvar .lt. 50.d0) then
      call DBESI (besvar,alpha1,kode,n,y,nz)
      y(1)=y(1)*sqrt(pi/(2.d0*besvar))
      dm=dm*y(1)
      arg_exp=-(r1*r1+r2*r2)*gc(l)
      dm=dm*exp(arg_exp)
      else
      arg_exp=besvar-(alpha2*(alpha2+1.d0)/(2.d0*besvar))
     + -(r1*r1+r2*r2)*gc(l)
      dm=dm*exp(arg_exp)
      dm=dm*asbes(besvar,alpha2)
      endif
      dmf=4.d0*pi*r2*r2*fc(l-1)
      besvar=2.d0*gc(l-1)*r2*r2
      if (besvar .lt. 50.d0) then
      call DBESI (besvar,alpha1,kode,n,y,nz)
      y(1)=y(1)*sqrt(pi/(2.d0*besvar))
      dmf=dmf*y(1)
      arg_exp=-(2.d0*r2*r2)*gc(l-1)
      dmf=dmf*exp(arg_exp)
      else
      arg_exp=besvar-(alpha2*(alpha2+1.d0)/(2.d0*besvar))
     + -(2.d0*r2*r2)*gc(l-1)
      dmf=dmf*exp(arg_exp)
      dmf=dmf*asbes(besvar,alpha2)
      endif
!interaction term
      reff=pw(j,j,h)/dmf
      dm=dm*reff*exp(-tau(l-1)*azp(r1))
      pw2(i,j,h)=dm
      enddo 
      enddo
      

!resetting values after squaring
      do i =1,ng
      do j =1,ng
      pw(i,j,h)=pw2(i,j,h)
      pw2(i,j,h)=0.d0
      enddo
      enddo

!end of loop for the waves
      enddo
!end of loop for the squares
      enddo



!--------------------------------------------------------------------------
!summation of partial waves
      do a = nad,nau
      darg=a*da
      do i = 1,ng/rn
      do j = 1,ng/rn
      awi=i*rn
      awj=j*rn
      r12=grid(awi)*grid(awj)
      do h=0,nw
      dm=pw(awi,awj,h)*sqrt(2.d0*(2.d0*h + 1.d0))
      dm=dm/(4.d0*pi*r12)
      call DXNRMP (h,MU1,MU2,DARG,MODE,DPN,IPN,ISIG,IERROR)
      dm=dm*dpn(1)
      pws(awi,awj,a)=pws(awi,awj,a)+dm
      enddo
      enddo
      enddo
      enddo      

!cos(theta)=1
      do i = 1,ng/rn
      do j = 1,ng/rn
      awi=i*rn
      awj=j*rn
      r12=grid(awi)*grid(awj)
      do h = 0,nw
      dm = pw(awi,awj,h)*(2.d0*h + 1.d0)
      dm = dm/(4.d0*pi*r12)
      pws(awi,awj,nau+1)=pws(awi,awj,nau+1)+dm
      enddo
      enddo
      enddo


      open(2,file='pdm7.txt')
      do i=1,ng/rn
      do j=1,ng/rn
      do a=nad,nau
      awi=i*rn
      awj=j*rn
      r1=grid(awi)
      r2=grid(awj)
      r12=r1*r2
      darg=a*da
      if (pws(awi,awj,a) .ge. pwc) then
      write(2,*) r1,r2,darg,pws(awi,awj,a)
      endif
      enddo
      enddo
      enddo

      open(3,file='d7.txt')
      do i = 1,ng/rn
      do j = 1,ng/rn
      awi=i*rn
      awj=j*rn
      r1=grid(awi)
      r2=grid(awj)
      darg=1.d0
      if (pws(awi,awj,nau+1) .ge. pwc) then
      write(2,*) r1,r2,darg,pws(awi,awj,nau+1)
      write(3,*) r1,r2,darg,pws(awi,awj,nau+1)
      endif
      enddo
      enddo




!----------------------------------------------------------------------------




















!-------------------------------------------------------------

      end program




!-------------------------------------------------------------
!function to calculate asymptotic bessel funcs
       real*8 function asbes(x,y)
       implicit none
       real*8 x,y

       asbes=(1.d0/(2.d0*x))*(1.d0 - y*(y+1.d0)/((2.d0*x)**2.d0) +
     + y*(y+1.d0)*(y-2.d0)*(y+3.d0)/(3.d0*((2.d0*x)**3.d0)) +
     + y*(y+1.d0)*(5.d0*y*y+5.d0*y-12.d0)/(2.d0*((2.d0*x)**4.d0)))

       end function

!------------------------------------------------------------
!function to calculate interaction term
      real*8 function ilj(x,y,t)
      implicit none
      real*8 sig,eps,x,y,t,potail,rc,apt
      parameter(sig=2.556d0,eps=10.22d0)
      
      rc = 8.d0
      apt = (sig/rc)**6.d0
      potail=4.d0*eps*apt*(apt - 1.d0)


      ilj= (t/(x-y))*((4.d0*eps*(sig**12.d0)/11.d0)*(x**(-11.d0)
     + - y**(-11.d0)) + (4.d0*eps*(sig**6.d0)/5.d0)*(
     + y**(-5.d0) - x**(-5.d0))) + potail*t
       
      end function

!------------------------------------------------------------
!function for end-point approximation
      real*8 function epa(x,t)
      implicit none
      real*8 sig,eps,lamb,x,t,potail,rc,apt,a
      parameter(sig=2.556d0,eps=10.22d0)
      parameter(lamb=12.1192393d0)

      rc=8.d0
      apt = (sig/rc)**6.d0
      potail=4.d0*eps*apt*(apt - 1.d0)

      a = (sig/x)**6.d0
      epa=t*4.d0*eps*a*(a - 1.d0) - potail*t

      end function

!--------------------------------------------------------------------
!function for lennard-jones potential
      real*8 function ljp(x)
      implicit none
      real*8 sig,eps,x,potail,rc,apt,a
      parameter(sig=2.556d0,eps=10.22d0)
    
      rc=8.d0
      apt = (sig/rc)**6.d0
      potail=4.d0*eps*apt*(apt - 1.d0)
   
      a = (sig/x)**6.d0
      ljp=4.d0*eps*a*(a - 1.d0) - potail

      end function
     

!-----------------------------------------------------------------------------------------
!function to calculate aziz interaction numerically
      real*8 function iaznum(x,y,t)
      implicit none
      real*8 x,y,t
      real*8 rm,eps,alp,a,c6,c8,c10,d,rc
      real*8 xtail,utail,vtail
      real*8 dt,tp,rt,r,psum,u,v
      integer*4 nstep,i
      parameter(nstep = 100000)
      parameter(rm=2.9673d0, eps=10.8d0)
      parameter(a=0.54458046d0*10.d0**6.d0)
      parameter(alp=13.353384d0)
      parameter(c6=1.3732412d0, c8=0.4253758d0)
      parameter(c10=0.1781d0, d=1.241314d0)
      parameter(rc=8.d0)

      xtail = rc/rm
      utail = a*exp(-alp*xtail) - (c6/xtail**6.d0 + c8/xtail**8.d0
     + + c10/xtail**10.d0)
      vtail = eps*utail


      dt = t/real(nstep)
      psum = 0.d0


      do i = 0,nstep
      tp = i*dt

      rt = x + (y - x)*tp/t
      r = rt/rm

      if ( r .lt. d) then
      u = a*exp(-alp*r) - (c6/r**6.d0 + c8/r**8.d0
     +  + c10/r**10.d0)*exp(-(d/r - 1.d0)**2.d0)
      else
      u = a*exp(-alp*r) - (c6/r**6.d0 + c8/r**8.d0
     +  + c10/r**10.d0)
      endif
      v = eps*u - vtail

      psum = psum + v*dt

      enddo

      iaznum = psum


      end function


!------------------------------------------------------------------------
      real*8 function azp(x)
      real*8 x,vtail,utail,r,u,v
      real*8 rm,eps,alp,a,c6,c8,c10,d,rc

      parameter(rm=2.9673d0, eps=10.8d0)
      parameter(a=0.54458046d0*10.d0**6.d0)
      parameter(alp=13.353384d0)
      parameter(c6=1.3732412d0, c8=0.4253758d0)
      parameter(c10=0.1781d0, d=1.241314d0)
      parameter(rc=8.d0)

      xtail = rc/rm
      utail = a*exp(-alp*xtail) - (c6/xtail**6.d0 + c8/xtail**8.d0
     + + c10/xtail**10.d0)
      vtail = eps*utail

      r = x/rm
      if (r .lt. d) then
      u = a*exp(-alp*r) - (c6/r**6.d0 + c8/r**8.d0
     +  + c10/r**10.d0)*exp(-(d/r - 1.d0)**2.d0)
      else
      u = a*exp(-alp*r) - (c6/r**6.d0 + c8/r**8.d0
     +  + c10/r**10.d0)
      endif
      v = eps*u - vtail

      azp = v

      end function
















