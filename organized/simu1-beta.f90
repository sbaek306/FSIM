!simulation

!simus: number of simulations
!n: sample size
!p1: number of rows in Z
!pt1: number of free parameters in beta1, pt1=p1-1
!p2: number of columns in Z. = m1 + d (Bspline w/o intercept)
!d2: reduced number of columnss in Z
!pt2: number of free parameters in beta2, pt2=(p2-d2)*d2
!pt: total number of parameters, pt=pt1+pt2
!m1: number of knots
!d: degree of B-spline basis
!q: number of intervals in Composite Simpson's rule 
!num: number of different estimation procedures included

module commondata
  implicit none
  save
  integer,parameter :: simus=1000,n=500,p1=9,m1=3,d=3,q=100
  integer,parameter :: pt1=p1-1,p2=m1+d+1,pt2=p2-1,pt=pt1+pt2,m=m1+2*(d+1),iseed=1,lw=pt*(pt*3+13),num=3
  double precision,dimension(n,q+1,p1) :: x
  double precision,dimension(n,p1,p2) :: z
  double precision,dimension(n) :: y
  double precision,dimension(q+1,p2) :: bf
  double precision :: pi,hx,hxl,hxe,hy,b,hl,delta,mxb,scl
end module commondata

!start of program
program cs
  use commondata
  integer :: iters,eflag,i,j
  double precision,dimension(pt1) :: betat1
  double precision,dimension(pt2) :: betat2
  double precision,dimension(pt):: betat,zeros
  double precision,dimension(pt,num) :: betaest
  double precision,dimension(pt,pt,num) :: covest
  double precision,dimension(q+1):: alpha 
  double precision,dimension(q+1):: t 
  double precision,dimension(p2-1,p2-1):: bftbf, bftbfInv 
  double precision:: tmp
  integer :: sumo=0,suml=0,sume=0
  zeros=0
  pi=3.14159265358979
  delta=1e-3
  betat1=[1.0,1.2,1.5,0.5,-0.5,-1.5,-1.2,-1.0]  

! obtain spline coefficient estimates - betat2 - using LS
  call gendata1(p2,pt2,q,t,bf,alpha,betat2)
  betat(1:pt1)=betat1
  betat(pt1+1:pt)=betat2

! the bandwidth hx for kernel estimator regress
     hx=(1.0*n)**(-1./5.)*sqrt(sum(betat**2.+1.))*2.5
     hxl=(1.0*n)**(-1./5.)*sqrt(sum(betat**2.+1.))*0.5
     hxe=(1.0*n)**(-1./5.)*sqrt(sum(betat**2.+1.))*1.5
     hy=(1.0*n)**(-1./6.)*sqrt(sum(betat**2.+1.))*2.5
! the bandwidth b
     b=(1.0*n)**(-1./7.)*sqrt(sum(betat**2.+1.))*2.5
! initialize random number generator at iseed
  tmp=rand(iseed)
  iters=1

open(1,file='simu1-beta.txt') 
  write(1,*) betat
10 if (iters.le.simus) then
     print*,'iters=',iters
     call gendata2(n,p1,p2,pt1,pt2,betat1,betat2,q,t,bf,alpha,x,y,z) 

! the ideal estimator used by god
     call getcsgod(betat,betaest(:,1),covest(:,:,1),eflag)
     if (eflag.ne.1) then
        print*,'failed cs_god, eflag=',eflag
        sumo=sumo+1
        goto 10
     end if

! the local estimator assuming wrong model
     call getcslocal(betat,betaest(:,2),covest(:,:,2),eflag)
     if (eflag.ne.1) then
        print*,'failed cs_local, eflag=',eflag
        suml=suml+1
        goto 10
     end if

! our efficient estimator
     call getcseff(betat,betaest(:,3),covest(:,:,3),eflag)
     if (eflag.ne.1) then
        print*,'failed cs_eff, eflag=',eflag
        sume=sume+1
        goto 10
     end if

     do j=1,3
        write(1,*) betaest(:,j)
        write(1,*) (covest(i,i,j),i=1,pt)
     end do
     iters=iters+1
     goto 10

  end if
  print*,'oracle error=',sumo 
  print*,'lcoal error=',suml
  print*,'efficient error=',sume
  close(1)
  return 
end program cs


subroutine gendata1(p2,pt2,q,t,bf,alpha,betat2) 
  use commondata, only: m1,d,pi   
  integer,intent(in) :: p2,pt2,q 
  integer,parameter :: s=m1+d+1,m=m1+2*(d+1),ic=0 
  double precision,dimension(q+1,s):: basis  
  double precision,dimension(pt2),intent(out):: betat2
  double precision,dimension(q+1,p2),intent(out) :: bf
  double precision,dimension(q+1),intent(out):: alpha
  double precision,dimension(q+1),intent(out):: t
  double precision,dimension(p2-1,p2-1):: bftbf, bftbfInv 
  double precision,dimension(q+1):: one
  integer :: i,j,k
  double precision innerknots(m1), boundaryknots(2) 

  innerknots = [0.25, 0.5, 0.75]
  boundaryknots = [0.0, 1.0]

  do i=1,q+1
    t(i)=(dble(i)-1.0) / 100.0
  end do 

  alpha=sin(pi*t)+1.
  
! Obtain spline basis functions
  call splinebasis(d, q+1, m1, m, s, t, innerknots, boundaryknots, basis, ic)
  bf=basis
! Least square 
  bftbf=matmul(transpose(bf(1:q+1,2:p2)),bf(1:q+1,2:p2))
  call inv(p2-1,bftbf,bftbfInv)
  betat2=matmul(matmul(bftbfInv,transpose(bf(1:q+1,2:p2))),(alpha-bf(1:q+1,1)))

  return
end subroutine gendata1 


subroutine gendata2(n,p1,p2,pt1,pt2,beta1,beta2,q,t,bf,alpha,x,y,z) 
  use commondata, only: m1,d,pi   
  integer,intent(in) :: n,p1,pt1,p2,pt2,q 
  double precision,dimension(n,q+1,p1) :: x
  double precision,dimension(n),intent(out) :: y
  double precision,dimension(n,p1,p2),intent(out) :: z
  double precision,dimension(pt1),intent(in):: beta1
  double precision,dimension(pt2),intent(in):: beta2
  double precision,dimension(q+1),intent(in) :: t
  double precision,dimension(q+1,p2),intent(in) :: bf
  double precision,dimension(q+1),intent(in):: alpha
  double precision,dimension(p1):: betal1
  double precision,dimension(p2):: betal2
  double precision :: a,b,a1,a2,b1,b2,a3,a4,a5,b3,b4,b5,dinvnr,integral,zb
  double precision :: a6,a7,a8,a9,b6,b7,b8,b9
  double precision :: psi1,psi2,psi3,psi4
  integer :: i,j,k,h,kflag
  
  betal1(1:pt1)=beta1
  betal1(p1)=1.0
  betal2(1)=1.0
  betal2(2:p2)=beta2 

  do i=1,n
     do h=1,q+1
        a1=rand(0)*10.-5.
        a2=rand(0)*10.-5.
        a3=rand(0)*10.-5.
        a4=rand(0)*10.-5.
        a5=rand(0)*10.-5.
        a6=rand(0)*10.-5.
        a7=rand(0)*10.-5.
        a8=rand(0)*10.-5.
        a9=rand(0)*10.-5.
        b1=1.-a1
        b2=1.-a2
        b3=1.-a3
        b4=1.-a4 
        b5=1.-a5
        b6=1.-a6
        b7=1.-a7
        b8=1.-a8
        b9=1.-a9 
        x(i,h,1)=a1
        x(i,h,2)=a2
        x(i,h,3)=a3
        x(i,h,4)=a4
        x(i,h,5)=a5
        x(i,h,6)=a6
        x(i,h,7)=a7
        x(i,h,8)=a8
        x(i,h,9)=a9
    end do 

    do j=1,p1
      do k=1,p2  
        call compositesimpson(x(i,:,j),bf(:,k),integral)
        z(i,j,k)=integral 
      end do
    end do
    call compositesimpson(matmul(x(i,:,:),betal1),alpha,integral)
    a=rand(0)
    b=1.0-a
    y(i) = integral + dinvnr(a,b)
  end do 
  return
end subroutine gendata2

subroutine getcsgod(betat,betaest,covest,eflag)
  use commondata, only: pt,lw
  double precision,dimension(pt),intent(in) :: betat
  double precision,dimension(pt),intent(out) :: betaest
  double precision,dimension(pt,pt),intent(out) :: covest
  integer,intent(out) :: eflag
  double precision,dimension(pt) :: fval
  double precision :: tol
  double precision,dimension(lw) :: wv
  external estcsgod
  tol=1e-8
  betaest=betat 
  call hybrd10(estcsgod,pt,betaest,fval,tol,eflag,wv,lw)
  if (eflag.ne.1) then
     print*,'failed in hybrid god, eflag=',eflag
     return
  end if
  call covcsgod(betaest,covest) 
  return
end subroutine getcsgod


subroutine estcsgod(pt,betaest,fval,eflag)
  use commondata, only: n
  integer,intent(in) :: pt
  double precision,dimension(pt),intent(in) :: betaest
  integer,intent(in) :: eflag
  double precision,dimension(pt),intent(out) :: fval
  double precision,dimension(n,pt) :: yout
  integer :: i
  call csgod(betaest,yout)
  fval=sum(yout,1)/n
  return 
end subroutine estcsgod

subroutine csgod(betaest,yout)
  use commondata, only: pt,n,z,y,hx,p1,p2,pt1,pt2 !,tzb
  double precision,dimension(pt),intent(in) :: betaest
  double precision,dimension(n,pt),intent(out) :: yout
  double precision,dimension(n,pt) :: zv,zv0
  double precision,dimension(p1) :: betal1
  double precision,dimension(p2) :: betal2
  double precision,dimension(n) :: zb
  double precision,dimension(pt) :: mz
  double precision,dimension(p1,p2) :: my
  double precision,dimension(pt1) :: z1
  double precision,dimension(pt2) :: z2
  double precision,dimension(n) :: fb,mfb
  double precision :: z0,tmp,t1  
  integer :: i,j,k,one
  one=1
  betal1(1:pt1)=betaest(1:pt1)
  betal1(p1)=1.
  betal2(1)=1.
  betal2(2:p2)=betaest(pt1+1:pt)

  do i=1,n
     z1=matmul(z(i,1:pt1,:),betal2)
     z2=matmul(transpose(z(i,:,2:p2)),betal1)
     zv(i,1:pt1)=z1
     zv(i,pt1+1:pt)=z2
     zb(i)=dot_product(matmul(betal1,z(i,:,:)),betal2)
  end do
 
  do i=1,n

    z0=zb(i)

! kernel estimation,
    call regress(n,pt,z0,hx,zb,zv,mz) 
    zv0(i,:)=zv(i,:)-mz
  end do

  call truedlogeta(n,zb,y,fb)
  do i=1,n
     z0=zb(i)
     call regress(n,one,z0,hx,zb,fb,tmp) 
     mfb(i)=fb(i)-tmp
  end do
  do i=1,n
     yout(i,:)=zv0(i,:)*mfb(i)
  end do
  return
end subroutine csgod

subroutine getcslocal(betat,betaest,covest,eflag)
  use commondata, only: pt,lw
  double precision,dimension(pt),intent(in) :: betat
  double precision,dimension(pt),intent(out) :: betaest
  double precision,dimension(pt,pt),intent(out) :: covest
  integer,intent(out) :: eflag
  double precision,dimension(pt) :: fval
  double precision :: tol
  double precision,dimension(lw) :: wv
  external estcslocal
  tol=1e-5
  betaest=betat
  call hybrd10(estcslocal,pt,betaest,fval,tol,eflag,wv,lw)
  if (eflag.ne.1) then
     print*,'failed in hybrid local, eflag=',eflag
     return
  end if
  call covcslocal(betaest,covest)
  return
end subroutine getcslocal

! the local estimating equation
subroutine estcslocal(pt,betaest,fval,eflag)
  use commondata, only: n
  integer,intent(in) :: pt
  double precision,dimension(pt),intent(in) :: betaest
  integer,intent(in) :: eflag
  double precision,dimension(pt),intent(out) :: fval
  double precision,dimension(n,pt) :: yout
  integer :: i
  call cslocal(betaest,yout)
  fval=sum(yout,1)/n
  return
end subroutine estcslocal

! the local estimating function components
subroutine cslocal(betaest,yout)
  use commondata, only: pt,n,z,y,hx,hxl,p1,p2,pt1,pt2
  double precision,dimension(pt),intent(in) :: betaest
  double precision,dimension(n,pt),intent(out) :: yout
  double precision,dimension(n,pt) :: zv,zv0
  double precision,dimension(p1) :: betal1
  double precision,dimension(p2) :: betal2
  double precision,dimension(n) :: zb
  double precision,dimension(pt) :: mz
  double precision,dimension(pt1) :: z1
  double precision,dimension(pt2) :: z2
  double precision,dimension(n) :: fb,mfb
  double precision :: z0,tmp 
  integer :: i,one
  one=1
  betal1(1:pt1)=betaest(1:pt1)
  betal1(pt1+1:p1)=1
  betal2(1)=1
  betal2(2:p2)=betaest(pt1+1:pt)
  
  do i=1,n
     z1=matmul(z(i,1:pt1,2:p2),betal2(2:p2))+z(i,1:pt1,1)
     z2=matmul(betal1(1:pt1),z(i,1:pt1,2:p2))+z(i,2,2:p2) 
     zv(i,1:pt1)=z1
     zv(i,pt1+1:pt)=z2
     zb(i)=dot_product(matmul(betal1,z(i,:,:)),betal2)
  end do
  
  do i=1,n
     z0=zb(i)
! kernel estimation,
     call regress(n,pt,z0,hxl,zb,zv,mz) 
     zv0(i,:)=zv(i,:)-mz
  end do
  call localdlogeta(n,zb,y,fb)
  do i=1,n
     z0=zb(i)
     call regress(n,one,z0,hxl,zb,fb,tmp) 
     mfb(i)=fb(i)-tmp
  end do
  do i=1,n
     yout(i,:)=zv0(i,:)*mfb(i)
  end do
  return
end subroutine cslocal

subroutine getcseff(betat,betaest,covest,eflag)
  use commondata, only: pt,lw
  double precision,dimension(pt),intent(in) :: betat
  double precision,dimension(pt),intent(out) :: betaest
  double precision,dimension(pt,pt),intent(out) :: covest
  integer,intent(out) :: eflag
  double precision,dimension(pt) :: fval
  double precision :: tol
  double precision,dimension(lw) :: wv
  external estcseff
  tol=1e-5
  betaest=betat
  call hybrd10(estcseff,pt,betaest,fval,tol,eflag,wv,lw)
  if (eflag.ne.1) then
     print*,'failed in hybrid eff, eflag=',eflag
     return
  end if
  call covcseff(betaest,covest)
  return
end subroutine getcseff

! the efficient estimating equation
subroutine estcseff(pt,betaest,fval,eflag)
  use commondata, only: n
  integer,intent(in) :: pt
  double precision,dimension(pt),intent(in) :: betaest
  integer,intent(in) :: eflag
  double precision,dimension(pt),intent(out) :: fval
  double precision,dimension(n,pt) :: yout
  integer :: i
  call cseff(betaest,yout)
  fval=sum(yout,1)/n
  return
end subroutine estcseff

! the efficient estimating function components
subroutine cseff(betaest,yout)
  use commondata, only: pt,n,z,y,hx,hxe,p1,p2,pt1,pt2 
  double precision,dimension(pt),intent(in) :: betaest
  double precision,dimension(n,pt),intent(out) :: yout
  double precision,dimension(n,pt) :: zv,zv0
  double precision,dimension(p1) :: betal1
  double precision,dimension(p2) :: betal2
  double precision,dimension(n) :: zb
  double precision,dimension(pt) :: mz
  double precision,dimension(p1,p2) :: my
  double precision,dimension(pt1) :: z1
  double precision,dimension(pt2) :: z2
  double precision,dimension(n) :: fb,mfb
  double precision :: z0,tmp,t1  
  integer :: i,j,k,one
 one=1
  betal1(1:pt1)=betaest(1:pt1)
  betal1(p1)=1.
  betal2(1)=1.
  betal2(2:p2)=betaest(pt1+1:pt)

  do i=1,n
     z1=matmul(z(i,1:pt1,:),betal2)
     z2=matmul(transpose(z(i,:,2:p2)),betal1)
     zv(i,1:pt1)=z1
     zv(i,pt1+1:pt)=z2
     zb(i)=dot_product(matmul(betal1,z(i,:,:)),betal2)
  end do

   do i=1,n
     z0=zb(i)
! kernel estimation,
     call regress(n,pt,z0,hxe,zb,zv,mz) 
     zv0(i,:)=zv(i,:)-mz
   end do

  call dlogeta(n,zb,y,fb)
  do i=1,n
     z0=zb(i)
    call locallinear0(n,one,z0,hxe,zb,fb,tmp) 
     mfb(i)=fb(i)-tmp
  end do
  do i=1,n
     yout(i,:)=zv0(i,:)*mfb(i)
  end do
  return
end subroutine cseff

subroutine truelogeta(n,t,y,f)  
  use commondata, only: pi
  integer,intent(in) :: n
  double precision,dimension(n),intent(in) :: t
  double precision,dimension(n),intent(in) :: y
  double precision,dimension(n),intent(out) :: f
  double precision,dimension(n):: m,v
  m=1.0*t 
  v=1.    
  f=-(y-m)**2./2./v-log(2.*pi)/2.-log(v)/2.
  return
end subroutine truelogeta

!true dlog(\eta_2(y,t))/dt function
subroutine truedlogeta(n,t,y,f)
  integer,intent(in) :: n
  double precision,dimension(n),intent(in) :: t
  double precision,dimension(n),intent(in) :: y
  double precision,dimension(n),intent(out) :: f
  double precision,dimension(n):: m,v,t1,t2,dm,dv
  m=1.0*t ! mean 
  v=1.0   ! variance
  dm=1.0
  dv=0.0
  t1=(y-m)/v
  t2=t1**2./2.-0.5/v
  f=t1*dm+t2*dv
  return
end subroutine truedlogeta

subroutine localdlogeta(n,t,y,f)
  integer,intent(in) :: n
  double precision,dimension(n),intent(in) :: t
  double precision,dimension(n),intent(in) :: y
  double precision,dimension(n),intent(out) :: f
  double precision,dimension(n):: m,v,t1,t2,dm,dv
  m=1.0*t 
  v=1./log(2.+t**2.)
  dm=1. 
  dv=-1./(log(2.+t**2.))**2*2.*t/(2+t**2.)
  t1=(y-m)/v
  t2=t1**2./2.-0.5/v
  f=t1*dm+t2*dv
  return
end subroutine localdlogeta

! the nonparametric estimation of dlog(\eta_2(y,t))/dt function
subroutine dlogeta(n,t,y,f)
  use commondata, only: hx,hy,b
  integer,intent(in) :: n
  double precision,dimension(n),intent(in) :: t
  double precision,dimension(n),intent(in) :: y
  double precision,dimension(n),intent(out) :: f
  double precision,dimension(n) :: yy,ys,f1,f2,ft1
  double precision :: eps,z0
  integer :: i,one
  one=1
  eps=1e-8
  do i=1,n
     ys=y-y(i)
     call kh(n,ys,b,yy)
     z0=t(i)
     call locallinear(n,one,z0,hx,hy,t,yy,f1(i),f2(i))
  end do
  f=f2/(f1+eps)
  return
end subroutine dlogeta

! regress y on x, univariate x
subroutine regress(n,p,x0,h,x,y,my)
  integer,intent(in) :: n,p
  double precision,intent(in) :: x0,h
  double precision,dimension(n),intent(in) :: x
  double precision,dimension(n,p),intent(in) :: y
  double precision,dimension(p),intent(out) :: my
  double precision,dimension(n) :: t,s
  double precision,dimension(p) :: t2
  double precision :: t1
  integer :: i
  t=x-x0
  call kh(n,t,h,s)
  t1=sum(s)
  t2=0
  do i=1,n
     t2=t2+y(i,:)*s(i)
  end do
  my=t2/t1
  return
end subroutine regress

!local linear estimation of p nonparametric functions and their derivatives with respect to covariate, evaluated at x0, 
subroutine locallinear0(n,p,z0,hx,z,y,f)
  integer,intent(in) :: n,p
  double precision,intent(in) :: z0,hx
  double precision,dimension(n),intent(in) :: z
  double precision,dimension(n,p),intent(in) :: y
  double precision,dimension(p),intent(out) :: f
  double precision,dimension(n) :: u,w,tmp1,tmp2
  double precision :: a,b,d,deno,t1,t2,num1,num2,eps
  integer :: i,j
  eps=1e-8
  u=z-z0
  call kh(n,u,hx,w)
  tmp1=w*u
  tmp2=tmp1*u
  a=sum(w)
  b=sum(tmp1)
  d=sum(tmp2)
  deno=a*d-b**2.
  do j=1,p
     tmp1=w*y(:,j)
     tmp2=tmp1*u
     t1=sum(tmp1)
     t2=sum(tmp2)
     num1=d*t1-b*t2
     num2=-b*t1+a*t2
     f(j)=num1/(deno+eps)
  end do
  return
end subroutine locallinear0

!local linear estimation of p nonparametric functions and their derivatives with respect to covariate, evaluated at x0, difference bandwidth for f and f'
subroutine locallinear(n,p,z0,hx,hy,z,y,f,df)
  integer,intent(in) :: n,p
  double precision,intent(in) :: z0,hx,hy
  double precision,dimension(n),intent(in) :: z 
  double precision,dimension(n,p),intent(in) :: y
  double precision,dimension(p),intent(out) :: f,df
  double precision,dimension(n) :: u,w,tmp1,tmp2
  double precision :: a,b,d,deno,t1,t2,num1,num2,eps
  integer :: i,j
  eps=1e-8
  u=z-z0
! get function
  call kh(n,u,hx,w)
  tmp1=w*u
  tmp2=tmp1*u
  a=sum(w)
  b=sum(tmp1)
  d=sum(tmp2)
  deno=a*d-b**2.
  do j=1,p
     tmp1=w*y(:,j)
     tmp2=tmp1*u
     t1=sum(tmp1)
     t2=sum(tmp2)
     num1=d*t1-b*t2
     num2=-b*t1+a*t2
     f(j)=num1/(deno+eps)
!     df(j)=num2/(deno+eps)
  end do
! get derivative
  call kh(n,u,hy,w)
  tmp1=w*u
  tmp2=tmp1*u
  a=sum(w)
  b=sum(tmp1)
  d=sum(tmp2)
  deno=a*d-b**2.
  do j=1,p
     tmp1=w*y(:,j)
     tmp2=tmp1*u
     t1=sum(tmp1)
     t2=sum(tmp2)
     num1=d*t1-b*t2
     num2=-b*t1+a*t2
!     f(j)=num1/(deno+eps)
     df(j)=num2/(deno+eps)
  end do
  return
end subroutine locallinear

!quatic kernel function kh evaluated at x (n vector), bandwidth h
subroutine kh(n,x,h,y)
  integer,intent(in) :: n
  double precision,dimension(n),intent(in) :: x
  double precision,intent(in) :: h
  double precision,dimension(n),intent(out) :: y
  double precision,dimension(n) :: tmp
  integer :: i
  tmp=x/h
  y=0
  do i=1,n
     if (abs(tmp(i)).lt.1) then
        y(i)=(1.-tmp(i)**2.)**2./h*15./16.
     end if
  end do
  return
end subroutine kh

subroutine kh2(n,x,h,y)
  use commondata, only: pi 
  integer,intent(in) :: n
  double precision,dimension(n),intent(in) :: x
  double precision,intent(in) :: h
  double precision,dimension(n),intent(out) :: y
  double precision,dimension(n) :: tmp
  integer :: i
  tmp=x/h
  y=0
  do i=1,n
     if (abs(tmp(i)).lt.1) then
        y(i)=1./sqrt(2.*pi)*exp(-1./2.*(tmp(i)**2))/h  ! Gauss Kernel 
!        y(i)=(1.-tmp(i)**2.)**2./h*15./16. !Quartic
!        y(i)=(1.-tmp(i)**2.)/h*3./4. !Epanechnikov 
!        y(i)=pi/4.*cos(pi*tmp(i)/2.)/h ! cosine Kernel
!        y(i)=(1-tmp(i)**2)**3*35./32./h ! Triweight
!        y(i)=70./81./h*(1-abs(tmp(i))**3)**3
     end if
  end do
  return
end subroutine kh2

!variance calculation of the god' estimator
subroutine covcsgod(betaest,covest) !,A,B)
  use commondata, only: pt,n,delta
  double precision,dimension(pt),intent(in) :: betaest
  double precision,dimension(pt,pt),intent(out) :: covest
  double precision,dimension(n) :: xb,fb
  double precision,dimension(pt) :: x0,betal,betar,fl,fr
  double precision,dimension(n,pt) :: fval
  double precision,dimension(pt,pt) :: B,A,Ainv
  integer :: i,j
  call csgod(betaest,fval)
  B=matmul(transpose(fval),fval)/n
  do j=1,pt
     betal=betaest
     betal(j)=betaest(j)*(1-delta)
     call csgod(betal,fval)
     fl=sum(fval,1)/n
     betar=betaest
     betar(j)=betaest(j)*(1+delta)
     call csgod(betar,fval)
     fr=sum(fval,1)/n
     A(:,j)=(fr-fl)/2/delta/betaest(j)
  end do
  call inv(pt,A,Ainv)
  covest=matmul(matmul(Ainv,B),transpose(Ainv))
  return
end subroutine covcsgod

!variance calculation of the local estimator
subroutine covcslocal(betaest,covest)
  use commondata, only: pt,n,delta
  double precision,dimension(pt),intent(in) :: betaest
  double precision,dimension(pt,pt),intent(out) :: covest
  double precision,dimension(n) :: zb,fb
  double precision,dimension(pt) :: z0,betal,betar,fl,fr
  double precision,dimension(n,pt) :: fval
  double precision,dimension(pt,pt) :: B,A,Ainv
  integer :: i,j
  call cslocal(betaest,fval)
  B=matmul(transpose(fval),fval)/n
  do j=1,pt
     betal=betaest
     betal(j)=betaest(j)*(1-delta)
     call cslocal(betal,fval)
     fl=sum(fval,1)/n
     betar=betaest
     betar(j)=betaest(j)*(1+delta)
     call cslocal(betar,fval)
     fr=sum(fval,1)/n
     A(:,j)=(fr-fl)/2/delta/betaest(j)
  end do
  call inv(pt,A,Ainv)
  covest=matmul(matmul(Ainv,B),transpose(Ainv))  
  return
end subroutine covcslocal

!variance calculation of the efficient estimator
subroutine covcseff(betaest,covest)
  use commondata, only: pt,n,delta
  double precision,dimension(pt),intent(in) :: betaest
  double precision,dimension(pt,pt),intent(out) :: covest
  double precision,dimension(n) :: xb,fb
  double precision,dimension(pt) :: x0,betal,betar,fl,fr
  double precision,dimension(n,pt) :: fval
  double precision,dimension(pt,pt) :: B,A,Ainv
  integer :: i,j
  call cseff(betaest,fval)
  B=matmul(transpose(fval),fval)/n
  do j=1,pt
     betal=betaest
     betal(j)=betaest(j)*(1-delta)
     call cseff(betal,fval)
     fl=sum(fval,1)/n
     betar=betaest
     betar(j)=betaest(j)*(1+delta)
     call cseff(betar,fval)
     fr=sum(fval,1)/n
     A(:,j)=(fr-fl)/2/delta/betaest(j)
  end do
  call inv(pt,A,Ainv)
  covest=matmul(matmul(Ainv,B),transpose(Ainv))  
  return
end subroutine covcseff

subroutine compositesimpson(g,h,integral)
! Integration of f(x) on [0,1]
! Method: Composite Simpson rule for q intervals  
! g,h - Product Functions to integrate 
! q   - number of intervals
! integral - Result of integration

use commondata, only:q
double precision,dimension(q+1):: g,h 
double precision integral,s,t 
integer i

t = 1/dble(q)
s = 0.0
do i=1, q/2
   s = s + g(2*i-1)*h(2*i-1) + 4.0*g(2*i)*h(2*i) + g(2*i+1)*h(2*i+1) 
end do
integral = t/3.0*s 
return
end subroutine compositesimpson

subroutine splinebasis(d, q, m1, m, s, t, innerknots, boundaryknots, basis,ic)
! This subroutine generates Bspline basis functions.
! t(q) is a q by 1 input vector for which B-spline basis
! function() will be evaluated.
! innerknots(m1) set of m1 innerknot points.
! newknots is the entire set of knots, of length m=m1+2(d+1)
! where d is the degree of the splines. 
! order is d+1
! s=number of spline basis=m1+d+1
! ic=1 means there is no Intercept 

integer :: s,q,m1,d  
double precision t(q), innerknots(m1), boundaryknots(2)
double precision newknots(m), basis(q,s), result
integer :: i1, i, j, l 

do i1=1, (d+1)
  newknots(i1)=boundaryknots(1)
end do
do i1=(d+2), (m1+d+1)
  newknots(i1)=innerknots(i1-d-1)
end do
do i1=(m1+d+2), m
  newknots(i1)=boundaryknots(2) 

end do

do i=1, q
  if(t(i).eq.boundaryknots(2)) then
    basis(i, s)=1.d0
       do j=1, (s-1)
          basis(i, j)=0.d0
       end do
  else
        do j=1, s
          call b(m, j, (d+1), t(i), newknots, result, b)
          basis(i, j)=result  
        end do
  endif
end do

if (ic.eq.1) then
  basis=basis(1:q,2:s)
end if 

return
end subroutine splinebasis


subroutine b(i1, i2, i3, y, newknots, result, dumsub)
! This subroutine calculates i2 th basis of spline of degree (i3-1).

implicit none
integer :: i1, i2, i3
double precision y, newknots(i1), temp1, temp2, result, result1, result2
external dumsub
if(i3.eq.1) then
 if((y.ge.newknots(i2)).and.(y.lt.newknots(i2+1))) then
  result=1.d0
 else
  result=0.d0
 end if 
end if 

if(i3.ne.1) then 
 call dumsub(i1, i2, (i3-1), y, newknots, result1, dumsub)
 temp1=(y-newknots(i2))*result1/(newknots(i2+i3-1)-newknots(i2))
 if(temp1.ne.temp1) then 
  temp1=0.d0
 end if
 
 call dumsub(i1, (i2+1), (i3-1), y, newknots, result2, dumsub)
 temp2=(newknots(i2+i3)-y)*result2/(newknots(i2+i3)-newknots(i2+1))
 if(temp2.ne.temp2) then
  temp2=0.d0
 end if 
 result=temp1+temp2
end if
return
 
end subroutine b 


! inner product of two vectors: y=s't
subroutine inprod(l,s,t,y)
  integer,intent(in) :: l
  double precision,dimension(l),intent(in) :: s,t
  double precision,intent(out) :: y
  integer :: i
  y=0
  do i=1,l
     y=y+s(i)*t(i)
  end do
  return
end subroutine inprod

subroutine regress2(n,p1,p2,x0,h,x,z,my)
  integer,intent(in) :: n,p1,p2 
  double precision,intent(in) :: x0,h
  double precision,dimension(n),intent(in) :: x
  double precision,dimension(n,p1,p2),intent(in) :: z
  double precision,dimension(p1,p2),intent(out) :: my    
  double precision,dimension(n) :: t,s
  double precision,dimension(p1,p2) :: t2
  double precision :: t1
  integer :: i
  t=x-x0
  call kh(n,t,h,s)
  t1=sum(s)
  t2=0
  do i=1,n
     t2=t2+z(i,:,:)*s(i)
  end do
  my=t2/t1
  return
end subroutine regress2



