!     resolution par volumes finis des eq de st venant
!     ce programme resout le pb de riemann
!     a completer...

program venant

  implicit none
  real*8 :: ximin,ximax,dv,xi,h,u
  integer :: nmax,i

  !     vitesses min et max pour le trace
  ximin=-6
  ximax=6

  !     nombre de points
  nmax=1000
  dv=(ximax-ximin)/nmax

  !     test du solveur de riemann
  open(1,file='h')
  open(2,file='u')

  do i=0,nmax
     xi = ximin+i*dv
     call riemann(2.d0,0.d0,1.d0,0.d0,xi,h,u)
     write(1,*) xi,h
     write(2,*) xi,u
  enddo

  close(1)
  close(2)

  !     a vous de jouer...

end program venant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cccccc
!     fonction Z(h1,h2)

function z(h1,h2)
  implicit none
  real*8 :: g,h1,h2,z,zsol

  g=9.81

  if (h1.gt.h2) then
     zsol=sqrt(2.d0)*sqrt(g*(h1+h2)/h1/h2)/2.d0
  endif

  if (h1.le.h2) then
     zsol=2.d0*dsqrt(g)/(sqrt(h1)+sqrt(h2))
  endif

  z=zsol

end function z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cccccccc
!     derivee de z par rapport a h1

function dzdh1(h1,h2)
  implicit none
  real*8 :: g,zsol,h1,h2,dzdh1,dzsol

  g=9.81

  if (h1.gt.h2) then
     dzsol= sqrt(2.d0)/sqrt(g*(h1+h2)/h1/h2)* &
          (g/h1/h2-g*(h1+h2)/h1**2/h2)/4.d0
  endif

  if (h1.le.h2) then
     dzsol=-2.d0*sqrt(g)/(sqrt(h1)+sqrt(h2))**2/sqrt(h1)/2
  endif

  dzdh1=dzsol


end function dzdh1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ycccccccc
!     resolution du pb de riemann

subroutine riemann(hg,ug,hd,ud,xi,h,u)
  implicit none
  real*8 :: g,hg,hd,xi,h,u,hstar,dh,res
  real*8 :: cg,cd,cstar,ustar,ug,ud
  real*8:: xlambda1m,xlambda1p,xlambda2m,xlambda2p
  real*8 :: z,dzdh1

  !     accel pesanteur
  g=9.81d0

  !     valeur initiale de hstar pour newton
  hstar=.5*(hg+hd)

  !     debut des iterations
1 continue

  res=1

  do while((dabs(dh).gt.1.d-10).or.(dabs(res).gt.1.d-10))

     !     calcul du residu (on veut que res=0)
     res= ud-ug+(hstar-hd)*z(hstar,hd) &
          + (hstar-hg)*z(hstar,hg)

     !     calcul de l'accroissement
     dh=z(hstar,hd) +(hstar-hd)*dzdh1(hstar,hd) &
          + z(hstar,hg) +(hstar-hg)*dzdh1(hstar,hg) 

     dh=-res/dh

     !     mise a jour
     hstar=hstar+dh

     !     test d'arret
  enddo
  !     write(*,*) hstar

  ustar=ud+(hstar-hd)*z(hstar,hd)
  !     write(*,*) ustar,ug-(hstar-hg)*z(hstar,hg),res

  !     celerite du "son" a gauche droite et milieu
  cg=dsqrt(g*hg)
  cd=dsqrt(g*hd)
  cstar=dsqrt(g*hstar)

  !     calcul de lambda(i,plus) et lambda(i,moins)
  !     onde 1
  if (hstar.gt.hg) then
     !     choc
     xlambda1m=(hg*ug-hstar*ustar)/(hg-hstar)
     xlambda1p=xlambda1m
  else
     !     detente
     xlambda1m=ug-cg
     xlambda1p=ustar-cstar
  endif

  !     onde 2
  if (hstar.gt.hd) then
     !     choc
     xlambda2m=(hd*ud-hstar*ustar)/(hd-hstar)
     xlambda2p=xlambda2m
  else
     !     detente
     xlambda2p=ud+cd
     xlambda2m=ustar+cstar
  endif

  !     write(*,*) xlambda1m,xlambda1p,xlambda2m,xlambda2p

  !     en fonction de xi on evalue h et u

  if (xi.lt.xlambda1m) then
     h=hg
     u=ug
  endif

  if ((xi.ge.xlambda1m).and.(xi.lt.xlambda1p)) then
     !     1-detente
     h=(2.d0*dsqrt(hg)/3.d0+(ug-xi)/3.d0/dsqrt(g))**2;
     u=ug-(h-hg)*z(h,hg)
  endif

  if ((xi.ge.xlambda1p).and.(xi.lt.xlambda2m)) then
     !     etat milieu
     h=hstar
     u=ustar
  endif

  if ((xi.ge.xlambda2m).and.(xi.lt.xlambda2p)) then
     !     2-detente
     h=(2.d0*dsqrt(hd)/3.d0+(xi-ud)/3.d0/dsqrt(g))**2;
     u=ud+(h-hd)*z(h,hd)
  endif

  if (xi.ge.xlambda2p) then
     h=hd
     u=ud
  endif

  return

end subroutine riemann




