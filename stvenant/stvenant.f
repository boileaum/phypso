c234567
c resolution par volumes finis des eq de st venant
c ce programme resout le pb de riemann
c a completer...

      program venant

      implicit double precision (a-h,o-z)

c vitesses min et max pour le trace
      ximin=-6
      ximax=6

c nombre de points
      nmax=1000
      dv=(ximax-ximin)/nmax

c test du solveur de riemann
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

c a vous de jouer...

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fonction Z(h1,h2)

      function z(h1,h2)
      implicit double precision (a-h,o-z)

      g=9.81

      if (h1.gt.h2) then
       zsol=sqrt(2.d0)*sqrt(g*(h1+h2)/h1/h2)/2.d0
      endif

      if (h1.le.h2) then
       zsol=2.d0*dsqrt(g)/(sqrt(h1)+sqrt(h2))
      endif

      z=zsol
      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c derivee de z par rapport a h1

      function dzdh1(h1,h2)
      implicit double precision (a-h,o-z)

      g=9.81

      if (h1.gt.h2) then
      dzsol= sqrt(2.d0)/sqrt(g*(h1+h2)/h1/h2)*
     #(g/h1/h2-g*(h1+h2)/h1**2/h2)/4.d0
      endif

      if (h1.le.h2) then
       dzsol=-2.d0*sqrt(g)/(sqrt(h1)+sqrt(h2))**2/sqrt(h1)/2
      endif

      dzdh1=dzsol


      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c resolution du pb de riemann

      subroutine riemann(hg,ug,hd,ud,xi,h,u)
      implicit double precision (a-h,o-z)

c accel pesanteur
      g=9.81d0

c valeur initiale de hstar pour newton
      hstar=.5*(hg+hd)

c debut des iterations
 1    continue

c calcul du residu (on veut que res=0)
      res= ud-ug+(hstar-hd)*z(hstar,hd)
     #     + (hstar-hg)*z(hstar,hg)

c calcul de l'accroissement
      dh=z(hstar,hd) +(hstar-hd)*dzdh1(hstar,hd)
     # + z(hstar,hg) +(hstar-hg)*dzdh1(hstar,hg) 

      dh=-res/dh

c mise a jour
      hstar=hstar+dh

c test d'arret
      if ((dabs(dh).gt.1.d-10).and.(dabs(res).gt.1.d-10)) then 
      goto 1
      endif

c      write(*,*) hstar

      ustar=ud+(hstar-hd)*z(hstar,hd)
c      write(*,*) ustar,ug-(hstar-hg)*z(hstar,hg),res

c celerite du "son" a gauche droite et milieu
      cg=dsqrt(g*hg)
      cd=dsqrt(g*hd)
      cstar=dsqrt(g*hstar)

c calcul de lambda(i,plus) et lambda(i,moins)
c onde 1
      if (hstar.gt.hg) then
c choc
         xlambda1m=(hg*ug-hstar*ustar)/(hg-hstar)
         xlambda1p=xlambda1m
      else
c detente
         xlambda1m=ug-cg
         xlambda1p=ustar-cstar
      endif

c onde 2
      if (hstar.gt.hd) then
c choc
         xlambda2m=(hd*ud-hstar*ustar)/(hd-hstar)
         xlambda2p=xlambda2m
      else
c detente
         xlambda2p=ud+cd
         xlambda2m=ustar+cstar
      endif

c      write(*,*) xlambda1m,xlambda1p,xlambda2m,xlambda2p

c en fonction de xi on evalue h et u

      if (xi.lt.xlambda1m) then
         h=hg
         u=ug
      endif

      if ((xi.ge.xlambda1m).and.(xi.lt.xlambda1p)) then
c 1-detente
         h=(2.d0*dsqrt(hg)/3.d0+(ug-xi)/3.d0/dsqrt(g))**2;
         u=ug-(h-hg)*z(h,hg)
      endif
      
      if ((xi.ge.xlambda1p).and.(xi.lt.xlambda2m)) then
c etat milieu
         h=hstar
         u=ustar
      endif

      if ((xi.ge.xlambda2m).and.(xi.lt.xlambda2p)) then
c 2-detente
         h=(2.d0*dsqrt(hd)/3.d0+(xi-ud)/3.d0/dsqrt(g))**2;
         u=ud+(h-hd)*z(h,hd)
      endif

      if (xi.ge.xlambda2p) then
         h=hd
         u=ud
      endif

      return
 
      end


      

