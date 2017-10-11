c234567
      program burgers

c ce programme calcule la solution de l'équation de burgers
c par la méthode des volumes finis en 1D
c sur un intervalle [a1,a2]
c pour le tester, on compare à une solution exacte

c les variables commençant par a-h ou o-z sont des réels
c ATTENTION: les autres sont des entiers !
      implicit real*8 (a-h,o-z)

c en FORTRAN 77 on n'a pas le droit de dimensionner
c un tableau avec des variables qui ne sont pas des
c "parameter".
c les lignes du programme ne peuvent pas faire plus de 
c 72 colonnes. Elles doivent commencer en colonne >  7
c Une ligne de commentaire commence par un "c" en colonne 1
c Si une ligne est trop longue on continue sur la ligne 
c d'après en mettant un caractère quelconque en colonne 6
      parameter(nmax=100)

c solution numérique aux pas n et n+1
c et exemple de ligne coupée !
c234567
      dimension wn(0:nmax+1),
     $     wnp1(0:nmax+1)
c milieux des cellules
      dimension xmil(0:nmax+1)

c initialisations
c bornes en espace
      a1=-1.d0
      a2=2.d0
c borne en temps
      write(*,*) 'entrer le temps final'
      read(*,*) tmax
c pas d'espace
      dx=(a2-a1)/nmax
c coefficient de CFL
      cfl=0.8d0

c vecteur initiaux
c on triche en remplaçant la valeur moyenne
c par la valeur au centre de la cellule
c on initialise aussi les deux cellules fictives
c i=0 et i=nmax+1
      do i=0,nmax+1
         xmil(i)=a1+(i-0.5d0)*dx
         call solexacte(xmil(i),0.d0,w)
         wn(i)=w
         wnp1(i)=wn(i)
      enddo

      t=0.d0
      do while(t.lt.tmax)
c recherche de la vitesse max
         vmax=0.d0
         do i=0,nmax+1
            if (dabs(wn(i)).gt.vmax) vmax=dabs(wn(i))
         enddo
         dt=cfl*dx/vmax
c boucle sur les cellules internes i de 1 à nmax
         do i=1,nmax
c flux à droite
            call riemann(wn(i),wn(i+1),0.d0,w)
            wnp1(i)=wnp1(i)-dt/dx*w*w*0.5d0
c flux à gauche
            call riemann(wn(i-1),wn(i),0.d0,w)
            wnp1(i)=wnp1(i)+dt/dx*w*w*0.5d0
         enddo
c actualisations
         do i=1,nmax
            wn(i)=wnp1(i)
         enddo
         t=t+dt
         write(*,*) 't=',t
      enddo

c écriture de la solution exacte et de la
c solution numérique dans 2 fichiers "exact" et "godu"
c pour tracé dans "gnuplot"
c taper "gnuplot" puis, dans gnuplot
c "  plot 'godu', 'exact'  "

      open(1,file='godu')
      open(2,file='exact')

      do i=1,nmax
         write(1,*) xmil(i),wn(i)
         call solexacte(xmil(i),t,w)
         write(2,*) xmil(i),w
      enddo

      close(1)
      close(2)

      end


c ce sous-programme resout le problème de riemann
c234567
      subroutine riemann(wl,wr,xi,w)
      implicit real*8 (a-h,o-z)

c cas d'un choc
      if (wl.gt.wr) then

c vitesse du choc
         sigma=0.5d0*(wl+wr)

         if (xi.lt.sigma) w=wl
         if (xi.ge.sigma) w=wr

      else

         if (xi.le.wl) w=wl
         if (xi.ge.wr) w=wr
         if (xi.gt.wl.and.xi.lt.wr) w=xi

      endif

      return
      end
         

c ce sous-programme calcule une solution exacte
c de l'équation de burgers (afin de vérifier)
c234567
      subroutine solexacte(x,t,w)

      implicit real*8 (a-h,o-z)

      if (t.le.1.d0) then

         if (x.le.0.d0) w=1.d0
         if (x.ge.1.d0) w=0.d0
         if (x.gt.t.and.x.lt.1) w=(1.d0-x)/(1.d0-t)

      else

         if ((x-1.d0).le.0.5d0*(t-1.d0)) then
            w=1.d0
         else
            w=0.d0
         endif
         
      endif

      return
      end

