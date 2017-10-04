! programme à compléter
! chercher les triples "!!!" et modifier la partie
! du programme correspondante


! module fortran pour gérer un maillage volumes finis
module vf

  implicit none 

  ! nx: nombre de vf horizontaux
  ! ny: nombre de vf verticaux
  integer  :: nx,ny


  ! nbvol: nombre de vf
  ! nbsom: nombre de sommets
  ! nbaret: nombre d'arêtes
  integer :: nbvol,nbsom,nbaret
  ! xmin,xmax:dimensions en x
  ! ymin,ymin:dimensions en y
  real*8   :: xmin,xmax,ymin,ymax

  ! tmax: durée de la simulation
  real*8   :: tmax

  ! cfl: cfl
  ! g: accélération de la pesanteur
  real*8  :: cfl=0.8d0
  real*8,parameter :: g=9.81d0

  ! imiroir = -1 --> obstacle
  ! imiroir =  1 --> on recopie (debugage)
  integer :: miroir=-1


  ! tableaux
  ! xs,ys: coordonnées des sommets
  ! xv,yv  : milieux des volumes
  ! topo   ::  altitude des milieux des volumes
  real*8,dimension(:),allocatable  :: xs,ys,xv,yv,topo

  ! iaret: tableaux des arêtes
  ! iaret(1,ia) et le n° du premier volume voisin de ia
  ! iaret(2,ia) et le n° du second volume voisin de ia
  ! si iaret(2,ia)=0, ia est une arête de bord
  integer,dimension(:,:),allocatable :: iaret

  ! iext: permet de construire le maillage (inutile après)
  ! along, xnorm : tableaux des longueurs et coordonnées
  ! des  normales aux arêtes
  real*8,dimension(:),allocatable  :: along
  real*8,dimension(:,:),allocatable :: xnorm  ! xnorm(1:2,1:nbaret)

  ! ivol(ii,iv): n° du sommet ii (ii=1..4) du volume iv
  integer,dimension(:,:),allocatable   :: ivol

  ! vol(iv): volume du vf n° iv
  ! perim(iv): périmètre du vf n° iv
  real*8,dimension(:),allocatable  :: vol,perim

  ! wn et wnp1: stockage de w à l'étape n et n+1
  integer   :: mc=3 ! nb de variables conservées (3 pour shallow water)
  real*8,dimension(:,:),allocatable :: wn,wnp1


contains

  ! sous-programme de construction du maillage
  subroutine init_vf(px,py,xmi,xma,ymi,yma)

    integer,intent(in)  :: px,py
    integer   :: i,j,is,ia,icomptar,iloc1,iloc2,is1,is2,iv
    integer  :: init=0,ignuplot=1,ii
    real*8 :: dx,dxx,dlong,dy,dyy,x1,x2,xmil,xn,y1,y2,ymil,yn,coef
    real*8  :: xmi,xma,ymi,yma
    integer,dimension(:,:),allocatable :: iext

    nx=px
    ny=py

    xmin=xmi
    xmax=xma
    ymin=ymi
    ymax=yma

    write(*,*) 'xmin=',xmin,' xmax=',xmax
    write(*,*) 'ymin=',ymin,' ymax=',ymax

    nbvol=nx*ny   ! nombre de vf
    nbsom=(nx+1)*(ny+1) ! nombre de sommets
    nbaret=nx*(ny+1)+(nx+1)*ny  ! nombre d'arêtes

! tableaux du mmodule
    allocate(xs(nbsom))
    allocate(ys(nbsom))
    allocate(xv(nbvol))
    allocate(yv(nbvol))
    allocate(iaret(2,nbaret))
    allocate(ivol(4,nbvol))
    allocate(along(nbaret))
    allocate(xnorm(2,nbaret))
    allocate(vol(nbvol))
    allocate(wn(2,nbvol))
    allocate(wnp1(2,nbvol))
    allocate(perim(nbvol))
    allocate(topo(nbvol))

! tableau de travail
    allocate(iext(2,nbaret))

    !     variables de base

    write(*,*) 'nx=',nx,' ny=',ny
    write(*,*) 'nb de volumes=',nbvol
    write(*,*) 'nb de sommets=',nbsom
    write(*,*) 'nb d''arêtes=',nbaret

    !     pas en x et en y

    dx=(xmax-xmin)/nx
    dy=(ymax-ymin)/ny

    !     remplissage des tableaux de connectivité

    write(*,*) 'construction du maillage'

    !     les sommets

    do i=0,nx
       do j=0,ny
          is=i+j*(nx+1)+1
          xs(is)=i*dx+xmin
          ys(is)=j*dy+ymin
       end do
    end do


    !     les volumes
    xv=0
    yv=0

    do i=0,nx-1
       do j=0,ny-1
          iv=i+j*nx+1
! connectivité volumes --> sommets
          ivol(1,iv)=i+j*(nx+1)+1
          ivol(2,iv)=i+1+j*(nx+1)+1
          ivol(3,iv)=i+1+(j+1)*(nx+1)+1
          ivol(4,iv)=i+(j+1)*(nx+1)+1
          ! coordonnées des milieux des mailles
          do ii=1,4
             xv(iv)=xv(iv)+xs(ivol(ii,iv))/4
             yv(iv)=yv(iv)+ys(ivol(ii,iv))/4
          end do
          topo(iv) =0     !!!  à changer pour le fond variable
! volume et périmètre
          vol(iv)=dx*dy
          perim(iv)=2*dx+2*dy
       end do
    end do

    !     les arêtes horizontales

    icomptar=0

    do i=0,nx-1
       do j=0,ny
          ia=i+j*nx+1
! connectivité arêtes --> volumes
          iaret(1,ia)=i+j*nx+1
          iaret(2,ia)=i+(j-1)*nx+1
          iext(1,ia)=1
          iext(2,ia)=3

          !     correction en haut et en bas pour que
          !     le vecteur normal soit toujours sortant
          if (j.eq.ny) then
             iaret(1,ia)=i+(j-1)*nx+1
             iext(1,ia)=3
             iext(2,ia)=0
             iaret(2,ia)=0
          endif
          if (j.eq.0) then
             iaret(2,ia)=0
             iext(2,ia)=0
          endif

          icomptar=icomptar+1
       end do
    end do

    !     les arêtes verticales

    do i=0,nx
       do j=0,ny-1
! connectivité arêtes --> volumes
          ia=i+j*(nx+1)+1+icomptar
          iaret(1,ia)=i-1+j*nx+1
          iaret(2,ia)=i+j*nx+1
          iext(1,ia)=2
          iext(2,ia)=4
          !     correction à gauche et à droite pour que
          !     le vecteur normal soit toujours sortant
          if (i.eq.0) then
             iaret(1,ia)=i+j*nx+1
             iext(1,ia)=4
             iext(2,ia)=0
             iaret(2,ia)=0
          endif
          if (i.eq.nx) then
             iaret(2,ia)=0
             iext(2,ia)=0
          endif
       end do
    end do


    !     calcul des variables d'arêtes:
    !     longueur, vecteurs normaux
    do ia=1,nbaret
       iv=iaret(1,ia)
       iloc2=iext(1,ia)
       iloc1=iloc2+1
       if(iloc1.eq.5) iloc1=1
       is1=ivol(iloc1,iv)
       is2=ivol(iloc2,iv)
       x1=xs(is1)
       y1=ys(is1)
       x2=xs(is2)
       y2=ys(is2)
       xmil=0.5*(x1+x2)
       ymil=0.5*(y1+y2)
       dlong=dsqrt((x1-x2)**2+(y1-y2)**2)
       along(ia)=dlong
       dxx=(x2-x1)/dlong
       dyy=(y2-y1)/dlong
       xn=-dyy
       yn=dxx
       xnorm(1,ia)=xn
       xnorm(2,ia)=yn
    enddo

    !    dessin du maillage
    !     pour visualiser le maillage et les arêtes dans gnuplot
    !     taper: plot 'mail' with line
    if (init.eq.0.and.ignuplot.eq.1) then
       open(2,file='mail.plt')
       write(2,*) 'plot ''mail'' w l'
       close(2)
       open(1,file='mail')
       do iv=1,nbvol
          write(1,*) xs(ivol(1,iv)),ys(ivol(1,iv))
          write(1,*) xs(ivol(2,iv)),ys(ivol(2,iv))
          write(1,*) xs(ivol(3,iv)),ys(ivol(3,iv))
          write(1,*) xs(ivol(4,iv)),ys(ivol(4,iv))
          write(1,*) xs(ivol(1,iv)),ys(ivol(1,iv))
          write(1,*)
       end do
       do ia=1,nbaret
          iv=iaret(1,ia)
          iloc2=iext(1,ia)
          iloc1=iloc2+1
          if(iloc1.eq.5) iloc1=1
          is1=ivol(iloc1,iv)
          is2=ivol(iloc2,iv)
          x1=xs(is1)
          y1=ys(is1)
          x2=xs(is2)
          y2=ys(is2)
          xmil=0.5*(x1+x2)
          ymil=0.5*(y1+y2)
          coef=.25*dsqrt(dx*dx+dy*dy)*.5
          write(1,*) xmil,ymil
          write(1,*) xmil+xnorm(1,ia)*coef,ymil+xnorm(2,ia)*coef
          write(1,*)
       enddo
       close(1)
       call system('gnuplot mail.plt')   ! compatible uniquement avec gfortran
    endif

    deallocate(iext)

  end subroutine init_vf

! sous-programme d'initialisation des varaibles conservatives
  subroutine init_w()

    implicit none
    integer :: i,ii
    integer,parameter :: m=2

    write(*,*) 'Application de la condition initiale'

    do i=1,nbvol
       wn(1,i) = xv(i)**2+yv(i)**2    !!! changer cette valeur pour st venant
       wn(2,i) = 0    !!! changer cette valeur pour st venant
    end do

    
    

  end subroutine init_w


! sous-programme pour les affichage en 3D gnuplot
  subroutine plot_iso()

    implicit none

    integer  :: i,j,iv

    write(*,*) 'Affichage des résultats'

    open(1,file='iso_w')
    do i=0,nx-1
       do j=0,ny-1
          iv=i+j*nx+1
          write(1,*) xv(iv),yv(iv),wn(1,iv)
       end do
       write(1,*)
    end do
    close(1)

    open(2,file='iso_w.plt')

    write(2,*) 'set hidden3d'
    write(2,*) 'set xlabel "x"'
    write(2,*) 'set ylabel "y"'
    write(2,*) 'set zlabel "z"'
    write(2,*) 'set grid xtics ytics ztics'
    write(2,*) 'splot ''iso_w'' w l'

    close(2)

    call system('gnuplot iso_w.plt')  ! compatible uniquement avec gfortran


  end subroutine plot_iso

end module vf

! fin de la définition du module

! sous-programme de calcul du flux numérique

subroutine fluxnum(wL,wR,xnorm,flux)

  implicit none
  integer,parameter   :: m=2
  real*8   :: wR(m),wL(m),xnorm(m),flux(m)


  flux=0    !!! à modifier bien sûr

end subroutine fluxnum

! programme principal de calcul du schéma volumes finis
program shallow2d
  
  use vf
  
  implicit none
  
  real*8    :: t,vmax,dt
  integer  :: iv
  
  ! construction du maillage
  call init_vf(20,20,0.d0,1.d0,0.d0,1.d0)
  
  ! remplissage de la condition initiale
  call init_w()
  
  t=0
  tmax=0
  cfl=0.8
  
  do while(t < tmax) 
     
! calcul de la vitesse d'onde max
!!! à compléter
     vmax=0
     do iv=1,nbvol
        vmax=max(vmax,1.d0)    !!! à modifier
     end do

     dt=vol(iv)/perim(iv)/vmax*cfl

!!!    etc.
!!!   à compléter
!!! copier wn dans wnp1
!!! boucle sur les arêtes
!!! calcul des flux et envoi dans les volumes voisins
!!! pour les arêtes frontière extrapolation de wR pour le calcul du flux
!!! boucle sur les volumes pour le terme source dû à la topographie


  end do


! sauvegarde et affichage du résultat
  call plot_iso()


end program shallow2d



