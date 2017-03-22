
!PDM 2016
!This program computes the noise power spectrum for a
!cylindrical array like CHIME. experimental details can
!be specified below. Cosmology only affects projection and
!is limited to a flat Universe with a cosmological constant and
!matter as the two main components
!The code closely follows arXiv: 1410.7794
!!!!!!!!!!!
!The hardest part of the computation is the loop over
!unique baselines and this has been openmp parallized. The larger the number
!of antenna, the longer this computation takes depending on how well you want
!to sample the visability function. Current settings are a little overkill
!and the computation takes about 4 minutes on 16 cores.
!!!!!!!!!!!
!the senseitivity along the length of the cylinder (L) is now assumed to
!be reflection with a dish set by the length of the cylinder devided
!by the # of antenna. This is not correct, since the surface of the cylinder is
!only curved in the (W) direction. If you want to be more realistic you
!should lower this number by roughlt the wavelength of 21cm (i.e. 21cm :-))
!This then will require higher resolution in the sampling of the visability function
!since, if the # of anteanna is not closely packed, will lead to 'blind' spots
!in the array. If you plot the visablity function, you will see little holes in it.
!This ofcourse will result in a more bumpy noise curve. I found though that when you average
!(which is what I do anyways) this is not a big effect. Currently, I actually
!did set it to 21cm (see below). 
!!!!!!!!!!
!For the experiment you also have to specify the resolition in k-space and in redshift space.
!Anyways, modification should be relativily straightforward
!!!!!!!!!
!output of the code:
!----normalized visability function in the U_perp plane (or U-V plane)
!----Noise matrix in redshift-k space in mK^2 Mpc^3.
!Only use the noise for spectra with the SAME res
!settings. Otherwise you will be undercountin/overcounting information. 


program visability

  implicit none

  integer :: i,j,l,t
  integer, parameter :: dl= KIND(1.d0)

  !experiment:
  real(dl), parameter :: WidthCyl = 20.d0 !m
  real(dl), parameter :: LengthCyl = 100.d0 !m
  integer, parameter :: Nantenna = 256 !# of antenna p cylinder
  integer, parameter :: Ncyl = 4
  real(dl), parameter :: years = 1
  real(dl), parameter :: fsky = 0.3
  real(dl) :: zmin  = 0.7


  !cosmology:
  real(dl), parameter :: Omegam = 0.28
  real(dl), parameter :: h0 = 0.68

  !params
  real(dl), parameter :: l21 = 0.21d0 ! meters
  real(dl), parameter :: pi = 3.14159265359
  real(dl) :: dk = 0.005/h0 !in Mpc^-1
  real(dl) :: dz = 0.5d0

  real(dl) :: kmax, zmax  
  integer :: nbase

  !some arrays:
  integer :: Nk = 100
  real(dl), allocatable :: kloglist(:)
  real(dl), allocatable :: klinlist(:)
  integer :: NLambda = 5
  real(dl), allocatable :: llist(:)
  integer :: Nphi = 200
  integer :: Nmu  = 100
  real(dl), allocatable :: mulist(:)
  real(dl), allocatable :: philist(:)


  real(dl), allocatable :: UniquePairs(:,:)
  real(dl), allocatable :: MyGrid(:,:)
  real(dl), allocatable :: VisFunction(:,:)
  integer :: Ndetect
  integer :: s
  real(dl) :: DeltauL, DeltauW, uW0, uL0
  real(dl) :: lab, uL, uW, uLstep, uWstep, labstep
  integer :: Wsteps, Lsteps
  real(dl) :: res, sumX
  real(dl) :: tempNorm, Norm
  real(dl) :: phi,mu
  real(dl) :: y
  real(dl) :: k
  real(dl) :: shellphi, shellmu



  real(dl) :: secyears
  real(dl) :: Nkz
  real(dl) :: InstrumentalNoise

  !2d vs 1D (2D needed for H/D_A constraints)
  !if set to false, the code produces a file the has P(k,lambda)
  !this is useful for simple forecasts of the power spectrum (you can compare your P(k,z) to this one)
  !if set to true code produces an unintegrated P(k,mu,lambda)
  !this is useful if you would like to do Fisher forecast which rely on volume and integrate over mu (or kperp and kpar)
  logical :: want2D = .True. 


  Ndetect = Nantenna*Ncyl
  nbase = Ndetect*(Ndetect-1)/2
  DeltauW = WidthCyl
  !this is only true for fully packed. Otherwise, it has to be specified
  DeltauL = LengthCyl/(Nantenna-1)

  allocate(MyGrid(Ndetect,2))
  s = 1
  do i = 1, Nantenna
     do j = 1, Ncyl
        MyGrid(s,1) = (i-1)*DeltauL
        MyGrid(s,2) = (j-1)*DeltauW
        s = s + 1
     enddo
  enddo

  allocate(UniquePairs(nbase,2))
  s = 1
  do i = 1, Ndetect-1
     do j = i+1,Ndetect
        UniquePairs(s,1) = MyGrid(j,1)-MyGrid(i,1)
        UniquePairs(s,2) = MyGrid(j,2)-MyGrid(i,2)
        s = s + 1
     enddo
  enddo

  !check that the number of unique pairs is computed correctly
  if(s-1 .ne. nbase) print*, "number of baselines not computed correctly"

  !time to deallocate
  deallocate(MyGrid)

  !this is probably sufficient:
  !uW0 = -1.d0*Ncyl*WidthCyl
  !uL0 = -LengthCyl
  !symmetry works better:
  uW0 = -LengthCyl
  uL0 = -LengthCyl

  uLstep = 1.d0
  uWstep = 1.d0

  Wsteps = int(2*abs(uW0)/uLstep)+1
  Lsteps = int(2*abs(uL0)/uWstep)+1 

  write(*,*) "# steps in W", Wsteps
  write(*,*) "# steps in L", LSteps
  norm = 0
  tempNorm = 0
  allocate(VisFunction(Wsteps,Lsteps))

  !$OMP PARALLEL DO DEFAUlT(SHARED), &
  !$OMP PRIVATE(l,j,i,res,sumX, uL,uW, tempnorm) &
  !$OMP REDUCTION(+:norm)

  do l = 1, Wsteps
     uW = uW0 + (l-1)*Uwstep
     uL = uL0
     tempNorm = 0.d0

     do j = 1, Lsteps

        sumX = 0.d0
        res = 0.d0

        do i = 1, nbase
           !I now picked the res along the cylinder to be one wavelnegth. Can change this.
           res = Lambda((uL-UniquePairs(i,1))/(1.d0*l21))* &
                Lambda((uW-UniquePairs(i,2))/DeltauW) + &
                Lambda((uL+UniquePairs(i,1))/(1.d0*l21))* &
                Lambda((uW+UniquePairs(i,2))/DeltauW) 

           sumX = sumX + res
        enddo

        VisFunction(l,j) = sumX
        !integrating visability function
        tempNorm = tempNorm + sumX*uLstep
        uL = uL + uLstep

     enddo
     norm  = norm + uWstep*tempNorm
     write(*,*) Uw

  enddo


  !$OMP END PARALLEL DO
  write(*,*) "normalization:", norm

  !normalize visabilty function and store in file
  open(unit=21,file = 'nu_CHIMElikeExp.dat', status='replace')

  do l = 1, Wsteps
     uW = uW0 + (l-1)*Uwstep
     uL = uL0
     do j = 1, Lsteps
        write(21,"(F12.2,2X,F12.2,2X,F12.2)") uW, uL, nbase*VisFunction(l,j)/norm
        uL = uL + uLstep
     enddo
  enddo

  close(21)

  !time to deallocate
  deallocate(UniquePairs)


  allocate(kloglist(Nk))
  allocate(klinlist(Nk))
  allocate(llist(NLambda))
  allocate(philist(Nphi))
  allocate(mulist(Nmu))
  !convert U-V to k-lambda:

  zmin = 0.7d0 !also specified above: overwrite
  secyears = years*24.*60.*60.*365. !in sec

  zmax = zmin +  dz*(Nlambda-1)

  do i = 1, Nk
     !kloglist(i) = 10**(-2+(i-1)*0.004) !in Mpc^-1
     klinlist(i) = dk+dk*(i-1)*2
  enddo
  do i = 1, Nlambda
     llist(i) = l21*(1.+ zmin +(i-1)*dz) ! in m
  enddo
  do i = 1, Nphi
     philist(i) = 2.d0*pi/(Nphi-1)*(i-1)
  enddo
  do i = 1, Nmu
     mulist(i) = -1.d0+ 2.d0/(Nmu-1)*(i-1)
  enddo


  !This loop converts from W-L (or U-V) plane to k-lambda by taking shells in 3D and
  !projecting them on the 2D surface. This is done by suming over angles and deviding by the number of samples.

  if (want2D) then
     open(unit=21,file = 'Noise_CHIME_2D.dat', status='replace')
     do j = 1, Nlambda
        do i = 1, Nk
           !here you can choose if you want log or lin sampling in k
           k = klinlist(i)
           shellmu = 0.d0
           !averaging over phi/mu:
           do t = 1, Nmu 
              shellphi = 0.d0
              do l = 1, Nphi

                 uL = k*sqrt(1-mulist(t)**2)*Cos(philist(l))*llist(j)/(2.d0*pi)* DCapprox(llist(j)/L21-1.,Omegam,h0)
                 uW = k*sqrt(1-mulist(t)**2)*Sin(philist(l))*llist(j)/(2.d0*pi)* DCapprox(llist(j)/L21-1.,Omegam,h0)

                 !outside this range, we know the visability is 0. Interpolation will give out of bounds.
                 !there is also no zero mode
                 if (abs(uW) > abs(uW0)-1 .or. abs(uL) > abs(uL0)-1 .or. (abs(uL) .eq. 0.d0 .and. abs(uW) .eq. 0d0)) then
                    y = 0.d0
                 else                
                    call Simple2DLinInterp(uW,uL,uW0,uL0,uWstep,uLstep,VisFunction,y)
                 endif
                 !write(*,*) k, uL, uW, y
                 !add shell:
                 shellphi = shellphi + nbase*y/norm

              enddo !philoop
              !compute the # of modes in this shell:
              Nkz = fsky*k**2*dk* &
                   (DCapprox((llist(j)/l21 + dz/2.d0-1),OmegaM,h0)**3 - &
                   DCapprox((llist(j)/l21 - dz/2.d0-1),OmegaM,h0)**3)*1.d0/2.d0/pi**2
              !write(*,*) k/h0, llist(j),shellmu/Nmu/Nphi, IntensityNoise(llist(j)/l21-1,shellmu/Nmu/Nphi,secyears)/sqrt(Nkz)*h0**3
              !write k [1/Mpc], lambda, P_N [Mpc^3] mK^2
              if (shellphi < 0) shellphi = 0.d0
              InstrumentalNoise = IntensityNoise(llist(j)/l21-1,(shellphi/Nphi),secyears)/sqrt(Nkz)
              !make sure there are no infinities:
              if (InstrumentalNoise .ge. 1.E10) InstrumentalNoise = 1.E10
              write(21,"(E16.7,2X,E16.7,2X,F16.7,2X,E16.7)") llist(j), k, mulist(t),  &
                   InstrumentalNoise
              !write(*,"(E16.7,2X,E16.7,2X,F16.7,2X,E16.7, 2X, E16.7)") llist(j), k, mulist(t), shellphi/Nphi,  &
              !     InstrumentalNoise
              !shellmu/Nmu/Nphi
           enddo !muloop

        enddo !kloop
     enddo !lambdaloop
     close(21)
  else
     open(unit=21,file = 'Noise_CHIME_1D.dat', status='replace')

     do i = 1, Nk

        !here you can choose if you want log or lin sampling in k 
        k = klinlist(i)
        do j = 1, Nlambda
           shellmu = 0.d0
           !averaging over phi/mu:
           do t = 1, Nmu
              shellphi = 0.d0
              do l = 1, Nphi

                 uL = k*sqrt(1-mulist(t)**2)*Cos(philist(l))*llist(j)/(2.d0*pi)* DCapprox(llist(j)/L21-1.,Omegam,h0)
                 uW = k*sqrt(1-mulist(t)**2)*Sin(philist(l))*llist(j)/(2.d0*pi)* DCapprox(llist(j)/L21-1.,Omegam,h0)

                 !outside this range, we know the visability is 0. Interpolation will give out of bounds.
                 !there is also no zero mode
                 if (abs(uW) > abs(uW0)-1 .or. abs(uL) > abs(uL0)-1 .or. (abs(uL) .eq. 0.d0 .and. abs(uW) .eq. 0d0)) then
                    y = 0.d0
                 else                
                    call Simple2DLinInterp(uW,uL,uW0,uL0,uWstep,uLstep,VisFunction,y)
                 endif
                 !write(*,*) k, uL, uW, y
                 !add shell:
                 shellphi = shellphi + nbase*y/norm

              enddo
              if (shellphi < 0) shellphi = 0.d0
              shellmu = shellmu + shellphi
           enddo
           !compute the # of modes in this shell:
           Nkz = fsky*k**2*dk* &
                (DCapprox((llist(j)/l21 + dz/2.d0-1),OmegaM,h0)**3 - &
                DCapprox((llist(j)/l21 - dz/2.d0-1),OmegaM,h0)**3)*1.d0/2.d0/pi**2
           !write(*,*) k/h0, llist(j),shellmu/Nmu/Nphi, IntensityNoise(llist(j)/l21-1,shellmu/Nmu/Nphi,secyears)/sqrt(Nkz)*h0**3
           !write k [1/Mpc], lambda, P_N [Mpc^3] mK^2
           write(21,"(E16.7,2X,E16.7,2X,F16.7)") k, llist(j), &
                IntensityNoise(llist(j)/l21-1,(shellmu/Nmu/Nphi),secyears)/sqrt(Nkz)
           !shellmu/Nmu/Nphi
        enddo
     enddo
     close(21)
  endif


  !deallocate what is left:
  deallocate(VisFunction)
  deallocate(kloglist)
  deallocate(llist)
  deallocate(klinlist)
  deallocate(philist)
  deallocate(mulist)

contains
  !define wedge function that describes visability of cylindrical array:
  real(dl) function Lambda(x)
    real(dl) :: x

    if (abs(x) .lt. 1.d0) then
       Lambda = 1.d0-abs(x)
    else
       Lambda = 0.d0
    endif

  end function Lambda

  !simple linear 2D interpolation scheme:
  subroutine Simple2DLinInterp(x,y,xmin, ymin, dx,dy,f, fout)
    real(dl), intent(in)  :: x

    real(dl), intent(in)  :: y
    real(dl), intent(in)  :: dx, dy
    real(dl), intent(in)  :: f(:,:)
    real(dl), intent(out) :: fout
    real(dl), intent(in)  :: xmin, ymin
    integer :: s1p, s2p, s1m, s2m
    real(dl) :: slope1, slope2, slope3
    real(dl) :: f1, f2

    !determine the box, given x and y

    s1m=floor((x-xmin)/dx)+1
    s1p = s1m+2
    s2m=floor((y-ymin)/dy)+1   
    s2p = s2m+2

    !How does it work. Very simple:
    ! (s1m,s2p)       (s1p,s2p)      
    !
    !        ?.(x,y)
    !
    !
    !
    ! (s1m,s2m)       (s1p,s2m)

    !fac is to normalize visabilty function, i.e int du n(u) = nbase
    slope1 = (f(s1p,s2m)-f(s1m,s2m))/dx
    slope2 = (f(s1p,s2p)-f(s1m,s2p))/dx

    f1 = f(s1m,s2m) + slope1*(x-(xmin+s1m*dx))
    f2 = f(s1m,s2p) + slope2*(x-(xmin+s1m*dx))

    slope3 = (f2-f1)/dy

    fout = f1 + slope3*(y-(ymin+s2m*dy))

  end subroutine Simple2DLinInterp

  !comving radial distance 
  !below approximation works pretty well in our Universe   
  real(dl) function DCapprox(z,Omegam,h)
    real(dl), intent(in) :: z
    real(dl), intent(in) :: Omegam
    real(dl) :: a, b
    real(dl), intent(in)  :: h
    a = 1.718*Omegam
    b = 0.315*sqrt(Omegam)
    Dcapprox = 3000.*z/sqrt(1+a*z+b*z**2)/h !in Mpc

  end function DCapprox

  !Hubble (what need I say more)
  real(dl) function Hubble(z)
    real(dl), intent(in) :: z

    Hubble = 100*h0*sqrt(OmegaM*(z+1)**3 + 1 - OmegaM) !in km/Mpc/s
  end function Hubble

  !covert visability function to noise curve. 
  real(dl) function IntensityNoise(z,visab,int_time)
    real(dl), intent(in) :: visab
    real(dl), intent(in) :: z
    real(dl), intent(in) :: int_time
    real(dl) :: Tsys
    real(dl) :: Aeff
    real(dl) :: OmegaFoV
    real(dl) :: ra, yz

    Aeff = WidthCyl*LengthCyl !in m^2. Aeff is entire area of dish. Not per anntenna (Xin)
    ra = DCapprox(z,Omegam,h0) !in Mpc
    OmegaFOV = (l21*(z+1))**2*nantenna/Aeff 
    Tsys = 50000 !(*mK*) This I think is a fairly general assumption 
    yz = l21/Hubble(z)*(z+1)**2/1000. ! Mpc^-1 s
    !talked to Xin, not clear why there is no l^2/Aeff here 
    IntensityNoise =2*4.d0*pi*ra**2*yz*fsky*Tsys**2*(l21**2*(z+1)**2) &
         /Aeff/OmegaFOV/int_time/visab
    !in Mpc^3
  end function IntensityNoise

end program visability
