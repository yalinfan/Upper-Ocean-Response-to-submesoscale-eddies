      subroutine get_ic_bv(ind,n,m,l,r,e,u,v,zw,elon,alat)
c---- subroutine to set up the initial field for simulating a baroclinic vortex based on Penven et al. (2006) 
c---- and including a Doppler shift to propagate the eddy
c
c       ind = flag to indicate if perturbation for internal wave is
c             to be included:  =0 no; =1 yes.
c       n,m,l = grid dimensions.  l = number of layers + 1.
c       r     = T and S fields.
c
      implicit none
c
c  declare passed variables.
      integer ind,n,m,l
      real*4  r(n,m,l-1,2),e(n,m),u(n,m,l-1),v(n,m,l-1)
      real*4  elon(n,m),alat(n,m),zw(l)
c
c  declare local temporary variables.
c
      integer i,j,k
      real*4 Nfreq,rho0,f0,lam,umax,h1,P0,Ome2,g
      real*4 a_f,fij,lam2,Fz,dFz,rhoij,x,y,eddy_sp
c
      g     = 9.81
      Ome2  = 2*7.292e-5                 !2*Omega rad/sec
      Nfreq = 0.003                       !Brunt-Vaissala freq (1/sec), .003
      rho0  = 1025.0                      !mean density kg/m^3
      f0    = Ome2*sin(10.*3.1415/180.0) !f at 10.N
      lam   = 4.0e3/111.2e3              !e-folding scale in deg
      umax  = 1.                         !max surface geostrophic velocity m/sec
      h1    = 300                        !level of no-motion meter
      eddy_sp = 0.0    !background advection of eddy
c
      lam2 = lam**2
      P0 = rho0*f0*umax*(lam*111.2e3)*sqrt(exp(1.0))
c
      if (ind .eq. 1) then
        do j=1,m
          do i=1,n
            x=elon(i,j)-(150.4-(300/111.2)) ! X km from right boundary
            y=alat(i,j)-10.  ! centered in model domain
            a_f=exp( -(x**2+y**2)/(2.0*lam2) )
            fij=Ome2*sin(alat(i,j)*3.1415/180.0)
c Penven et al (2006)
c            e(i,j)=P0*a/(rho0*g-P0*a*(1.0-exp(-h1))/(h1-1.0+exp(-h1)))
c
c Yu et al. (2020)
            e(i,j)=P0*a_f/( rho0*g-P0*a_f*(exp(1.0)-1.0)/h1 )
            do k=1,l-1
              if ( zw(k) > -h1 ) then
c Penven et al(2006)
c               Fz=P0*( h1-1.0+zw(k)+exp( -(zw(k)+h1) ) )
c    &             /( h1-1.0+exp(-h1) )
c               dFz=P0*( 1.0-exp(-(zw(k)+h1)) )/( h1-1.0+exp(-h1) )
c 
c Yu et al. (2020) 
                Fz =P0*(zw(k)/h1+exp( -(1.0+zw(k)/h1) ))*exp(1.0)
                dFz=P0*exp(1.0)/h1*(1.0-exp( -(1.0+zw(k)/h1)))
                rhoij = rho0*( 1.0-Nfreq**2*zw(k)/g )-dFz*a_f/g
                u(i,j,k)=-1.0/fij/rhoij*( -y/lam2/111.2e3 )*Fz*a_f
                u(i,j,k)= eddy_sp+u(i,j,k) ! add Doppler shift
                v(i,j,k)= 1.0/fij/rhoij*( -x/lam2/111.2e3 )*Fz*a_f
                r(i,j,k,1)=(1030-rhoij)/0.28
                r(i,j,k,2)=33.0-0.008*zw(k)
                if ( r(i,j,k,2) > 35.0) then
		    r(i,j,k,2)=35.0
                endif
              else
                rhoij = rho0*( 1.0-Nfreq**2*(-h1)/g )
                u(i,j,k)=0.0
                u(i,j,k)= eddy_sp+u(i,j,k) ! add Doppler shift
                v(i,j,k)=0.0
                r(i,j,k,1)= (1030-rhoij)/0.28
                r(i,j,k,2)=33.0-0.008*zw(k)
                if ( r(i,j,k,2) > 35.0) then
		    r(i,j,k,2)=35.0
                endif
              endif
            enddo
          enddo
        enddo
      endif

      if (ind .eq. 0) then  !get background T, S
        do j=1,m
          do i=1,n
          e(i,j)=0.0
            do k=1,l-1
              rhoij = rho0*( 1.0-Nfreq**2*zw(k)/g )
              r(i,j,k,1)= (1030-rhoij)/0.28
              r(i,j,k,2)=33-0.008*zw(k)
              if ( r(i,j,k,2) > 35.0) then
		  r(i,j,k,2)=35.0
              endif
              u(i,j,k)=0.0
              u(i,j,k)= eddy_sp+u(i,j,k) !add Doppler shift
              v(i,j,k)=0.0
            enddo
          enddo
        enddo
      endif

      return
      end

