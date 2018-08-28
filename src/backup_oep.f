      SUBROUTINE oep(coeff,vx,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,VV,DFTmo,Energy)

      Integer*4 nocc,modims,auxdims,natom
      Real*8 coeff(modims,nocc), vx(auxdims),Ek,Energy
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims) 
      Integer*4 i,j
      Real*8 nuc_cons,alpha,beta,tmp1d(modims),W
      Real*8 VV(modims,modims,modims,modims),DFTmo(modims,modims)

      alpha=1d0; beta=0d0
!     Compute the external potential energy (constant)
      nuc_cons=0d0
      do i=1,nocc
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               coeff(:,i),1,oneenuc,modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,-2d0,
     \               coeff(:,i),1,tmp1d,modims,1d0,nuc_cons,1)
      end do

!     Newton's Optimization to search for OEP
      Call newton(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,1d-12,1000)

!     Final Lagrangian and kinetic energy
      Call func(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,W)

      Call kinetic(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,Ek)

!     W-T is 0 for exact density matching
      Write(6,*)'#final W T',W,Ek

!     Compute DFT energy using grids
      open(unit=21,File='grid_dimension.txt')
      read(21,*)natom
      close(21)

      Call DFTEnergy(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,VV,DFTmo,natom,Energy)

      End subroutine oep

!     Compute the DFT energy
      Subroutine DFTEnergy(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,VV,DFTmo,natom,Energy)

      Integer*4 nocc,modims,auxdims,i,j,k,P,mm,nn
      Real*8 coeff(modims,nocc), vx(auxdims),nuc_cons,alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),Energy
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims),W
      Real*8 VV(modims,modims,modims,modims),rhosum,DFTmo(modims,modims)
      Real*8 aoeigv(modims,modims),drho(3),drhom
      Real*8,allocatable ::  Grid(:,:), aogrid(:,:), densgrid(:,:)
      Integer*4 stat, Npoints(natom),Nptot,Npsum, Nbasis, l,natom
      Real*8 rho,ex,ec,rs,X,Xx,Xx0,Q,A,b,c,exB88
      Real*8 ExB88sum,ExLDAsum,EcLDAsum,ExEXX

      Real*8, parameter :: PI = 4.d0*atan(1.d0)

!     Read numbers of grid points, atomic basis functions and atoms
      open(unit=21,File='grid_dimension.txt')
      read(21,*)natom,Npoints,Nbasis
      close(21)

      Nptot=0
      do i=1,natom
        Nptot=Nptot+Npoints(i)
      end do

!     Read grid coordinates and weights
      open(unit=21,File='grid.txt')
      allocate(Grid(Nptot,4))
      do i=1,Nptot
        read(21,*) Grid(i,:)
      end do
      close(21)

!     Read AO basis values on grids
      open(unit=21,File='ao_grid.txt')
      allocate(aogrid(Nptot*Nbasis,4))
      do i=1,Nptot*Nbasis
        read(21,*) aogrid(i,:)
      end do
      close(21)

!     Construct Fock matrix and diagonalize
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,k)

!     Compute EXX exchange energy
      F=0d0
      Do i=1,nocc
        Do j=1,modims; Do k=1,modims
          F=F+eigv(j,i)*eigv(k,i)*(-VV(:,j,k,:))
        end do; end do
      end do

      ExEXX=0d0
      Do i=1,nocc
        ExEXX=ExEXX+Dot_product(eigv(:,i),MatMul(F,eigv(:,i)))
      end do
      write(6,*)'ExEXX',ExEXX

!     Compute J
      F=0d0
      Do i=1,nocc
        Do j=1,modims; Do k=1,modims
          F=F+eigv(j,i)*eigv(k,i)*(2d0*VV(:,:,k,j)-0*VV(:,j,k,:))
        end do; end do
      end do

      Energy=0d0
      Do i=1,nocc
        Energy=Energy+Dot_product(eigv(:,i),
     \                   MatMul(F,eigv(:,i)))
      end do
      write(6,*)'J',Energy

!     Compute T+V
      Do i=1,nocc
        Energy=Energy+Dot_product(eigv(:,i),
     \                   MatMul(2d0*oneeint,eigv(:,i)))
      end do
      write(6,*)'T+V+J',Energy

!     Transfer orbitals to AO basis
      aoeigv=MatMul(DFTmo,eigv)

      ExB88sum=0d0; ExLDAsum=0d0; EcLDAsum=0d0
!     Compute spin (alpha or beta) electron density and density gradient on grids
      do nn=1,natom
        Npsum=0d0
        do i=1,nn-1
          Npsum=Npsum+Npoints(i)
        end do
        Do l=1,Npoints(nn)
          rho=0.d0
          do k=1,nocc
            do i=1,modims; do j=1,modims
              rho = rho + aoeigv(i,k)*aoeigv(j,k)
     \            *aogrid(Npsum*Nbasis+Npoints(nn)*(i-1)+l,1)
     \            *aogrid(Npsum*Nbasis+Npoints(nn)*(j-1)+l,1)
            end do; end do
          end do

          drho=0d0
          do mm=1,3
            do k=1,nocc
              do i=1,modims; do j=1,modims
                drho(mm) = drho(mm) + aoeigv(i,k)*aoeigv(j,k)
     \          *(aogrid(Npsum*Nbasis+Npoints(nn)*(i-1)+l,mm+1)
     \          *aogrid(Npsum*Nbasis+Npoints(nn)*(j-1)+l,1)
     \          +aogrid(Npsum*Nbasis+Npoints(nn)*(i-1)+l,1)
     \          *aogrid(Npsum*Nbasis+Npoints(nn)*(j-1)+l,mm+1))
              end do; end do
            end do
          end do

        !  Slater exchange energy density
          ex = -(3./2.)*(3./PI/4.)**(1./3.)*(rho)**(4./3.)*2

        ! B88 exchange energy density
          beta=0.0042
          drhom=sqrt(drho(1)**2+drho(2)**2+drho(3)**2)
          x=drhom/rho**(4./3.)
          if(rho.gt.1d-12.and.drhom.gt.1d-12)then
            exB88=ex-2.*beta*rho**(4./3.)*x**2/(1+6.*beta*x*asinh(x))
          else
            exB88=0.
          end if

        !LSDA correlation energy density (VWN5 for spin-compensated only)
          A = 0.0310907
          p = -0.10498
          c = 3.72744
          d = 12.9352
          Q = (4*d-c**2)**(1./2.)
          x = 0.25*(3**(1./6.))*(4**(5./6.))*(PI*rho*2.)**(-1./6.)
          Xx = x**2 + c*x + d
          Xp = p**2 + c*p + d

          if(rho.gt.1d-12)then
            ec = A*(log(x**2/Xx)+2*c*atan(Q/(2*x+c))*(Q**(-1))
     &       -c*p*(log((x-p)**2/Xx)+2*(c+2*p)*atan(Q/(2*x+c))*(Q**(-1)))
     &        *Xp**(-1))
          else
            ec=0.
          end if
          ec = 2.*rho*ec

        ! Compute integrated DFT energy
          Energy = Energy + Grid(Npsum+l,4)*(ex+ec)

          ExB88sum = ExB88sum + Grid(Npsum+l,4)*(exB88)
          ExLDAsum = ExLDAsum + Grid(Npsum+l,4)*(ex)
          EcLDAsum = EcLDAsum + Grid(Npsum+l,4)*(ec)

          rhosum=rhosum+Grid(Npsum+l,4)*rho*2

        end do
      end do
      write(6,*)'rhosum',rhosum
      write(6,*)'ExB88,ExLDA',ExB88sum,ExLDAsum
      write(6,*)'EcLDA',EcLDAsum

      deallocate(Grid)
      deallocate(aogrid)

      End Subroutine DFTEnergy

!     Compute the Lagragian
      Subroutine func(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,W)
      
      Integer*4 nocc,modims,auxdims,i,j,k,P
      Real*8 coeff(modims,nocc), vx(auxdims),nuc_cons,alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)

      alpha=1d0; beta=0d0
!     Construct Fock matrix and diagonalize
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,k)

      W=nuc_cons
      do i=1,nocc
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,oneeint,modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,2d0,
     \               eigv(:,i),1,tmp1d,modims,1d0,W,1)
      end do

      do i=1,nocc
        tmp1d2=0d0
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,alpha,
     \               eigv(:,i),1,tmp1d,modims,beta,tmp1d2(P),1)
        end do
        Call DGemm("N","N",1,1,auxdims,2d0,vx,1,tmp1d2,auxdims,1d0,W,1)
      end do

      do i=1,nocc
        tmp1d2=0d0
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               coeff(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,alpha,
     \               coeff(:,i),1,tmp1d,modims,beta,tmp1d2(P),1)
        end do
        Call DGemm("N","N",1,1,auxdims,-2d0,vx,1,tmp1d2,auxdims,1d0,W,1)
      end do

      W=-W

      End subroutine func


!     Compute the gradient of Lagrangian (density matching)
      Subroutine gradient(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ eigs,eigv,deriv)

      Integer*4 nocc,modims,auxdims,i,j,k,P
      Real*8 coeff(modims,nocc), vx(auxdims),alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)
      Real*8 deriv(auxdims)

      alpha=1d0; beta=0d0
      deriv=0d0
      do i=1,nocc
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,2d0,
     \               eigv(:,i),1,tmp1d,modims,1d0,deriv(P),1)
        end do
      end do

      do i=1,nocc
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               coeff(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,-2d0,
     \               coeff(:,i),1,tmp1d,modims,1d0,deriv(P),1)
        end do
      end do

      deriv=-deriv

      End subroutine gradient


!     Compute the hessian of Lagrangian (needed for Newton's method)
      Subroutine hessian(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ eigs,eigv,hess)

      Integer*4 nocc,modims,auxdims,i,j,k,P,Q,a
      Real*8 coeff(modims,nocc), vx(auxdims),alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)
      Real*8 hess(auxdims,auxdims),der2(modims,auxdims)
      Real*8 der1(modims,modims,auxdims)

      alpha=1d0; beta=0d0
      hess=0d0
!      do i=1,nocc
!        do P=1,auxdims
!          tmp2d=0d0
!          Call DGemm("N","N",modims,modims,modims,alpha,
!     \             auxmomo(:,:,P),modims,eigv,modims,beta,tmp2d,modims)
!          Call DGemm("N","N",1,modims,modims,alpha,
!     \             eigv(:,i),1,tmp2d,modims,beta,der2(:,P),1)
!        end do
!        do a=nocc+1,modims
!          do Q=1,auxdims; do P=1,auxdims
!            hess(P,Q)=hess(P,Q)+4d0*der2(a,P)*der2(a,Q)
!     $                    /(eigs(i)-eigs(a))
!          end do; end do
!        end do
!      end do

      do P=1,auxdims
        tmp2d=0d0
        Call DGemm("N","N",modims,modims,modims,alpha,
     \             auxmomo(:,:,P),modims,eigv,modims,beta,tmp2d,modims)
        Call DGemm("N","N",modims,modims,modims,alpha,transpose(eigv),
     \             modims,tmp2d,modims,beta,der1(:,:,P),modims)
      end do
      do i=1,nocc
        do a=nocc+1,modims
          do Q=1,auxdims; do P=1,auxdims
            hess(P,Q)=hess(P,Q)+4d0*der1(i,a,P)*der1(i,a,Q)
     $                    /(eigs(i)-eigs(a))
          end do; end do
        end do
      end do

      hess=-hess

      End subroutine hessian


!     Newton's method
      Subroutine newton(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,gtol,maxiter)

      Integer*4 nocc,modims,auxdims,i,j,k,P,maxiter
      Real*8 coeff(modims,nocc), vx(auxdims),alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)
      Real*8 hess(auxdims,auxdims),deriv(auxdims),gtol
      Real*8 norm_grad,dvx(auxdims),step,eps

!     Construct Fock matrix and diagonalize
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,k)

      Call gradient(vx,coeff,nocc,modims,auxdims,auxmomo,
     \              eigs,eigv,deriv)
      norm_grad=0d0
      do i=1,auxdims; norm_grad=norm_grad+deriv(i)**2; end do
      norm_grad=sqrt(norm_grad)

      jter=0
      do while (norm_grad.gt.gtol .and. jter.lt.maxiter)
          jter=jter+1

          Call hessian(vx,coeff,nocc,modims,auxdims,auxmomo,
     \                 eigs,eigv,hess)

!         Get the search direction
          if (norm_grad.gt.1d-3) then
              Call GetLSQ(auxdims,auxdims,deriv,hess,1d-6,dvx)
          else if (norm_grad.gt.1d-5.and.norm_grad.lt.1d-3) then
              Call GetLSQ(auxdims,auxdims,deriv,hess,1d-12,dvx)
          else
              Call GetLSQ(auxdims,auxdims,deriv,hess,1d-18,dvx)
          end if

!         Trust radius
          step=0d0
          do i=1,auxdims; step=step+dvx(i)**2; end do; step=sqrt(step)
          eps=0.5d0
          if(norm_grad.lt.1d-5)eps=2d0
          if (step.gt.eps) then
              dvx=dvx/step*eps
          end if

!         Update potential
          vx=vx+dvx

!         Update Fock matrix
          tmp2d=0d0
          do P=1,auxdims
              tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
          end do
          F=oneeint+tmp2d
          Call Diagonalize(modims,F,eigs,eigv,k)

          Call gradient(vx,coeff,nocc,modims,auxdims,auxmomo,
     \                  eigs,eigv,deriv)
          norm_grad=0d0
          do i=1,auxdims; norm_grad=norm_grad+deriv(i)**2; end do
          norm_grad=sqrt(norm_grad)

          if ((jter/100)*100.eq.jter) then
              write(6,*)'#jter,norm_grad,step',jter,norm_grad,step
          end if
      end do

      write(6,*)'#jter,norm_grad,step',jter,norm_grad,step

      End subroutine newton


!     Compute kinetic energy
      Subroutine kinetic(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,Ek)

      Integer*4 nocc,modims,auxdims,i,j,k,P
      Real*8 coeff(modims,nocc), vx(auxdims),nuc_cons,alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),Ek
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)

      alpha=1d0; beta=0d0
!     Construct Fock matrix and diagonalize
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,k)

      Ek=0d0
      do i=1,nocc
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,oneekin,modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,2d0,
     \               eigv(:,i),1,tmp1d,modims,1d0,Ek,1)
      end do

      Ek=-Ek

      End Subroutine kinetic

