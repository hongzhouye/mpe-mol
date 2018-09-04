      SUBROUTINE partition(nocc,modims,auxdims,PQ,PQ_inv,auxmomo,
     \ oneeint,oneekin,alpha,pair)

!      use omp_lib

      Integer*4 nocc,modims,auxdims
      Real*8 PQ(auxdims,auxdims), PQ_inv(auxdims,auxdims)
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), alpha, pair(modims,nocc)
      Integer*4 i,j
      Real*8 lamb(auxdims)

!      write ( *, '(a,i8)' ) '  The number of processors available = ',
!     \                       omp_get_num_procs()
!      write ( *, '(a,i8)' ) '  The number of threads available    = ',
!     \                       omp_get_max_threads()

      pair=0d0
      Do i=1,nocc
          pair(i,i)=1d0
      End do

      lamb=0d0
      Call SCF(pair,lamb,alpha,nocc,modims,auxdims,PQ,PQ_inv,
     \         auxmomo,oneeint,oneekin)

      end subroutine partition


      SUBROUTINE SCF(pair,lamb,alpha,nocc,modims,auxdims,PQ,PQ_inv,
     \         auxmomo,oneeint,oneekin)

      Integer*4 nocc,modims,auxdims
      Real*8 PQ(auxdims,auxdims), PQ_inv(auxdims,auxdims)
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), alpha, pair(modims,nocc)

      Integer*4 i,j,k,P,Q,jter
      Real*8 lamb(auxdims),daux(auxdims,nocc),F(modims,modims,nocc)
      Real*8 tmp2d(modims,modims),tmp1d(auxdims),alpha2,beta
      Real*8 denmat(modims,modims,nocc),denmatold(modims,modims,nocc)
      Real*8 ER,Ekin,obj,objold,denmatdiff(modims,modims,nocc)
      Real*8 eigs(modims),eigv(modims,modims),diff,diffE,tmp1d2(modims)

      alpha2=1d0; beta=0d0
!     RI coefficient for each pair density
      daux=0d0
      do i=1,nocc
        tmp1d=0d0
        do Q=1,auxdims
!          tmp1d(Q)=Dot_product(pair(:,i),
!     \             MatMul(pair(:,i),auxmomo(:,:,Q)))
          Call DGemm("N","N",1,modims,modims,alpha2,
     \               pair(:,i),1,auxmomo(:,:,Q),modims,beta,tmp1d2,1)
          Call DGemm("N","N",1,1,modims,alpha2,
     \               pair(:,i),1,tmp1d2,modims,beta,tmp1d(Q),1)
        end do
!        daux(:,i)=MatMul(tmp1d,PQ_inv)
        Call DGemm("N","N",1,auxdims,auxdims,alpha2,
     \            tmp1d,1,PQ_inv,auxdims,beta,daux(:,i),1)
      end do

!     Construct initial Fock matrix
      F=0d0
      do i=1,nocc
        tmp2d=0d0
        do P=1,auxdims
          tmp2d=tmp2d-2d0*alpha*daux(P,i)*auxmomo(:,:,P)
        end do
        F(:,:,i)=oneeint+tmp2d
      end do

!     Solve chemical potential lamb to match total density
      Call densmatch(F,lamb,alpha,nocc,modims,auxdims,PQ,PQ_inv,
     \               auxmomo,oneeint,oneekin)

!     Compute density matrix
      denmat=0d0
      do i=1,nocc
          denmat(:,:,i)=spread(pair(:,i),dim=1,ncopies=modims)
     \                  *spread(pair(:,i),dim=2,ncopies=modims)
      end do
      denmatold=denmat

!     Compute Edmiston_Ruedenberg Localization
      ER=0d0
      do i=1,nocc
          ER=ER-4d0*Dot_product(daux(:,i),MatMul(PQ,daux(:,i)))
      end do

!     Compute kinetic energy
      Ekin=0d0
      do i=1,nocc
          Ekin=Ekin+Dot_product(pair(:,i),MatMul(oneekin,pair(:,i)))
      end do

!     Objective function
      obj=Ekin+alpha/4d0*ER
      objold=obj

!     Construct Fock matrix
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+lamb(P)*auxmomo(:,:,P)
      end do
      do i=1,nocc
          F(:,:,i)=F(:,:,i)+tmp2d
          Call Diagonalize(modims,F(:,:,i),eigs,eigv,k)
          pair(:,i)=eigv(:,1)
      end do
      
!     Compute density matrix
      denmat=0d0
      do i=1,nocc
          denmat(:,:,i)=spread(pair(:,i),dim=1,ncopies=modims)
     \                  *spread(pair(:,i),dim=2,ncopies=modims)
      end do

!     Compute Edmiston_Ruedenberg Localization
      ER=0d0
      do i=1,nocc
          ER=ER-4d0*Dot_product(daux(:,i),MatMul(PQ,daux(:,i)))
      end do

!     Compute kinetic energy
      Ekin=0d0
      do i=1,nocc
          Ekin=Ekin+Dot_product(pair(:,i),MatMul(oneekin,pair(:,i)))
      end do

!     Objective function
      obj=Ekin+alpha/4d0*ER

!     Compare the old and new pair densities
      diffE=abs(obj-objold)
      denmatdiff=denmat-denmatold
      diff=0d0
      do i=1,nocc; do j=1,modims; do k=1,modims
          diff=diff+denmatdiff(k,j,i)**2
      end do; end do; end do
      diff=Sqrt(diff)

!     Loop until density matrix difference is small
      jter=0
!      Do while (diff.gt.1d-8 .and. diffE.gt.1d-14 .and. jter.lt.1000)
      Do while (diff.gt.1d-8 .and. diffE.gt.1d-14 .and. jter.lt.2000)
          jter=jter+1
      !   RI coefficient for each pair density
          daux=0d0
          do i=1,nocc
            tmp1d=0d0
            do Q=1,auxdims
              Call DGemm("N","N",1,modims,modims,alpha2,
     \               pair(:,i),1,auxmomo(:,:,Q),modims,beta,tmp1d2,1)
              Call DGemm("N","N",1,1,modims,alpha2,
     \               pair(:,i),1,tmp1d2,modims,beta,tmp1d(Q),1)
            end do
            Call DGemm("N","N",1,auxdims,auxdims,alpha2,
     \            tmp1d,1,PQ_inv,auxdims,beta,daux(:,i),1)
          end do
  
      !   Construct Fock matrix
          F=0d0
          do i=1,nocc
            tmp2d=0d0
            do P=1,auxdims
              tmp2d=tmp2d-2d0*alpha*daux(P,i)*auxmomo(:,:,P)
            end do
            F(:,:,i)=oneeint+tmp2d
          end do

      !   Solve chemical potential lamb to match total density
          Call densmatch(F,lamb,alpha,nocc,modims,auxdims,PQ,PQ_inv,
     \                   auxmomo,oneeint,oneekin)

      !   Save old parameters
          objold=obj
          denmatold=denmat

      !   Construct Fock matrix
          tmp2d=0d0
          do P=1,auxdims
              tmp2d=tmp2d+lamb(P)*auxmomo(:,:,P)
          end do
          do i=1,nocc
              F(:,:,i)=F(:,:,i)+tmp2d
              Call Diagonalize(modims,F(:,:,i),eigs,eigv,k)
              pair(:,i)=eigv(:,1)
          end do

      !   Compute density matrix
          denmat=0d0
          do i=1,nocc
              denmat(:,:,i)=spread(pair(:,i),dim=1,ncopies=modims)
     \                      *spread(pair(:,i),dim=2,ncopies=modims)
          end do

      !   Compute Edmiston_Ruedenberg Localization
          ER=0d0
          do i=1,nocc
              ER=ER-4d0*Dot_product(daux(:,i),MatMul(PQ,daux(:,i)))
          end do

      !   Compute kinetic energy
          Ekin=0d0
          do i=1,nocc
              Ekin=Ekin+Dot_product(pair(:,i),MatMul(oneekin,pair(:,i)))
          end do

      !   Objective function
          obj=Ekin+alpha/4d0*ER

      !   Compare the old and new pair densities
          diffE=abs(obj-objold)
          denmatdiff=denmat-denmatold
          diff=0d0
          do i=1,nocc; do j=1,modims; do k=1,modims
              diff=diff+denmatdiff(k,j,i)**2
          end do; end do; end do
          diff=Sqrt(diff)

          Write(6,*)'jter,diffE,diff,obj,ER,Ekin',
     /               jter,diffE,diff,obj,ER,Ekin

      end do

      End subroutine SCF


      SUBROUTINE densmatch(Fold,lamb,alpha,nocc,modims,auxdims,PQ,
     \               PQ_inv,auxmomo,oneeint,oneekin)

!      use omp_lib
      
      Integer*4 nocc,modims,auxdims
      Real*8 PQ(auxdims,auxdims), PQ_inv(auxdims,auxdims)
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), alpha, pair(modims,nocc)

      Integer*4 i,j,k,P,Q,jter,a
      Real*8 lamb(auxdims),daux(auxdims,nocc),Fold(modims,modims,nocc)
      Real*8 tmp2d(modims,modims),F(modims,modims),eigs(modims),norm
      Real*8 eigv(modims,modims),eigs_store(modims,nocc),W(auxdims)
      Real*8 eigv_store(modims,modims,nocc),tmp2dri(modims,auxdims)
      Real*8 Resid,deriv(auxdims,auxdims),der1(auxdims),dlamb(auxdims)
      Real*8 eps,t1,t2,t3,t4, der2(modims,auxdims),beta,alpha2
      Real*8 tmp1d2(modims)

      alpha2=1d0; beta=0d0
      pair=0d0; eigv_store=0d0; eigs_store=0d0
!     Construct full Fock matrix
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+lamb(P)*auxmomo(:,:,P)
      end do
      do i=1,nocc
          F=Fold(:,:,i)+tmp2d
          Call Diagonalize(modims,F,eigs,eigv,k)
          pair(:,i)=eigv(:,1)
          eigv_store(:,:,i)=eigv; eigs_store(:,i)=eigs
      end do

!     Compute the dnesity difference
      W=0d0
      do i=1,nocc
          do P=1,auxdims
!              W(P)=W(P)+2d0*Dot_product(pair(:,i),
!     \                MatMul(auxmomo(:,:,P),pair(:,i)))
              Call DGemm("N","N",modims,1,modims,alpha2,auxmomo(:,:,P),
     \                   modims,pair(:,i),modims,beta,tmp1d2,modims)
              Call DGemm("N","N",1,1,modims,2d0,
     \               pair(:,i),1,tmp1d2,modims,1d0,W(P),1)
          end do
          W=W-2d0*auxmomo(i,i,:)
      end do
      Resid=0d0
      Do i=1,auxdims; Resid=Resid+W(i)**2; End do; Resid=Sqrt(Resid)

      jter=0; norm=1d0
!      Do while (norm.gt.1d-5 .and. jter.lt.1000)
      Do while (norm.gt.1d-5 .and. jter.lt.2000)
          jter=jter+1
      !   Compute the response kernel dW/dlamb
          deriv=0d0
!          Call cpu_time(t1)
!          do i=1,nocc; do a=2,modims
!              do P=1,auxdims
!                  der1(P)=Dot_product(eigv_store(:,a,i),
!     \                    MatMul(auxmomo(:,:,P),pair(:,i)))
!                tmp1d2=0d0
!                Call DGemm("N","N",modims,1,modims,alpha2,
!     \                     auxmomo(:,:,P),modims,pair(:,i),
!     \                     modims,beta,tmp1d2,modims)
!                Call DGemm("N","N",1,1,modims,alpha2,
!     \                     eigv_store(:,a,:i),1,tmp1d2,
!     \                     modims,beta,der1(P),1)
!              end do
!              do Q=1,auxdims; do P=1,auxdims
!                  deriv(P,Q)=deriv(P,Q)+4d0*der1(P)*der1(Q)
!     $                        /(eigs_store(1,i)-eigs_store(a,i))
!              end do; end do
!          end do; end do
!          Call cpu_time(t2)
!          print '("Time = ",f9.3," seconds.")',t2-t1

!          wtime = omp_get_wtime()
!          Call cpu_time(t1)
!          alpha2=1d0; beta=0d0
!$OMP DO
          do i=1,nocc
              do P=1,auxdims
!                  der2(:,P)=MatMul(pair(:,i),
!     \                      MatMul(auxmomo(:,:,P),eigv_store(:,:,i)))
                tmp2d=0d0
                Call DGemm("N","N",modims,modims,modims,alpha2,
     \                     auxmomo(:,:,P),modims,eigv_store(:,:,i),
     \                     modims,beta,tmp2d,modims)
                Call DGemm("N","N",1,modims,modims,alpha2,
     \                     pair(:,i),1,tmp2d,modims,beta,der2(:,P),1)
              end do
              do a=2,modims
                do Q=1,auxdims; do P=1,auxdims
                  deriv(P,Q)=deriv(P,Q)+4d0*der2(a,P)*der2(a,Q)
     $                        /(eigs_store(1,i)-eigs_store(a,i))
                end do; end do
              end do
          end do
!$OMP END DO
!          Call cpu_time(t2)
!          print '("Time = ",f9.3," seconds.")',t2-t1
!          wtime = omp_get_wtime() - wtime
!          write(6,*)'wtime',wtime
!          write(6,*)'deriv',deriv(1,1)

      !   Get step in lamb
          if (Resid.gt.1d-4) then
              Call GetLSQ(auxdims,auxdims,W,deriv,0d0,dlamb)
          else
              Call GetLSQ(auxdims,auxdims,W,deriv,0d0,dlamb)
          end if
          norm=0d0
          Do i=1,auxdims; norm=norm+dlamb(i)**2; End do; norm=Sqrt(norm)

      !   Trust radius
          eps=5d-1
          if (norm.gt.eps) then
              dlamb=dlamb/norm*eps
          end if

      !   Update the potential and Fock matrix
          lamb=lamb+dlamb
          tmp2d=0d0
          do P=1,auxdims
              tmp2d=tmp2d+lamb(P)*auxmomo(:,:,P)
          end do
          do i=1,nocc
              F=Fold(:,:,i)+tmp2d
              Call Diagonalize(modims,F,eigs,eigv,k)
              pair(:,i)=eigv(:,1)
              eigv_store(:,:,i)=eigv; eigs_store(:,i)=eigs
          end do
      !   Update the density difference
          W=0d0
          do i=1,nocc
              do P=1,auxdims
!                  W(P)=W(P)+2d0*Dot_product(pair(:,i),
!     \                    MatMul(auxmomo(:,:,P),pair(:,i)))
                Call DGemm("N","N",modims,1,modims,alpha2,
     \                    auxmomo(:,:,P),
     \                   modims,pair(:,i),modims,beta,tmp1d2,modims)
                Call DGemm("N","N",1,1,modims,2d0,
     \               pair(:,i),1,tmp1d2,modims,1d0,W(P),1)
              end do
              W=W-2d0*auxmomo(i,i,:)
          end do
          Resid=0d0
          Do i=1,auxdims; Resid=Resid+W(i)**2; End do; Resid=Sqrt(Resid)

          if ((jter/100)*100.eq.jter) then
              write(6,*)'#jter,Resid,norm',jter,Resid,norm
          end if
      end do

      write(6,*)'#jter,Resid,norm',jter,Resid,norm
          
      End subroutine densmatch

      
