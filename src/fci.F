      Subroutine 
     $     FCIman(N,No,Nv,M,N2,Nstr,Max1,Max2,N0,h,htest,V,wi,X,P,Ecorr)
      Use FCIUtils
!     Solves for FCI ground state and returns the average values of the M
!     operators wi(:,;,i) in P(i)
C      Use Operations
      Implicit None
      Integer*4 N,No,Nv,N2,Nstr,Max1,Max2,N0,M
      Integer*4 i,j,k,l,mm,a,b,ij,kl,ab,ii,jj,ierr,iter,Info
      Integer*4 iSpin,iX
      Integer*4, Save :: ItTot
      Integer*4 Istr(N0)
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 Iocc(No),Isubst(No),Inda,Indb,Isigna,Isignb
      Real*8 Ecorr,q,htest(N,N),Ecorr2
      Real*8 Norm,DE,EX,EY,NormX,NormY,fac,thresh,DeltaE,scale,Tmp
      Real*8 X(Nstr,Nstr),X1(Nstr,Nstr),XH(Nstr,Nstr),X1H(Nstr,Nstr)
      Real*8 H0(N0*N0,N0*N0),E0(N0*N0),U0(N0*N0,N0*N0)
      Real*8 A0(N0*N0),B0(N0*N0),C0(N0*N0)
      Real*8 Hd(Nstr,Nstr)
      Real*8 Hm(2,2),eig(2),U(2,2),Ptest2(N,N)
      Real*8 h(N,N),V(N,N,N,N),wi(N,N,M),P(M),T1(N,N),Ttest2(N,N)
      Common /ZIndex/ Zptr
      Integer*4,Pointer :: Zptr(:,:)
      Integer*4,Target :: Zindex(N,N)
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

      Ecorr=zero; Zptr=>Zindex; Fac=One;
! Build indexing array for future use
      ZIndex=0
      Do k=1,No; Do l=k,N-No+k 
         If(k==No) Then
            ZIndex(k,l)=l-k
         Else
            Do mm=N-l+1,N-k
               ZIndex(k,l)=ZIndex(k,l)+NcR(mm,No-k)-NcR(mm-1,No-k-1)
            End Do
         End If
      End Do; End Do

!!----Determine which strings are connected by various operators.----!!
      Call IString(N,No,N2,Max1,Max2,Ex1,Ex2)

! Build Diagonal part of H
      Hd=Zero
      Call GetHd(N,No,Nstr,h,V,Iocc,Isubst,Hd,1)

! Get N0 lowest Strings
      Call GetIstr(N,No,Nstr,N0,Hd,Istr)

! Build + Diagonalize H0
      Call GetH0(N,No,N0,N2,Max1,Max2,Nstr,Ex1,Ex2,Istr,h,V,H0)
      Call Diagonalize(N0*N0,H0,E0,U0,Info)

! Initial Vector (ensure it is singlet)
      iSpin=-1; iX=0; Norm=MDot(X,X)
      If(Norm.lt.0.01d0) Then
         Do While(iSpin==-1)
            iX=iX+1
            X=Zero; ij=0
            Do i=1,N0; Do j=1,N0; ij=ij+1
               X(Istr(j),Istr(i))=U0(ij,iX)
            End Do; End Do
            X1=X-Transpose(X)
            Norm=MDot(X1,X1)/MDot(X,X)
            If(Norm<1d-2) iSpin=1
         End Do
      End If

      X=X/Sqrt(MDot(X,X))
      Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,XH)

      Ecorr=MDot(X,XH); DE=Ecorr/3; iter=0
      Thresh=1.0d-7
      Fac=1.0d0
      Do While(Abs(fac)>Thresh.and.iter<5000)
         iter=iter+1; DE=Ecorr
! Make the (orthogonal component of the) Davidson Update
         X1=-(XH-Ecorr*X)/(Hd-Ecorr)
! Olsen Update
!         DeltaE=MDot(X,X1)/MDot(X,X/(Hd-Ecorr))
!         X1=X1-DeltaE*X/(Hd-Ecorr)

! Efficiently Build the Davidson Update using H0
         ij=0; A0=0d0; B0=0d0; C0=0d0; X1H=XH-Ecorr*X
         Do i=1,N0; Do j=1,N0; ij=ij+1; 
            A0(ij)=X1H(Istr(j),Istr(i))
         End Do; End Do
         ij=0
         Do ij=1,N0*N0
            If(Abs(E0(ij)-Ecorr).lt.1d-4) Cycle
            Do kl=1,N0*N0
               B0(ij)=B0(ij)-A0(kl)*U0(kl,ij)
            End Do
            DeltaE=1d0/(E0(ij)-Ecorr)
            B0(ij)=B0(ij)*DeltaE
         End Do
         ij=0
         Do ij=1,N0*N0
            Do kl=1,N0*N0
               C0(ij)= C0(ij)+U0(ij,kl)*B0(kl)
         End Do; End Do
         ij=0
         Do i=1,N0; Do j=1,N0; ij=ij+1; 
            X1(Istr(j),Istr(i))=C0(ij)
         End Do; End Do

         X1=X1-MDot(X1,X)*X
         X1=(X1+Transpose(X1))/Two
         X1=X1/Sqrt(MDot(X1,X1))

         X1=X1-MDot(X1,X)*X
         X1=(X1+Transpose(X1))/Two
         X1=X1/Sqrt(MDot(X1,X1))

         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X1,h,V,X1H)

!Build Hm
         Hm(1,1)=MDot(X,XH)
         Hm(2,2)=MDot(X1,X1H)
         Hm(1,2)=MDot(X,X1H); Hm(2,1)=Hm(1,2)

! Diagonalize Hm
         Call GetEig2(Hm,Eig,U)

!  Keep Lowest Eigenvector
         If(iter<50.or.(iter/10)*10/=iter) Then
            fac=U(2,1)
            X=U(1,1)*X+U(2,1)*X1
            XH=U(1,1)*XH+U(2,1)*X1H
         Else
            fac=.50d0*U(2,1); Norm=One/Pythag(One,fac)
            X=Norm*(X+fac*X1)
            XH=Norm*(XH+fac*X1H)
         End If

         X=X/Sqrt(MDot(X,X)); XH=XH/Sqrt(MDot(X,X))

! Avoid accumulating roundoff error
         If((iter/4)*4==iter)
     $        Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,XH)

         Ecorr=MDot(X,XH)/MDot(X,X); DE=DE-Ecorr
!     If(iter==1)Write(6,*)'##Step    Energy             Predicted     
!     $        Delta_E          Residual'
!         If((iter/20)*20==iter-1)
!     $        Write(6,114)iter,Ecorr,Eig(1),DE,fac

 114     Format('##',I4,4f18.13)
 113     Format(I4,3f18.13)
      End Do

!     Compute Density
      Norm=MDot(X,X)
      Do i=1,M
         T1=0.5d0*wi(:,:,i)
         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,T1,V*0d0,XH)
         P(i)=MDot(X,XH)/Norm
      End Do

!      Do i=1,N; Do j=1,N
!         Ttest2=0d0; Ttest2(i,j)=0.5d0
!         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,Ttest2,V*0d0,XH)
!         Ptest2(i,j)=MDot(X,XH)/Norm
!      End Do; End do

!      P=0d0
!      Do i=1,N; Do j=1,N
!        P=P+Ptest2(j,i)*wi(i,j,:)
!      End do; End do

!      If(No<3) Then
!         T1=h; Do i=1,N; T1(i,i)=0d0; End Do
!         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,T1,V*0d0,XH)
!         Fac=MDot(X,XH)/Norm
!         Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h*0d0,V,XH)
!         Write(6,*)'#@@# FCI PE',MDot(X,XH)/Norm,Fac
!      End If

      ItTot=ItTot+iter

!     Get the true FCI energy without added potential
      Call HX(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,htest,V,XH)
      Ecorr=MDot(X,XH)/MDot(X,X)
!      Write(6,*)'####FCI Iterations',ItTot,Ecorr

      Return
      End Subroutine

      Subroutine Hx(N,No,N2,Max1,Max2,Nstr,Ex1,Ex2,X,h,V,Y)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N2,Max1,Max2,Nstr
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 i,j,k,l,ij,kl,ii,jj,I1,I2,iiMax,J1,J2,jjMax,ik,jl
      Integer*4 AEx1(Max1,N,N),IfZero(N,N),IfSym
      Real*8 X(Nstr,Nstr),h(N,N),V(N,N,N,N),Y(Nstr,Nstr)
      Real*8 Vtmp,VS,VSS,htmp,hS,Tmp,Spin
      Real*8 Xtmp(Max1,Nstr),Ytmp(Max1,NStr)
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

      Y=Zero; AEx1=Abs(Ex1)

!     Check Spin Symmetry of X
      Spin=One; Y=X-Transpose(X); Tmp=MDot(Y,Y)
      If(Tmp>0.1d0) Spin=-One
      Y=Zero;

! Check for Zero Blocks in V
      IfZero=0
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If(Abs(V(i,j,k,l))>1.0d-10)Go To 20
         End Do; End Do
         IfZero(i,k)=1
 20      Continue
      End Do; End Do

! Check Symmetry of V
      IfSym=1
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If (Abs(V(i,j,k,l)-V(j,i,l,k))>1.0d-10) Then
               IfSym=0
               Go To 40
            End If
         End Do; End Do
      End Do; End Do
 40   Continue

! One Electron Part
      Do i=1,N
         Do j=1,N
            htmp=h(i,j)
            If(Abs(htmp)<1.0d-10)Cycle
            If(i==j) Then
               iiMax=NcR(N-1,No-1)
            Else
               iiMax=NcR(N-2,No-1)
            End If
            Do ii=1,iiMax
               I1=AEx1(ii,i,j)
               I2=AEx1(ii,j,i)
               hS=htmp
               If(I1/=Ex1(ii,i,j))hS=-hS
               Call DaxPy(Nstr,hS,X(1,I1),1,Y(1,I2),1)
            End Do
         End Do
      End Do
 143  Format(4I4,3f14.8)

! Same Spin Two-Electron Part
      ij=0
      Do i=1,N; Do j=i+1,N; ij=ij+1; kl=0
         Do k=1,N; Do l=k+1,N; kl=kl+1
            Vtmp=(V(i,j,k,l)-V(i,j,l,k)-V(j,i,k,l)+V(j,i,l,k))/Two
            If(abs(Vtmp)<1.0d-10)Go To 100 
            If(i==k.and.j==l) Then
               iiMax=NcR(N-2,No-2)
            Else If(i==k) Then  !i=2,3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(i==l) Then  !i=1-3
               iiMax=NcR(N-3,No-2)
            Else If(j==k) Then  !i=3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(j==l) Then  !i=1,2,3,4,7
               iiMax=NcR(N-3,No-2)
            Else
               iiMax=NcR(N-4,No-2)
            End If
            Do ii=1,iiMax
               I1=Abs(Ex2(ii,ij,kl))
               I2=Abs(Ex2(ii,kl,ij))
               VS=Vtmp
               If(I1/=Ex2(ii,ij,kl))VS=-VS
               Call DaxPy(Nstr,VS,X(1,I1),1,Y(1,I2),1)
            End Do
 100        Continue
         End Do; End Do
      End Do; End Do

! Opposite Spin Two-Electron Part
      ik=0
      Do i=1,N; Do k=1,N; ik=ik+1; jl=0
         If(IfZero(i,k)==1)Cycle
         If(i==k) Then
            iimax=NcR(N-1,No-1)
         Else
            iimax=NcR(N-2,No-1)
         End If
         
! Gather together elements of Xtmp
         Do ii=1,iiMax
            I1=AEx1(ii,i,k)
            VS=One; 
            If(I1/=Ex1(ii,i,k)) VS=-VS
            Do jj=1,Nstr
               Xtmp(ii,jj)=X(I1,jj)*VS
            End Do
         End Do
            
! Collect Elements of Ytmp
         Ytmp=Zero
         Do j=1,N; Do l=1,N; jl=jl+1
            If(IfSym==1) Then
               If(ik<jl) Cycle
               If(ik==jl)Vtmp=V(i,j,k,l)/Two
               If(ik>jl)Vtmp=V(i,j,k,l)
            Else
               Vtmp=V(i,j,k,l)/Two
            End If

            If(Abs(Vtmp)<1.0d-10)Cycle
 
            If(j==l) Then
               jjmax=NcR(N-1,No-1)
            Else
               jjmax=NcR(N-2,No-1)
            End If

            Do jj=1,jjMax
               J1=AEx1(jj,j,l)
               J2=AEx1(jj,l,j)
               VS=Vtmp
               If(J1/=Ex1(jj,j,l)) VS=-VS
               Call DaxPy(iiMax,VS,Xtmp(1,J1),1,Ytmp(1,J2),1)
            End Do
         End Do; End Do

! Scatter Elements of Y
         Do ii=1,iiMax
            I1=AEx1(ii,k,i)
            Do jj=1,Nstr
               Y(I1,jj)=Y(I1,jj)+Ytmp(ii,jj)
            End Do
         End Do
      End Do; End Do
! Enforce MS=0
      Y=Y+Spin*Transpose(Y)

      Return
      End Subroutine

      Subroutine GetH0(N,No,N0,N2,Max1,Max2,Nstr,Ex1,Ex2,Istr,h,V,H0)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N0,N2,Max1,Max2,Nstr,m,If0,IY(N0),mmax
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2),Istr(N0),RIstr(Nstr)
      Integer*4 i,j,k,l,ij,kl,ii,jj,I1,I2,iiMax,J1,J2,jjMax,ik,jl
      Integer*4 AEx1(Max1,N,N),AEx2(Max2,N2,N2),IfZero(N,N)
      Real*8 h(N,N),V(N,N,N,N),H0(N0,N0,N0,N0),H0tmp(N0,N0,N0,N0)
      Real*8 Vtmp,VS,VSS,htmp,hS
      Real*8 Stmp(N0)
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

      AEx1=Abs(Ex1); AEx2=Abs(Ex2); H0=Zero
! Reverse String Ordering
      RIstr=0
      Do m=1,N0
         RIstr(Istr(m))=m
      End Do

! Remove terms that are not in H0
      Do i=1,N; Do j=1,N; Do ii=1,Max1
         If0=0
         Do m=1,N0
            If(Istr(m)==AEx1(ii,i,j))If0=1
         End Do
         If(If0==0) AEx1(ii,i,j)=0
      End Do; End Do; End Do

      Do i=1,N2; Do j=1,N2; Do ii=1,Max2
         If0=0
         Do m=1,N0
            If(Istr(m)==AEx2(ii,i,j))If0=1
         End Do
         If(If0==0) AEx2(ii,i,j)=0
      End Do; End Do; End Do

! Check for Zero Blocks in V
      IfZero=0
      Do i=1,N; Do k=1,N
         Do j=1,N; Do l=1,N
            If(Abs(V(i,j,k,l))>1.0d-10)Go To 20
         End Do; End Do
         IfZero(i,k)=1
 20      Continue
      End Do; End Do

! One Electron Part
      Do i=1,N
         Do j=1,N
            htmp=h(i,j)
            If(Abs(htmp)<1.0d-10)Cycle
            If(i==j) Then
               iiMax=NcR(N-1,No-1)
            Else
               iiMax=NcR(N-2,No-1)
            End If
            Do ii=1,iiMax
               I1=AEx1(ii,i,j)
               I2=AEx1(ii,j,i)
               If(I1==0.or.I2==0) Cycle
               hS=htmp
               If(I1/=Ex1(ii,i,j))hS=-hS
!               Call DaxPy(Nstr,hS,X(1,I1),1,Y(1,I2),1)
               Do m=1,N0
                  H0(m,RIstr(I2),m,RIstr(I1))=
     $                 H0(m,RIstr(I2),m,RIstr(I1))+hS
               End Do
            End Do
         End Do
      End Do
      
! Same Spin Two-Electron Part
      ij=0
      Do i=1,N; Do j=i+1,N; ij=ij+1; kl=0
         Do k=1,N; Do l=k+1,N; kl=kl+1
            Vtmp=(V(i,j,k,l)-V(i,j,l,k)-V(j,i,k,l)+V(j,i,l,k))/Two
            If(abs(Vtmp)<1.0d-10)Go To 100 
            If(i==k.and.j==l) Then
               iiMax=NcR(N-2,No-2)
            Else If(i==k) Then  !i=2,3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(i==l) Then  !i=1-3
               iiMax=NcR(N-3,No-2)
            Else If(j==k) Then  !i=3,4*,7
               iiMax=NcR(N-3,No-2)
            Else If(j==l) Then  !i=1,2,3,4,7
               iiMax=NcR(N-3,No-2)
            Else
               iiMax=NcR(N-4,No-2)
            End If
            Do ii=1,iiMax
               I1=AEx2(ii,ij,kl)
               I2=AEx2(ii,kl,ij)
               If(I1==0.or.I2==0) Cycle
               VS=Vtmp
               If(I1/=Ex2(ii,ij,kl))VS=-VS
!               Call DaxPy(Nstr,VS,X(1,I1),1,Y(1,I2),1)
               Do m=1,N0
                  H0(m,RIstr(I2),m,RIstr(I1))=
     $                 H0(m,RIstr(I2),m,RIstr(I1))+VS
               End Do
            End Do
 100        Continue
         End Do; End Do
      End Do; End Do
      
! Opposite Spin Two-Electron Part
      ik=0
      Do i=1,N; Do k=1,N; ik=ik+1; jl=0
         If(IfZero(i,k)==1)Cycle
         If(i==k) Then
            iimax=NcR(N-1,No-1)
         Else
            iimax=NcR(N-2,No-1)
         End If
         
! Gather together phases in Stmp
         Stmp=Zero; IY=0
         Do ii=1,iiMax
            I1=AEx1(ii,i,k)
            I2=AEx1(ii,k,i)
            If(I1==0.or.I2==0)Cycle
            VS=One 
            If(I1/=Ex1(ii,i,k)) VS=-VS
            IY(RIstr(I1))=RIstr(I2)
            Stmp(RIstr(I1))=VS
            Do jj=1,Nstr
!               Xtmp(ii,jj)=X(I1,jj)*VS
            End Do
         End Do
         mmax=0
         Do m=1,N0
            If(IY(m)/=0)mmax=m
            If(IY(m)==0)IY(m)=1
         End Do
         
! Collect Elements of H0
         Do j=1,N; Do l=1,N; jl=jl+1
            If(ik<jl) Cycle
            If(ik==jl)Vtmp=V(i,j,k,l)/Two
            If(ik>jl)Vtmp=V(i,j,k,l)
            If(Abs(Vtmp)<1.0d-10)Cycle
 
            If(j==l) Then
               jjmax=NcR(N-1,No-1)
            Else
               jjmax=NcR(N-2,No-1)
            End If
            
            Do jj=1,jjMax
               J1=AEx1(jj,j,l)
               J2=AEx1(jj,l,j)
               If(J1==0.or.J2==0)Cycle
               VS=Vtmp
               If(J1/=Ex1(jj,j,l)) VS=-VS
!                  Call DaxPy(iiMax,VS,Xtmp(1,J1),1,Ytmp(1,J2),1)
               Do m=1,mmax
                  H0(IY(m),RIstr(J2),m,RIstr(J1))=
     $                 H0(IY(m),RIstr(J2),m,RIstr(J1))+VS*Stmp(m)
               End Do
            End Do
         End Do; End Do
      End Do; End Do
      Do i=1,N0; Do j=1,N0; Do k=1,N0; Do l=1,N0
         H0tmp(i,j,k,l)=H0(i,j,k,l)+H0(j,i,l,k)
      End Do; End Do; End Do; End Do
      H0=H0tmp
      Return
      End Subroutine

!      Subroutine GetIstr(N,No,Nstr,N0,Hd,Istr)
!      Integer*4 N,No,Nstr,N0,Istr(N0),itmp,jtmp,a
!      Real*8 Hd(Nstr,Nstr),Min
!
!      Do a=1,N0
!         Min=5000.0d0; ITmp=1
!         Do i=1,NStr; Do j=1,NStr
!            If(Hd(i,j)< Min) Then
!               Min=Hd(i,j)
!               ITmp=i
!            End If            
!         End Do; End Do
!      
!         Istr(a)=Itmp
!         Do i=1,NStr
!            Hd(Istr(a),i)=Hd(Istr(a),i)+10000.0d0
!         End Do
!      End Do
!
!      Do a=1,N0
!         Do i=1,NStr
!            Hd(Istr(a),i)=Hd(Istr(a),i)-10000.0d0
!         End Do
!      End Do
!
!      Return
!      End Subroutine

      Subroutine GetIstr(N,No,Nstr,N0,Hd,Istr)
      Integer*4 N,No,Nstr,N0,Istr(N0),itmp,a,b
      Real*8 Hd(Nstr,Nstr),Min

      Do a=1,N0
         Min=500.0d0; Itmp=1
         Do i=1,Nstr
            If(Hd(i,i)< Min) Then
               Min=Hd(i,i)
               Itmp=i
            End If            
         End Do
      
         Istr(a)=Itmp
         Hd(Istr(a),Istr(a))=Hd(Istr(a),Istr(a))+1000.0d0

      End Do

      Do a=1,N0
         Hd(Istr(a),Istr(a))=Hd(Istr(a),Istr(a))-1000.0d0
      End Do

      Return
      End Subroutine

      Recursive Subroutine 
     $     GetHd(N,No,Nstr,h,V,Iocca,Ioccb,Hd,IRecur)
      Use FCIUtils
      Implicit None
      Integer*4 N,No,N2,Nstr,Iocca(No),Ioccb(No),IRecur
      Integer*4 i,j,imin,jmin,k,ka,kb,l,la,lb,Isigna,Isignb
      Real*8 h(N,N),V(N,N,N,N),Hd(Nstr,Nstr),tmp

      If(IRecur==1) Then
         imin=1
         jmin=1
      Else
         imin=Iocca(IRecur-1)+1
         jmin=Ioccb(IRecur-1)+1
      End If

      Do i=imin,N
         Iocca(IRecur)=i
         Do j=jmin,N
            Ioccb(IRecur)=j

            If(IRecur==No) Then
               tmp=0.0d0
               Do k=1,No
                  ka=Iocca(k); kb=Ioccb(k)
! Spin Contaminated Elements
                  tmp=tmp+h(ka,ka)+h(kb,kb)
                  Do l=1,No
                     la=Iocca(l); lb=Ioccb(l)
                     tmp=tmp+.50d0*(
     $                    +V(ka,la,ka,la)-V(ka,la,la,ka)
     $                    +V(ka,lb,ka,lb)+V(kb,la,kb,la)
     $                    +V(kb,lb,kb,lb)-V(kb,lb,lb,kb))
                  End Do
               End Do
!               Write(6,*)'@@tmp',tmp,Iocca,Ioccb
               Hd(Index(No,Iocca),Index(No,Ioccb))=tmp
               
            Else
               Call GetHd(N,No,Nstr,h,V,Iocca,Ioccb,Hd,IRecur+1)
            End If
         End Do
      End Do

      Return
      End Subroutine


      Subroutine IString(N,No,N2,Max1,Max2,Ex1,Ex2)
      Implicit None
      Integer*4 N,No,N2,Max1,Max2
      Integer*4 Ex1(Max1,N,N),Ex2(Max2,N2,N2)
      Integer*4 i,j,k,a,b,ij,ab,IRecur,m,Ind,Iocc(No)

! Find Strings that differ by a Single Excitation
      Ex1=0
      Do i=1,N; Do a=1,N
         m=1; Ind=1
         Call RecurEx1(N,No,i,a,Max1,Iocc,Ex1,0,m,Ind)
      End Do; End Do

! Find Strings that differ by a Double Excitation
      Ex2=0; ij=0
      Do i=1,N; Do j=i+1,N; ij=ij+1; ab=0; Do a=1,N; Do b=a+1,N
         ab=ab+1; m=1; Ind=1
         Call RecurEx2(N,No,N2,i,j,a,b,ij,ab,
     $        Max2,Iocc,Ex2,0,m,Ind)
      End Do; End Do; End Do; End Do

      Return
      End Subroutine

      Recursive Subroutine 
     $     RecurEx2(N,No,N2,i,j,a,b,ij,ab,Max2,Iocc,Ex2,IRecur,m,Ind)
      Implicit None
      Integer*4 N,No,N2,i,j,a,b,ij,ab,Max2,Iocc(No),Ex2(Max2,N2,N2)
      Integer*4 IRecur,m,Ind,ii,jj,IsiO,IsjO,IsaO,IsbO
      Integer*4 Isubst(No),iind,jind,aind,bind,Isign,iimin

      If(IRecur==0)iimin=1
      If(IRecur/=0)iimin=Iocc(IRecur)+1

      Do ii=iimin,N
         Iocc(IRecur+1)=ii
         
         If(IRecur==No-1) Then
! End Recursion Condition
            IsiO=0; IsjO=0; IsaO=0; IsbO=0

! Is abj*i* |Iocc> Non-Zero?
            Do jj=1,No
               If(Iocc(jj)==i) Then
                  IsiO=1
               Else If(Iocc(jj)==j) Then
                  IsjO=1
               Else If(Iocc(jj)==a) Then
                  IsaO=1
               Else If(Iocc(jj)==b) Then
                  IsbO=1
               End If
            End Do
            If(IsiO==1 .and. IsjO==1 .and. IsaO==0 .and. IsbO==0) Then
               Isubst=Iocc
               Do jj=1,No
                  If(Isubst(jj)==i) Then
                     Isubst(jj)=a
                  Else If(Isubst(jj)==j) Then
                     Isubst(jj)=b
                  End If
               End Do

               Call ISort(No,Isubst,Isign)

               Ex2(m,ab,ij)=Isign*Ind
               m=m+1
            End If
            Ind=Ind+1
         Else
            Call RecurEx2(N,No,N2,i,j,a,b,ij,ab,
     $           Max2,Iocc,Ex2,IRecur+1,m,Ind)
         End If

      End Do

      Return
      End Subroutine

      Recursive Subroutine 
     $     RecurEx1(N,No,i,a,Max1,Iocc,Ex1,IRecur,m,Ind)
      Implicit None
      Integer*4 N,No,i,a,Max1,Iocc(No),Ex1(Max1,N,N),IRecur,m,Ind,ii
      Integer*4 jj,IsiO,IsaO,Isubst(No),aind,iind,Isign,iimin

      If(IRecur==0)iimin=1
      If(IRecur/=0)iimin=Iocc(IRecur)+1

      Do ii=iimin,N
         Iocc(IRecur+1)=ii
         
         If(IRecur==No-1) Then
! End Recursion Condition
            IsiO=0; IsaO=0
            Do jj=1,No
               If(Iocc(jj)==i) Then 
                  IsiO=1
               Else If(Iocc(jj)==a) Then
                  IsaO=1
               End If
            End Do

            If(IsiO==1 .and. IsaO==0) Then
               Isubst=Iocc
               Do jj=1,No
                  If(i==Isubst(jj)) Isubst(jj)=a
               End Do

               Call Isort(No,Isubst,Isign)

               Ex1(m,a,i)=Isign*Ind
               m=m+1
            End If
            Ind=Ind+1
         Else
            Call RecurEx1(N,No,i,a,Max1,Iocc,Ex1,IRecur+1,m,Ind)
         End If

      End Do

      Return
      End Subroutine

      Subroutine GSOrtho(N,C)
! Takes a set of Non-Orthogonal, Non-Normalized vectors and makes
!   them orthonormal by the Gramm-Schmidt procedure
!   Vectors are C(:,p)
      Implicit None
      Integer*4 N,p,q
      Real*8 C(N,N),pNorm,DDot

      Do p=1,N
         Do q=1,p-1
            C(:,p)=C(:,p)-DDot(N,C(:,p),1,C(:,q),1)*C(:,q)
         End Do
         pNorm=0.0d0
         Do q=1,N
            pNorm=pNorm+C(q,p)*C(q,p)
         End Do

         C(:,p)=C(:,p)/Sqrt(pNorm)
      End Do

      Return
      End Subroutine

      Subroutine ISort(N,X,Isign)
!  Sorts an Integer arrary X by straight insertion
!     Isign is the sign of the Permutation needed to bring it in order
      Implicit None
      Integer*4 N,X(N),i,j,tmp,Isign

      Isign=1
      Do i=2,N
         tmp=X(i)
         Do j=i-1,1,-1
            If(X(j).le.tmp)Go To 100
! Transposing indicies ...
            X(j+1)=X(j)
! Gives us a minus sign.
            Isign=-Isign
         End Do
         j=0
 100     X(j+1)=tmp
      End Do

      Return
      End Subroutine
      Subroutine Sort(N,X,Isign)
!  Sorts a Real arrary X by straight insertion
!     Isign is the sign of the Permutation needed to bring it in order
      Implicit None
      Integer*4 N,i,j,Isign
      Real*8 X(N),tmp

      Isign=1
      Do i=2,N
         tmp=X(i)
         Do j=i-1,1,-1
            If(X(j).le.tmp)Go To 100
! Transposing indicies ...
            X(j+1)=X(j)
! Gives us a minus sign.
            Isign=-Isign
         End Do
         j=0
 100     X(j+1)=tmp
      End Do

      Return
      End Subroutine

      Subroutine GetEig2(H,E,U)
      Implicit None
! Diagonalize a 2x2 in High Precision
      Real*8 H(2,2),E(2),U(2,2)
      Real(Kind=10) HQ(2,2),EQ(2),UQ(2,2)
      Real(Kind=10) Two,fac

      Two=2.0q0
      HQ=H
      
      fac=(HQ(1,1)-HQ(2,2))*(HQ(1,1)-HQ(2,2))/4.0q0 + HQ(1,2)*HQ(1,2)
      fac=sqrt(fac)
      EQ(1)=(HQ(1,1)+HQ(2,2))/Two-fac
      EQ(2)=(HQ(1,1)+HQ(2,2))/Two+fac
      
      UQ=0q0; UQ(1,1)=1q0; UQ(2,2)=1q0
      If(Abs(EQ(1)-HQ(1,1)).gt.1d-12) Then
         fac= 1.0q0 + HQ(1,2)/(HQ(1,1)-EQ(1))*HQ(1,2)/(HQ(1,1)-EQ(1))
         fac=sqrt(fac)
         
         UQ(1,1)= -HQ(1,2)/fac/(HQ(1,1)-EQ(1))
         UQ(2,1)= 1.0q0/fac
         
         UQ(1,2)=-UQ(2,1)
         UQ(2,2)=UQ(1,1)
      End If

      E=EQ; U=UQ

      Return 
      End
