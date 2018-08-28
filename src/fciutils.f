      Module FCIUtils
      Contains
      Function MDot(X,Y)
      Implicit None
      Real*8 MDot,X(:,:),Y(:,:)
      Integer*4 i,j
      
      If(Size(X,1)/=Size(Y,1) .or. Size(X,2)/=Size(Y,2)) 
     $     Write(6,*) 'Array Dimensions incompatible in MDot',
     $     Size(X),Size(Y)      

      MDot=0.0d0
      Do i=1,Size(X,1)
         Do j=1,Size(X,2)
            MDot=MDot+X(i,j)*Y(i,j)
         End Do
      End Do

      End Function


      Function Index(No,Iocc)
!
! Returns the correct index for a _pre-sorted_ String
!
      Implicit None
      Integer*4 Index,No,Iocc(No),i,Isign
      Common /ZIndex/ ZIndex
      Integer*4,Pointer :: ZIndex(:,:)

      Index=1
      Do i=1,No
         Index=Index+ZIndex(i,Iocc(i))
      End Do

      End Function

      Function Factorial(N)
      Implicit None
      Integer*4 N,Factorial,i

      Factorial=1
      Do i=2,N
         Factorial=Factorial*N
      End Do

      End Function

      Function NcR(N,R)
! Uses a temporary integer*8 to make sure it doesn't overflow
!  Works for N!/R! < 2^64
      Implicit None
      Integer*4 NcR,N,R,i
      Integer*8 Ntmp

      If(R>N .or. R<0) Then
         NcR=0; Return
      End If
      
      Ntmp=1
      Do i=N,N-R+1,-1
         Ntmp=Ntmp*i
      End Do

      Do i=2,R
         Ntmp=Ntmp/i
      End Do

      NcR=Ntmp

      End Function

      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

      End Module
