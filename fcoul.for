c    **************************************************************
        real*4 function fcoul(z,e,a)
        implicit none
        real*4 gmma
        complex clngam
        complex zi,clg
	real*4 z1,z2
        real*4 xke,e,a,z,rad,gam1,eta,alpha,xme,rat
	real*4 expo
c	real*4 rad1
c     calculates coulomb function
c     translate energy to MeV
	e=xme*e
        data alpha/.007297203/,xme/.510999/
        xke=sqrt(e*e-xme*xme)
        rad=(1.2/197.3)*a**(1./3.)
c        rad=(.0056908*a**(1./3.)+.01191977*a**(-1./3.)-.01049119/a)
        gam1=sqrt(1.-(alpha*z)**2)
        eta=z*alpha*e/xke
        zi=(0.,1.)
        clg=clngam(gam1+zi*eta)
        z1=cabs(cexp(clg))
        z2=gmma(2*gam1+1.)
        rat=z1/z2
	expo=exp(3.1415926*eta)
        fcoul=rat*rat*4.*(2.*xke*rad)**(2.*gam1-2.)*expo
        return
        end
c    **************************************************************
      function gmma(y)
	implicit none
      double precision c
	real*4 y,f,gmma,z
	integer n,ii,i
      dimension c(14)
      data c/                          0.99999 99999 99990 44d0,
     1   0.42278 43351 02334 79d0,     0.41184 03301 66781 29d0,
     2   0.08157 69261 24155 46d0,     0.07424 89154 19444 74d0,
     3  -0.00026 61865 94953 06d0,     0.01114 97143 35778 93d0,
     4  -0.00283 64625 30372 82d0,     0.00206 10918 50225 54d0,
     5  -0.00083 75646 85135 17d0,     0.00037 53650 52263 07d0,
     6  -0.00012 14173 48706 32d0,     0.00002 79832 88993 83d0,
     7  -0.00000 30301 90810 28d0/
      n=0
      f=1.
    1 if (y+n-2.) 2,5,3
    2 f=f*(y+n)
      n=n+1
      go to 1
    3 if (y+n-3.) 6,4,4
    4 n=n-1
      f=f*(y+n)
      go to 3
    5 gmma=1.
      go to 8
    6 gmma=c(14)
      z=y+n-2
      do 7 ii=1,13
      i=13-ii+1
    7 gmma=gmma*z+c(i)
    8 if (n) 11,12,9
    9 if (abs(f).lt.1.e-30) go to 10
      gmma=gmma/f
      go to 12
   10 gmma=1.e30
      go to 12
   11 gmma=gmma*f
   12 return
      end
c  **********************************************************
      COMPLEX FUNCTION CLNGAM(X)
C
C          COMPUTES LN GAMMA(X) TO ABOUT 14 FIGURES
C                      X COMPLEX,  MAX(7, RE X)  GE  ABS(IM X)
C
	implicit none
	real*8 C,D,A1,A2,A3,A4,A5,A6,A7
        COMPLEX X,Y,G,Z
c
      Y=X
      G=1.
   21 IF(REAL(Y).GT.7.)  GO TO 23
      G=G*Y
      Y=Y+1.
      GO TO 21
   23 Z=-Y**(-2)
      C=.39894 22804 01433 D0
      D=.29550 65359 47712 D-1
      A1=.28200 16588 33287 D1
      A2=.94000 55294 44291 D-1
      A3=.26857 30084 12655 D-1
      A4=.20142 97563 09491 D-1
      A5=.28485 01604 37664 D-1
      A6=.64889 49259 20085 D-1
      A7=.21692 43529 48683 D0
      CLNGAM = (Y-.5)*CLOG(Y)-Y-CDLOG(C*G)+(A1+Z*(A2+Z*(A3+Z*(A4+Z*(A5
     1   +Z*(A6+Z*(A7+Z)))))))*D/Y
      RETURN
      END
