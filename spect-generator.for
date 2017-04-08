c  comes from spectrum-2.for
c  calculates expected spectrum, includes statistical fluctuations; then it uses it as data;
c  fits for overall normalization and bFierz.
c  does not include wall effect on efficiency.
c  includes gaussian spread via sigmaEovEinp,sigmaEovEanlys
c  or alternatively a low-energy tail.
	implicit none
	integer Nmax
	integer Nfreq
	integer NBfield
	parameter(Nmax=10000)
	parameter(Nfreq=6)
	parameter(NBfield=12)
	real*8 me
	real*8 clight
	real*8 TotArea
	real*8 Pi
	real*8 Bfield,DeltaBfield
	real*8 freq2energ
	real*8 area
	real*8 E(Nmax)
	real*8 E0
	real*8 dNdE(Nmax),DeltadNdE(Nmax)
	real*8 dNdEgas(Nmax)
	real*8 standard(Nmax),fierz(Nmax)
	real*8 standgas(Nmax),fierzgas(Nmax)
	real*8 resid(Nmax)
	real*8 fmin(Nfreq),fmax(Nfreq)
	real*8 Emin(Nfreq),Emax(Nfreq)
	real*8 spectrum_data,spectrum_fit,spectrum2
	real*8 spectrum_WM
	real*8 chisq,chisqPerNu
	real*8 fluct
	real*4 gasdev,ran4
	real*8 sum11,sum12,sum22,sum1y,sum2y
	real*8 covar(2,2),beta(2)
	real*8 bFierz,dbFierz
	real*8 relativeerrorSq
	real*8 bFinput,bWMinput
	real*8 bWMdistort,wrongnessWM
	real*8 areadata,areatot,funct,fract
	real*8 radius
	real*8 wrongnessEgain,wrongnessEshift
	real*8 Enrg
	real*8 sigmaEovEinp,sigmaEovEanlys
	real*8 lambdatail,fractail
	real*8 counts
	real*8 Ejunk,random,maximum
	parameter(wrongnessEgain=1.0d0,wrongnessEshift=0.0d0,wrongnessWM=1.00d0)
	real*4 fcoul
	integer n1,i,nb
	integer idum,iseed1
	integer idum0(NBfield)
	integer Nqvalue
	integer ncounter
	integer NmaxHi
	integer ibinjunk,NjunkMax
	logical flag_inwindow,flag_ignore
	character*14 file_aux_1
	character*13 file_aux_2
	character*15 file_aux_3

	data idum0/-512863494,-152863494,-112863494,-121683494,-158263494,
     &             -128613494,
     &-225863494,-252863494,-312863494,-131683494,-758263494,-928613494/
	common/bfactors/bFinput,bWMinput
c
	flag_ignore=.true.
	me=5.11d5
	E0=3508.d3/me+1.d0
c	E0=2216.4d3/me				!19Ne
c	E0=1808.44d3/me+1.d0			!14O
	Nqvalue=3500
	clight=3.0d8
	TotArea=20.d9
	Pi=3.14159265359d0
	bFinput=+0.02d0
	bWMinput=1.41d-3
	bWMdistort=bWMinput*wrongnessWM
	DeltaBfield=0.5d0
	sigmaEovEinp=1.0d-3
	sigmaEovEanlys=1.0d-3
	NjunkMax=100000
	fractail=0.006d0
	lambdatail=14.d0
        open(unit=50,file='events-vs-B.out',status='unknown')
        open(unit=20,file='bFierz-vs-B.out',status='unknown')
	do nb=1,NBfield
c          open(unit=35,name='/dev/urandom',access='keyed',
c     &         form='unformatted',recl=1)
c          read(35)iseed1
c	  idum=-iabs(iseed1)
c	  write(*,*)'idum',idum
c          close(unit=35)
	idum=idum0(nb)
	Bfield=1.0+float(nb-1)*DeltaBfield
	freq2energ=clight*clight*Bfield/(2*Pi*1.d9*me)		!GHz to 1/me
        write(file_aux_1,'(''spec-B'',f4.2,''.dat'')')Bfield
        write(file_aux_2,'(''fit-B'',f4.2,''.dat'')')Bfield
        write(file_aux_3,'(''resid-B'',f4.2,''.dat'')')Bfield
        open(unit=30,file=file_aux_1,status='unknown')
        open(unit=40,file=file_aux_2,status='unknown')
        open(unit=60,file=file_aux_3,status='unknown')
c set frequencies (in GHz) and corresponding energies
c reset for 6 GHz coverage starting at 18 GHz
	fmin(1)=18.00
	do n1=2,Nfreq
	  fmin(n1)=fmin(n1-1)+1.d0
	enddo
	do n1=1,Nfreq
	  fmax(n1)=fmin(n1)+1.d0
	  Emin(n1)=freq2energ/fmax(n1)
	  Emax(n1)=freq2energ/fmin(n1)
	enddo
	write(*,*)fmin
c  in preparation for smearing calculate energies over whole range of Nmax
	do i=1,Nmax
	  E(i) = 1.d3*float(i)/me+1.d0		!in 1/me
	enddo
c  now calculate spectrum over up to qvalue
	do i=1,Nqvalue
	  Enrg=E(i)*wrongnessEgain+wrongnessEshift
	  if(Enrg .ge. E0)Enrg=E0
	  dNdE(i) = spectrum_data(Enrg,E0)
	enddo
c  to check smearing: produce delta fn.
c	if(nb .eq. 1)then
c	  do i=1,2*Nqvalue
c	    dNdE(i)=0.d0
c	  enddo
c	  dNdE(3500)=1.d0
c	endif
c	if(nb .eq. 1)then
c	  do i=1,Nmax
c	    write(35,*)E(i),dNdEgas(i)
c	  enddo
c	endif
c	pause
c  call gaussian smearing (if using check using dNdEgas below)
	if(sigmaEovEinp .ne. 0.d0)then
	  call gas_smearing(dNdE,E,Nmax,Nqvalue,sigmaEovEinp,dNdEgas)
	else
	  do i=1,Nqvalue
	    dNdEgas(i)=dNdE(i)
	  enddo
	endif
c  get area of dNdEgas array
	NmaxHi=1
	area=0.d0
	do i=1,Nmax
	  counts=dNdEgas(i)
	  if(counts .gt. 0.d0)then
	    NmaxHi=i
	    area=area+counts
	  endif
	enddo
	write(*,*)'NmaxHi',NmaxHi
	write(*,*)'area after gauss',area
c  renormalize so area=TotArea and get maximum
	maximum=0.d0
	do i=1,Nmax
	  dNdEgas(i)=dNdEgas(i)*TotArea/area 
	  if(dNdEgas(i) .gt. maximum)maximum=dNdEgas(i)
	enddo
c  randomize spectrum
	do i=1,NmaxHi
	  DeltadNdE(i)=dsqrt(dNdEgas(i))
	  if(DeltadNdE(i) .lt. 1.d0)DeltadNdE(i)=1.d0
	  fluct=DeltadNdE(i)*dble(gasdev(idum))
	  dNdEgas(i)=dNdEgas(i)+fluct
	enddo
c  now distribute an additional Njunk counts from a random place in the spectrum 
c  above the energy window covered by the frequency window.
c  get random energy (between Emax(Nfreq) and endpoint)
	do i=1,NjunkMax
	  Ejunk=Emax(Nfreq)+(E(Nqvalue)-Emax(Nfreq))*ran4(idum)
	  random=ran4(idum)*maximum
	  if(random .lt. maximum)then
c  take this and add it to the spectrum within the frequency band of interest
	    Ejunk=Emin(1)+ran4(idum)*(Emax(Nfreq)-Emin(1))
	    ibinjunk=(Ejunk-1.d0)*me/1.d3
	    dNdEgas(ibinjunk)=dNdEgas(ibinjunk)+1.d0  
	  endif
	enddo
c  write out spectrum
	areadata=0.d0
	do i=1,NmaxHi
	  flag_inwindow=.false.
	  do n1=1,Nfreq
	    if(E(i) .lt. Emax(n1) .and. E(i) .ge. Emin(n1))then
	      flag_inwindow=.true.
	    else
	    endif
	  enddo
	  if(flag_inwindow .eqv. .true.)then   
	    radius=100.d0*dsqrt(E(i)**2-1.d0)*me/(clight*Bfield)	!in cm
	    write(30,*)E(i)-1.d0,dNdEgas(i),DeltadNdE(i),radius
	    areadata=areadata+dNdEgas(i)
	  else
	  endif
	enddo
c
c
c compute matrix elements for minimization
c  first compute spectrum for analysis
	do i=1,Nqvalue
	  standard(i)=spectrum_fit(E(i),E0)+bWMdistort*spectrum_WM(E(i),E0)
	  fierz(i)=spectrum2(E(i),E0)
	enddo
c  call gaussian smearing (if using check using dNdEgas below)
	if(sigmaEovEanlys .ne. 0.d0)then
	  call gas_smearing(standard,E,Nmax,Nqvalue,sigmaEovEanlys,standgas)
	  call gas_smearing(fierz,E,Nmax,Nqvalue,sigmaEovEanlys,fierzgas)
	else
	  do i=1,Nqvalue
	    standgas(i)=standard(i)
	    fierzgas(i)=fierz(i)
	  enddo
	endif
	sum11=0.d0
	sum12=0.d0
	sum22=0.d0
	sum1y=0.d0
	sum2y=0.d0
	do i=1,NmaxHi
	  flag_inwindow=.false.
	  do n1=1,Nfreq
	    if(E(i) .lt. Emax(n1) .and. E(i) .ge. Emin(n1))then
	      flag_inwindow=.true.
	    else
	    endif
	  enddo
	  if(flag_inwindow .eqv. .true.)then
	    sum11=sum11+(standgas(i)/DeltadNdE(i))**2.d0
	    sum12=sum12+standgas(i)*fierzgas(i)/DeltadNdE(i)**2.d0
	    sum22=sum22+(fierzgas(i)/DeltadNdE(i))**2.d0
	    sum1y=sum1y+standgas(i)*dNdEgas(i)/DeltadNdE(i)**2.d0
	    sum2y=sum2y+fierzgas(i)*dNdEgas(i)/DeltadNdE(i)**2.d0
	  else
	  endif
	enddo
c produce matrix for inversion
	covar(1,1)=sum11
	covar(1,2)=sum12
	covar(2,1)=sum12
	covar(2,2)=sum22
	beta(1)=sum1y
	beta(2)=sum2y
	if(covar(1,1) .eq. 0.d0)goto 400
	call gaussj(covar,2,2,beta,1,1)
	bFierz=beta(2)/beta(1)
	relativeerrorSq=covar(1,1)/beta(1)**2.0d0+covar(2,2)/beta(2)**2.0d0
	dbFierz=dabs(bFierz)*dsqrt(relativeerrorSq)
	areatot=0.d0
	chisq=0.d0
	ncounter=0
	do i=1,Nqvalue
	  funct=beta(1)*standgas(i)+beta(2)*fierzgas(i)
	  write(40,*)E(i)-1,funct
	  areatot=areatot+funct
	  flag_inwindow=.false.
	  do n1=1,Nfreq
	    if(E(i) .lt. Emax(n1) .and. E(i) .ge. Emin(n1))then
	      flag_inwindow=.true.
	    else
	    endif
	  enddo
	  if(flag_inwindow .eqv. .true.)then
	    resid(i)=(funct-dNdEgas(i))/DeltadNdE(i)
	    write(60,*)E(i)-1,resid(i)
	    chisq=chisq+resid(i)**2
	    ncounter=ncounter+1
	  else
	  endif
	enddo
	write(*,*)Bfield,areatot,areadata
	close(unit=30)
	close(unit=40)
	close(unit=60)
	fract=areadata/areatot
	chisqPerNu=chisq/float(ncounter)
	write(50,500)Bfield,TotArea*(dbFierz/1.d-3)**2/1.d9,fract,chisqPerNu
	write(20,510)Bfield,bFierz,dbFierz,chisqPerNu
400	enddo
	close(unit=50)
500	format(1x,f4.2,2x,f7.4,2x,f7.4,2x,f9.4)
510	format(1x,f4.2,2x,f9.6,2x,f9.6,2x,f9.4)
	stop
	end
c   ************************************************************************************
	real*8 function spectrum_data(E,E0)
	real*8 E
	real*8 E0
	real*8 bFinput,bWMinput
	real*4 fcoul
	common/bfactors/bFinput,bWMinput
	spectrum_data=fcoul(3.,real(E),6.)*(E-E0)*(E-E0)*E*sqrt(E*E-1.d0)*
     &		(1.d0+bFinput/E+bWMinput*(2.d0*E-E0-1.d0/E))
	return
	end
c
	real*8 function spectrum_fit(E,E0)
	implicit none
	real*8 E
	real*8 E0
	real*4 fcoul
	spectrum_fit=fcoul(3.,real(E),6.)*(E-E0)*(E-E0)*E*sqrt(E*E-1.d0)
	return
	end
c
	real*8 function spectrum2(E,E0)
	implicit none
	real*8 E
	real*8 E0
	real*4 fcoul
	spectrum2=fcoul(3.,real(E),6.)*(E-E0)*(E-E0)*sqrt(E*E-1.d0)
	return
	end
c
	real*8 function spectrum_WM(E,E0)
	implicit none
	real*8 E
	real*8 E0
	real*4 fcoul
	spectrum_WM=fcoul(3.,real(E),6.)*(E-E0)*(E-E0)*E*sqrt(E*E-1.d0)*
     &		(2.d0*E-E0-1.d0/E)
	return
	end
c  *****************************************************************************
	subroutine gas_smearing(dNdE,E,Nmax,Nqvalue,sigmaEovE,dNdEgas)
	integer Nmax
	integer Nqvalue
	integer i
	real*8 dNdE(Nmax),dNdEgas(Nmax)
	real*8 E(Nmax)
	real*8 sigmaEovE
	real*8 sigma,exponent,deltaE
	do i=1,Nmax
	  dNdEgas(i)=0.d0
	  do j=1,Nqvalue
	    deltaE=E(i)-E(j)
	    sigma=sigmaEovE*E(j)
	    if(sigma .le. 0.d0)sigma=0.01
	    exponent=((deltaE/sigma)**2.d0)/2.d0
	    dNdEgas(i)=dNdEgas(i)+dNdE(j)*dexp(-exponent)
	  enddo
	enddo
	return
	end
c
