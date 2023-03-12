C     EQUILIBRIUM RECONSTRUCTION
C     0!
C     Lluís Torres 
C     TFG
C     FORTRAN 2003

      PROGRAM EQUILIBRIUM_RECONSTRUCTION
            
      USE SPIN_MODEL

C-------------------------------------------------------------
C     NODES, EDGES, CONNECTIVITY
      INTEGER N, M, z
C     SAMPLE SIZE (# OF SPIN CONFIGURATIONS)
      INTEGER C
C     TEMPERATURE (TEMP = k_B·T)
      REAL*8 TEMP
C-------------------------------------------------------------
C     NUMBER OF GRAPHS TO SIMULATE FOR EVERY P VALUE
      INTEGER NSEEDS
C     SPIN CONFIGURATION SAVING VARIABLES
C     (STORE SPIN CONFIGURATION AS N/zip_size INTEGERS)
      INTEGER zip_size
      CHARACTER(:), ALLOCATABLE:: bin
      INTEGER, ALLOCATABLE:: decimal(:)
C-------------------------------------------------------------
C     MAXIMUM NUMBER OF NEIGHBORS A NODE CAN HAVE
      INTEGER zmax
c     NUMBER OF MONTE-CARLO STEPS FOR RECONSTRUCTION
      INTEGER TAU
C-------------------------------------------------------------
C     (FOR PLOTING PL AND GAMMA AS A FUNCTION OF TIME)
C     INITIAL MONTE-CARLO AVERAGE STEP
      INTEGER MCINI2
      parameter (MCINI2 = 1)
C     RATE OF MONTE-CARLO AVERAGES
      INTEGER MCD2
      parameter (MCD2 = 1)
C-------------------------------------------------------------
C     +1 -1 EDGES RATIO (1 => all +1), (0 => all -1)
      REAL*8 p
C-------------------------------------------------------------
      INTEGER i, IMC, IPAS
C     SEED NUMBER
      INTEGER SEED, SEEDini
C     PRINCIPAL ARRAYS
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR_0(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR_0(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ_0(:)

      INTEGER, ALLOCATABLE:: S_SET(:,:)
      
      LOGICAL valid
      REAL*8 TF_STEP
      REAL*8 PL, DPL, SUMPL, SUMg
      REAL*4 TIME1, TIME2, time
C     FICTICIOUS TEMPERATURE (TEMP_F = k_B·TF)
      REAL*8 TEMP_F
C-------------------------------------------------------------
C     NUMBER OF p VALUES TO SIMULATE
      INTEGER npvalues, np
      REAL*8 pvalues(1:1)
C-------------------------------------------------------------
      CHARACTER(4) str
      CHARACTER(3) str1, str2

      INTEGER, ALLOCATABLE:: LAMBDA(:,:)
      REAL*8, ALLOCATABLE:: funct(:,:)
C-------------------------------------------------------------
C     (DUMMY VARIABLE)
      INTEGER SC
C-------------------------------------------------------------

C-----------------------------------------------------------------------
C     START
C-----------------------------------------------------------------------

      WRITE(*,*) '>>> EQUILIBRIUM RECONSTRUCTION <<<'

C     READ SIMULATION VARIABLES FROM INPUT FILE
      CALL READ_INPUT(N,z,TEMP,pvalues,C,NSEEDS,SC,zip_size,TAU)
      M = z*N/2
      zmax = N-1
C     ALLOCATION 
      ALLOCATE(decimal(1:N/zip_size))
      bin = REPEAT(' ',N)
      ALLOCATE(LAMBDA(1:C,1:N))
      ALLOCATE(funct(-1:1,-zmax:zmax))

      SEEDini = 100
      
C     # of p VALUES
      npvalues = size(pvalues)

      call cpu_time(TIME1)

C     DO FOR ALL p VALUES
      DO np = 1,npvalues

      p = pvalues(np)

      ALLOCATE(S_SET(1:C,1:N))

      WRITE(str,'(f4.2)') p 
      str1 = str(1:1)//str(3:4)

c     GAMMA DATA FILE
      OPEN(UNIT=10,FILE='g_'//str1//'.dat')

C     DO FOR ALL SEEDS
      DO SEED = SEEDini,SEEDini+NSEEDS-1

      WRITE(str2,'(i3)') SEED

C	GENERATE THE ORIGINAL GRAPH FROM METROPOLIS ALGORITHM
      call init_genrand(SEED) 

C     ORIGINAL SYSTEM
      call IRS(N,M,p,NBR_0,INBR_0,JJ_0)

C     ORIGINAL GRAPH
      !CALL IRG(N,M,NBR_0,INBR_0,JJ_0)
C     ORIGINAL COUPLING ASSIGMENT
      !CALL RCA(N,M,p,NBR_0,INBR_0,JJ_0)

C     READ THE SAMPLE
      OPEN(UNIT=1,FILE='S_'//str1//'_'//str2//'.bin',FORM='UNFORMATTED')
      DO i = 1,C
            READ(1) decimal
            CALL DEC2BIN(N,zip_size,bin,decimal)
            CALL BIN2ARRAY(N,bin,S_SET(i,:))
      END DO
      CLOSE(1)

C     INITIAL FICTICIOUS TEMPERATURE  
      TEMP_F = -log(0.5d0*(1+tanh(1./TEMP)))/z
      TF_STEP = TEMP_F/(TAU)
      
C     INITIAL RANDOM SYSTEM
      CALL init_genrand(seed)

C     INITIAL RANDOM SYSTEM
      call IRS(N,M,p,NBR,INBR,JJ)

c     GENERATE SAME GRAPH STRUCTURE
      !CALL init_genrand(SEED)
      !CALL IRG(N,M,NBR,INBR,JJ)
c     WITH RANDOM COUPLING ASSIGMENT
      !CALL init_genrand(999)
      !CALL RCA(N,M,p,NBR,INBR,JJ)

      CALL CLASS_FUNCTION(zmax,TEMP,funct)
      CALL CLASS_LAMBDA(N,C,S_SET,NBR,JJ,LAMBDA)

C     INITIAL PSEUDOLIKELIHOOD
      PL = PSEUDO(N,C,S_SET,TEMP,NBR,JJ)

      SUMPL = 0
      SUMg = 0

C     MONTE-CARLO SIMULATION
      DO IMC = 1,TAU

            DO IPAS = 1,M
            CALL PSEUDOLIKELIHOOD(N,C,S_SET,valid,TEMP_F,
     .                        DPL,NBR,INBR,JJ,funct,zmax,LAMBDA)
            IF (valid) THEN
                  PL = PL + DPL
            END IF
            END DO

            SUMPL = SUMPL + PL
            SUMg = SUMg + GAMMAA(N,M,NBR,JJ,NBR_0,JJ_0)
C           AVERAGE VALUES EVERY MCD2 MONTE-CARLO STEPS
            IF ((IMC.GT.MCINI2).AND.(MCD2*(IMC/MCD2).EQ.IMC)) THEN      
                  SUMPL = 0
                  SUMg = 0
            END IF

            TEMP_F = TEMP_F - TF_STEP

      ENDDO

      WRITE(10,*) SEED, GAMMAA(N,M,NBR,JJ,NBR_0,JJ_0)

c     SAVE JJ
      CALL SAVE_COUPLINGS(N,str1,str2,NBR,JJ)

      DO i = 1,N
            DEALLOCATE(NBR(i)%v)
            DEALLOCATE(INBR(i)%v)
            DEALLOCATE(JJ(i)%v)
      END DO
      DEALLOCATE(NBR)
      DEALLOCATE(INBR)
      DEALLOCATE(JJ)

      DO i = 1,N
            DEALLOCATE(NBR_0(i)%v)
            DEALLOCATE(INBR_0(i)%v)
            DEALLOCATE(JJ_0(i)%v)
      END DO
      DEALLOCATE(NBR_0)
      DEALLOCATE(INBR_0)
      DEALLOCATE(JJ_0)

200   FORMAT (A,I3,A,I3,A,I3,A,I3,A,I3,A,I3)
      IF ((SEED.EQ.SEEDini).AND.(p.EQ.pvalues(1))) THEN
      CALL CPU_TIME(TIME2)
      time = (TIME2-TIME1)*npvalues*NSEEDS
      WRITE(*,200) "ESTIMATED TIME: ", INT(time/3600), ' h',
     . INT((time/3600-int(time/3600))*60), ' min', 
     . INT((time/60-int(time/60))*60), ' s'
      END IF

      END DO !SEED

      DEALLOCATE(S_SET)

      END DO !p

      CALL CPU_TIME(TIME2)
      time = (TIME2-TIME1)
      WRITE(*,200) "CPU TIME: ", INT(time/3600), ' h',
     . INT((time/3600-int(time/3600))*60), ' min', 
     . INT((time/60-int(time/60))*60), ' s'

      CLOSE(10)

      END PROGRAM EQUILIBRIUM_RECONSTRUCTION

c------------------------------------------------------------------
c  A C-program for MT19937, with initialization improved 2002/1/26.
c  Coded by Takuji Nishimura and Makoto Matsumoto.
c
c  Before using, initialize the state by using init_genrand(seed)  
c  or init_by_array(init_key, key_length).
c
c  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
c  All rights reserved.                          
c  Copyright (C) 2005, Mutsuo Saito,
c  All rights reserved.                          
c
c  Redistribution and use in source and binary forms, with or without
c  modification, are permitted provided that the following conditions
c  are met:
c
c    1. Redistributions of source code must retain the above copyright
c       notice, this list of conditions and the following disclaimer.
c
c    2. Redistributions in binary form must reproduce the above copyright
c       notice, this list of conditions and the following disclaimer in the
c       documentation and/or other materials provided with the distribution.
c
c    3. The names of its contributors may not be used to endorse or promote 
c       products derived from this software without specific prior written 
c       permission.
c
c  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
c  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c
c
c  Any feedback is very welcome.
c  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
c  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
c
c-----------------------------------------------------------------------
c  FORTRAN77 translation by Tsuyoshi TADA. (2005/12/19)
c
c     ---------- initialize routines ----------
c  subroutine init_genrand(seed): initialize with a seed
c  subroutine init_by_array(init_key,key_length): initialize by an array
c
c     ---------- generate functions ----------
c  integer function genrand_int32(): signed 32-bit integer
c  integer function genrand_int31(): unsigned 31-bit integer
c  double precision function genrand_real1(): [0,1] with 32-bit resolution
c  double precision function genrand_real2(): [0,1) with 32-bit resolution
c  double precision function genrand_real3(): (0,1) with 32-bit resolution
c  double precision function genrand_res53(): (0,1) with 53-bit resolution
c
c  This program uses the following non-standard intrinsics.
c    ishft(i,n): If n>0, shifts bits in i by n positions to left.
c                If n<0, shifts bits in i by n positions to right.
c    iand (i,j): Performs logical AND on corresponding bits of i and j.
c    ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
c    ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
c
c-----------------------------------------------------------------------
c     initialize mt(0:N-1) with a seed
c-----------------------------------------------------------------------
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
c
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*
     &          ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
c
      return
      end
c-----------------------------------------------------------------------
c     initialize by an array with array-length
c     init_key is the array for initializing keys
c     key_length is its length
c-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
c
      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)
     &           +init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
c
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0xffffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
c
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
c
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
c
      y=mt(mti)
      mti=mti+1
c
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
c
      genrand_int32=y
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0x7fffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1]-real-interval
c-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on (0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1) with 53-bit resolution
c-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
c-----------------------------------------------------------------------
c     initialize large number (over 32-bit constant number)
c-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
CC    TOPBIT_MASK = Z'80000000'
CC    ALLBIT_MASK = Z'ffffffff'
CC    UPPER_MASK  = Z'80000000'
CC    LOWER_MASK  = Z'7fffffff'
CC    MATRIX_A    = Z'9908b0df'
CC    T1_MASK     = Z'9d2c5680'
CC    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end