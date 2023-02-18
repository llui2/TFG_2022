C     SPIN MODEL MODULE
C     0!
C     Lluís Torres 
C     TFG
C     FORTRAN 2003

      MODULE SPIN_MODEL 

C     MULTI ARRAY TYPE
      TYPE :: MULTI_ARRAY
      INTEGER,ALLOCATABLE :: v(:)
      END TYPE MULTI_ARRAY

      CONTAINS

C-----------------------------------------------------------------------
C     METROPOLIS.F
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     READ INPUT FILE
      SUBROUTINE READ_INPUT(N,z,TEMP,pvalues,C,NSEEDS,SC,zip_size,TAU)

      INTEGER N,z
      REAL*8 TEMP
      REAL*8 pvalues(1:7)
      INTEGER C
      INTEGER NSEEDS
      INTEGER SC
      INTEGER zip_size
      INTEGER TAU

      OPEN(UNIT=0,FILE="input.txt")
      
      READ(0,*)
      READ(0,*) N,z
      READ(0,*)
      READ(0,*) TEMP
      READ(0,*)
      READ(0,*) pvalues
      READ(0,*)
      READ(0,*) C
      READ(0,*)
      READ(0,*) NSEEDS
      READ(0,*)
      READ(0,*) SC
      READ(0,*)
      READ(0,*) zip_size
      READ(0,*)
      READ(0,*) TAU
      
      CLOSE(0)

      RETURN
      END SUBROUTINE READ_INPUT
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     RANDOM ERDÖS-RÉNYI GRAPH WITH COUPLINGS GENERATOR
      SUBROUTINE IRS(N,M,p,NBR,INBR,JJ)
C     THIS SUBROUTINE GENERATES A RANDOM ERDÖS-RÉNYI GRAPH WITH p*M EDGES WITH
C     A WEIGHT OF 1 AND (1-p)*M EDGES WITH A WEIGHT OF 1
C     AND SAVES IT IN THE NBR, INBR AND JJ ARRAYS.

      INTEGER i, j, k
C     NODES, EDGES
      INTEGER N, M
C     FUNCTION TO GENERATE A U(0,1) RANDOM NUMBER
      REAL*8 genrand_real2 
C     FRACTION OF EDGES WITH VALUE 1
      REAL*8 p
      INTEGER edges_p, edges_n

      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
       
      edges_p = INT(p*M)
      edges_n = INT((1-p)*M)

      IF (edges_p+edges_n.NE.M) THEN
            PRINT*, 'ERROR: p VALUE NOT VALID'
            STOP
      END IF

      ALLOCATE(NBR(N))
      ALLOCATE(INBR(N))
      ALLOCATE(JJ(N))

      DO i = 1,N
            ALLOCATE(NBR(i)%v(0))
            ALLOCATE(JJ(i)%v(0))
      END DO

C     GENERATE M/2 EDGES OF WEIGHT 1
      k = 0
      DO WHILE(k<edges_p)
            i = INT(genrand_real2()*N) + 1 
            j = INT(genrand_real2()*N) + 1 
            IF ((ANY(NBR(i)%v == j).EQV..FALSE.).AND.(i.NE.j)) THEN
                  CALL ADDTOLIST(NBR(i)%v,j)
                  CALL ADDTOLIST(NBR(j)%v,i)
                  CALL ADDTOLIST(JJ(i)%v,1)
                  CALL ADDTOLIST(JJ(j)%v,1)
                  k = k + 1 
            END IF
      END DO
       
C     GENERATE M/2 EDGES OF WEIGHT -1 
      k = 0
      DO WHILE(k<edges_n)
            i = INT(genrand_real2()*N) + 1 
            j = INT(genrand_real2()*N) + 1 
            IF ((ANY(NBR(i)%v == j).EQV..FALSE.).AND.(i.NE.j)) THEN    
                  CALL ADDTOLIST(NBR(i)%v,j)
                  CALL ADDTOLIST(NBR(j)%v,i)
                  CALL ADDTOLIST(JJ(i)%v,-1)
                  CALL ADDTOLIST(JJ(j)%v,-1)
                  k = k + 1 
            END IF
      END DO

C	INBR GENERATION
      DO i = 1,N
            ALLOCATE(INBR(i)%v(SIZE(NBR(i)%v)))
      END DO
      DO i = 1,N 
            DO j = 1,SIZE(NBR(i)%v)
                  DO k = 1,SIZE(NBR(NBR(i)%v(j))%v)
                        IF ((NBR(NBR(i)%v(j))%v(k).EQ.i)) THEN
                              INBR(i)%v(j) = k
                        END IF
                  END DO
            END DO 
      END DO

      RETURN
      END SUBROUTINE IRS
C-----------------------------------------------------------------------

C------------------------------------------------------------------
C     RANDOM ERDÖS-RÉNYI GRAPH WITHOUT COUPLINGS GENERATOR
      SUBROUTINE IRG(N,M,NBR,INBR,JJ)
C     THIS SUBROUTINE GENERATES A RANDOM ERDÖS-RÉNYI GRAPH WITH M EDGES
C     AND SAVES IT IN THE NBR, INBR AND JJ ARRAYS.

      INTEGER i, j, k
C     NODES, EDGES
      INTEGER N, M
C     FUNCTION TO GENERATE A U(0,1) RANDOM NUMBER
      REAL*8 genrand_real2 

      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      ALLOCATE(NBR(N))
      ALLOCATE(INBR(N))
      ALLOCATE(JJ(N))

      DO i=1,N
            ALLOCATE(NBR(i)%v(0))
            ALLOCATE(JJ(i)%v(0))
      END DO

C     GENERATE M EDGES OF WEIGHT 0
      k = 0
      DO WHILE(k<M)
            i = INT(genrand_real2()*N) + 1 
            j = INT(genrand_real2()*N) + 1 
            IF ((ANY(NBR(i)%v == j).EQV..FALSE.).AND.(i.NE.j)) THEN
                  CALL ADDTOLIST(NBR(i)%v,j)
                  CALL ADDTOLIST(NBR(j)%v,i)
                  CALL ADDTOLIST(JJ(i)%v,0)
                  CALL ADDTOLIST(JJ(j)%v,0)
                  k = k + 1 
            END IF
      END DO

C     INBR GENERATION
      DO i = 1,N
            allocate(INBR(i)%v(size(NBR(i)%v)))
      END DO
      DO i = 1,N 
            DO j = 1,SIZE(NBR(i)%v)
                  DO k = 1,SIZE(NBR(NBR(i)%v(j))%v)
                        IF ((NBR(NBR(i)%v(j))%v(k).EQ.i)) THEN
                              INBR(i)%v(j) = k
                        END IF
                  END DO
            END DO
      END DO

      RETURN
      END SUBROUTINE IRG
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     RANDOM COUPLING ASSIGMENT
      SUBROUTINE RCA(N,M,p,NBR,INBR,JJ)
C     THIS SUBROUTINE ASSIGNS TO THE GRAPH P*M EDGES WITH
C     A WEIGHT OF 1 AND (1-P)*M EDGES WITH A WEIGHT OF -1
C     UPDATING THE ARRAY JJ.

      INTEGER i, j, k
      INTEGER N, M
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
C     FRACTION OF EDGES WITH VALUE 1
      REAL*8 p
      INTEGER edges_p, edges_n
C     FUNCTION TO GENERATE A U(0,1) RANDOM NUMBER
      REAL*8 genrand_real2 

      edges_p = INT(p*M)
      edges_n = INT((1-p)*M)

      IF (edges_p+edges_n.NE.M) THEN
            PRINT*, 'PROBLEM: p VALUE NOT VALID'
            STOP
      END IF

C     GENERATE M*p EDGES OF WEIGHT 1
      k = 0
      DO WHILE(k<edges_p)
1           i = INT(genrand_real2()*N) + 1 
            j = INT(genrand_real2()*SIZE(JJ(i)%v)) + 1 
            IF (SIZE(NBR(i)%v).EQ.0) GO TO 1
            IF (JJ(i)%v(j).EQ.0) THEN
                  JJ(i)%v(j) = -1
                  JJ(NBR(i)%v(j))%v(INBR(i)%v(j)) = -1
                  k = k + 1 
            END IF
      END DO

C     GENERATE M*(1-p) EDGES OF WEIGHT -1
      k = 0
      DO WHILE(k<edges_n)
2           i = INT(genrand_real2()*N) + 1 
            j = INT(genrand_real2()*SIZE(JJ(i)%v)) + 1 
            IF (SIZE(NBR(i)%v).EQ.0) GO TO 2
            IF (JJ(i)%v(j).EQ.0) THEN
                  JJ(i)%v(j) = 1
                  JJ(NBR(i)%v(j))%v(INBR(i)%v(j)) = 1
                  k = k + 1 
            END IF
      END DO

      RETURN
      END SUBROUTINE RCA
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      REAL*8 FUNCTION ENERG(S,N,NBR,JJ)
C     THIS FUNCTION CALCULATES THE ENERGY OF THE SYSTEM GIVEN AN
C     S CONFIGURATION

      INTEGER S(1:N)
      INTEGER i, N
      REAL*8 ENE

      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      ENE = 0.0d0
      
      DO i = 1,N
            DO k= 1,SIZE(NBR(i)%v)
                  ENE = ENE - S(i)*S(NBR(i)%v(k))*JJ(i)%v(k)
            END DO
      END DO

      ENERG =  ENE/2
      
      RETURN
      END FUNCTION ENERG
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     METROPOLIS ALGORITHM 
      SUBROUTINE METROPOLIS(S,N,valid,TEMP,DE,NBR,JJ)
C     THIS SUBROUTINE PROPOSES A CHANGE OF SPIN IN A RANDOM NODE,
C     CALCULATES THE ENERGY VARIATION OF THE SYSTEM (ΔH) DUE TO IT,
C     IF ΔH < 0 THEN THE CHANGE IS ACCPETED, ELSE IF ΔH > 0 THEN THE
C     CHANGE IS ACCEPTED WITH A PROBABILITY OF EXP(-ΔH/k_BT).

      INTEGER i, k, N
      INTEGER S(1:N)
      REAL*8 genrand_real2
      LOGICAL valid

      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      REAL*8 suma, TEMP, DE

      valid = .FALSE.

C     RANDOM NODE SELECTION
      i = INT(genrand_real2()*N) + 1

C     CALCULATION OF ΔH
      suma = 0
      DO k=1,SIZE(NBR(i)%v)
            suma = suma + JJ(i)%v(k)*S(NBR(i)%v(k))
      END DO
      DE = 2*S(i)*suma

C     CHECK IF THE CHANGE IS ACCEPTED
      IF (genrand_real2().LT.min(1.d0,exp(-DE/TEMP))) THEN
            S(i) = -S(i)
            valid = .true.
      END IF

      RETURN
      END SUBROUTINE METROPOLIS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     ARRAY TO BINARY
      SUBROUTINE ARRAY2BIN(N,binary,array)
C     THIS SUBRUTINE CONVERTS AN ARRAY OF -1 AND 1 TO A BINARY NUMBER

      INTEGER i, N, array(1:N)
      CHARACTER(N) binary

      DO i = 1,N
            binary(i:i) = '0'
      END DO
      DO i = 1,N
            IF (array(i) == 1) THEN
                  binary(i:i) = '1'
            END IF
      END DO

      RETURN
      END SUBROUTINE ARRAY2BIN
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     BINARY TO DECIMALS
      SUBROUTINE BIN2DEC(N,zip_size,binary,decimal)
C     THIS SUBRUTINE CONVERTS A BINARY NUMBER TO N/zip_size DECIMAL NUMBERS

      INTEGER i, j, N, zip_size, scale
      CHARACTER(N) binary
      INTEGER decimal(1:N/zip_size)

      DO j = 1,N/zip_size
            decimal(j) = 0
            scale = (j-1)*zip_size
            DO i = 1,zip_size
                  IF (binary(scale+i:scale+i).EQ.'1') THEN
                        decimal(j) = decimal(j) + 2**(zip_size-i)
                  END IF
            END DO
      END DO

      RETURN
      END SUBROUTINE BIN2DEC
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     DECIMALS TO BINARY
      SUBROUTINE DEC2BIN(N,zip_size,binary,decimal)
C     THIS SUBRUTINE CONVERTS N/zip_size DECIMAL NUMBERS TO A BINARY NUMBER

      INTEGER i, j, N, zip_size, scale
      CHARACTER(N) binary
      INTEGER decimal(1:N/zip_size)
      INTEGER decimal_copy(1:N/zip_size)

      decimal_copy = decimal
      DO j=1,N/zip_size
      scale = (j-1)*zip_size
      DO i = zip_size,1,-1
            IF (MOD(decimal_copy(j),2)==1) THEN
                  binary(scale+i:scale+i) = '1'
            ELSE
                  binary(scale+i:scale+i) = '0'
            END IF
            decimal_copy(j) = decimal_copy(j)/2
      END DO
      END DO

      RETURN
      END SUBROUTINE DEC2BIN
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     BINARY TO ARRAY
      SUBROUTINE BIN2ARRAY(N,binary,array)
C     THIS SUBRUTINE CONVERTS A BINARY NUMBER TO AN ARRAY OF -1 AND 1

      INTEGER i, N, array(1:N)
      CHARACTER(N) binary

      DO i = 1,N
            IF (binary(i:i).EQ.'1') THEN
                  array(i) = 1
            ELSE
                  array(i) = -1
            END IF
      END DO

      RETURN
      END SUBROUTINE BIN2ARRAY
C-----------------------------------------------------------------------

C------------------------------------------------------------------
C     ADD ELEMENT TO LIST
      SUBROUTINE ADDTOLIST(list, element)

      INTEGER i, isize
      INTEGER element
      INTEGER, DIMENSION(:), ALLOCATABLE:: list
      INTEGER, DIMENSION(:), ALLOCATABLE :: clist

      IF (ALLOCATED(list)) THEN
            isize = size(list)
            ALLOCATE(clist(isize+1))
            DO i = 1,isize          
                  clist(i) = list(i)
            END DO
            clist(isize+1) = element
            DEALLOCATE(list)
            CALL MOVE_ALLOC(clist, list)
      ELSE
            ALLOCATE(list(1))
            list(1) = element
      END IF

      RETURN
      END SUBROUTINE ADDTOLIST
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     REMOVE INDEX FROM LIST
      SUBROUTINE RMVOFLIST(list, index)

      INTEGER i, isize
      INTEGER index
      INTEGER, DIMENSION(:), ALLOCATABLE:: list
      INTEGER, DIMENSION(:), ALLOCATABLE:: clist

      IF (ALLOCATED(list)) THEN
            isize = SIZE(list)
            ALLOCATE(clist(isize-1))
            DO i = 1,index-1
                  clist(i) = list(i)
            END DO
            DO i = index,isize-1
                  clist(i) = list(i+1)
            END DO
            DEALLOCATE(list)
            CALL MOVE_ALLOC(clist, list)

      END IF

      RETURN
      END SUBROUTINE RMVOFLIST
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     PERFORMANCE.F
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     TABULATE LAMBDA FOR A GIVEN S_SET        
      SUBROUTINE CLASS_LAMBDA(N,C,S_SET,NBR,JJ,LAMBDA)
C     THIS SUBROUTINE TABULATES THE VALUES OF THE SUM INSIDE THE
C     TANH OF THE PSEUDOLIKELIHOOD FOR A GIVEN S_SET AND GRAPH
C     DEFINED BY NBR,INBR AND JJ.

      INTEGER i, N, k, m
      INTEGER C
      INTEGER sum
      INTEGER S_SET(1:C,1:N), LAMBDA(1:C,1:N)
      TYPE(multi_array),ALLOCATABLE:: NBR(:)
      TYPE(multi_array),ALLOCATABLE:: JJ(:)

      sum = 0
      DO m = 1,C
            DO i = 1,N
                  DO k = 1, SIZE(NBR(i)%v)
                        sum = sum + JJ(i)%v(k)*S_SET(m,NBR(i)%v(k))
                  ENDDO
                  LAMBDA(m,i) = sum
                  sum = 0
            ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CLASS_LAMBDA
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     TABULATE FUNCT FOR A GIVEN TEMP           
      SUBROUTINE CLASS_FUNCTION(zmax,TEMP,funct)
C     THIS SUBROUTINE TABULATES THE FUNCTION IN THE SUMATORIES OF
C     THE PSEUDOLIKELIHOOD.

      INTEGER i, j, zmax
      REAL*8 TEMP
      REAL*8 funct(-1:1,-zmax:zmax)

      DO i = -1,1,2
            DO j = -zmax,zmax
                  funct(i,j) = LOG(0.5d0*(1 + i*TANH(j/TEMP)))
            ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CLASS_FUNCTION
C-----------------------------------------------------------------------

c------------------------------------------------------------------
C     FIND LOCATION OF A VALUE IN AN ARRAY
      INTEGER FUNCTION FLOC(array,value)

      INTEGER,ALLOCATABLE:: array(:)
      INTEGER value, r
      LOGICAL out

      floc = 0
      DO r = 1, size(array)
            IF (array(r).EQ.value) THEN
                  floc = r
                  out = .TRUE.
            END IF
      END DO

      RETURN
      END FUNCTION FLOC
C-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     PSEUDOLIKELIHOOD ALGORITHM 
c-----------------------------------------------------------------------
      SUBROUTINE PSEUDOLIKELIHOOD(N,C,S_SET,
     . valid,TEMP_F,DPL,NBR,INBR,JJ,funct,zmax,LAMBDA)
C     THIS SUBROUTINE PROPOSES A CHANGE OF PAIRWISE COUPLING,
C     CALCULATES THE PSEUDOLIKELIHOOD VARIATION OF THE SYSTEM (ΔPL) DUE TO IT,
C     IF ΔPL > 0 THEN THE CHANGE IS ACCPETED, ELSE IF ΔPL < 0 THEN THE
C     CHANGE IS ACCEPTED WITH A PROBABILITY OF EXP(-ΔPL/k_BT'), WHERE
C     T' IS THE FICTICIOUS TEMPERATURE

      INTEGER N
      LOGICAL valid
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      INTEGER C, S_SET(1:C,1:N)
      REAL*8 TEMP_F
      LOGICAL change
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newNBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newINBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newJJ(:)
C     FUNCTION TO GENERATE A U(0,1) RANDOM NUMBER
      REAL*8 genrand_real2 
      INTEGER zmax, m
      
      REAL*8 sum
      INTEGER i1,i2,i3,i4
      REAL*8  DPL
      INTEGER ip,ip2

      REAL*8 funct(-1:1,-zmax:zmax)
      INTEGER LAMBDA(1:C,1:N), newLAMBDA(1:C,1:N)
      
      valid = .FALSE.
      change = .FALSE. 

      ALLOCATE(newNBR(N))
      ALLOCATE(newINBR(N))
      ALLOCATE(newJJ(N))

      DO i = 1,N
            ALLOCATE(newNBR(i)%v(SIZE(NBR(i)%v)))
            ALLOCATE(newINBR(i)%v(SIZE(NBR(i)%v)))
            ALLOCATE(newJJ(i)%v(SIZE(NBR(i)%v)))
      END DO

      newNBR = NBR
      newINBR = INBR
      newJJ = JJ 

      newLAMBDA = LAMBDA

C     RANDOM PAIRWISE COUPLING CHANGE
      DO WHILE (change.EQV..FALSE.)
      CALL JJ_CHANGE(N,NBR,INBR,JJ,newNBR,newINBR,newJJ,change,
     .              i1,i2,i3,i4)
      END DO

C     CALCULATE newLAMBDA
      DO m = 1,C

      ip = floc(NBR(i1)%v,i3)
      ip2 = floc(newNBR(i1)%v,i3)
      IF (ip.EQ.0) THEN
      newLAMBDA(m,i1) = newLAMBDA(m,i1)+newJJ(i1)%v(ip2)*S_SET(m,i3)
      ELSE IF (ip.NE.0) THEN
      newLAMBDA(m,i1) = newLAMBDA(m,i1)+2*newJJ(i1)%v(ip2)*S_SET(m,i3)
      END IF

      ip = floc(NBR(i3)%v,i1)
      ip2 = floc(newNBR(i3)%v,i1)
      IF (ip.EQ.0) THEN
      newLAMBDA(m,i3) = newLAMBDA(m,i3)+newJJ(i3)%v(ip2)*S_SET(m,i1)
      ELSE IF (ip.NE.0) THEN
      newLAMBDA(m,i3) = newLAMBDA(m,i3)+2*newJJ(i3)%v(ip2)*S_SET(m,i1)
      END IF

      ip = floc(NBR(i2)%v,i4)
      ip2 = floc(newNBR(i2)%v,i4)
      IF (ip2.EQ.0) THEN
      newLAMBDA(m,i2) = newLAMBDA(m,i2)-JJ(i2)%v(ip)*S_SET(m,i4)
      ELSE IF (ip2.NE.0) THEN
      newLAMBDA(m,i2) = newLAMBDA(m,i2)+2*newJJ(i2)%v(ip2)*S_SET(m,i4)
      END IF 

      ip = floc(NBR(i4)%v,i2)
      ip2 = floc(newNBR(i4)%v,i2)
      IF (ip2.EQ.0) THEN
      newLAMBDA(m,i4) = newLAMBDA(m,i4)-JJ(i4)%v(ip)*S_SET(m,i2)
      ELSE IF (ip2.NE.0) THEN
      newLAMBDA(m,i4) = newLAMBDA(m,i4)+2*newJJ(i4)%v(ip2)*S_SET(m,i2)
      END IF 
      
      END DO      

C     CALCULATE DPL
      sum = 0

      IF ((i1.NE.i4).AND.(i2.NE.i3).AND.(i3.NE.i4)) THEN
      DO m = 1,C
            sum = sum + 
     .      funct(S_SET(m,i1),newLAMBDA(m,i1))
     .    - funct(S_SET(m,i1),LAMBDA(m,i1)) +
     .      funct(S_SET(m,i2),newLAMBDA(m,i2))
     .     -funct(S_SET(m,i2),LAMBDA(m,i2)) + 
     .      funct(S_SET(m,i3),newLAMBDA(m,i3)) 
     .     -funct(S_SET(m,i3),LAMBDA(m,i3)) +
     .      funct(S_SET(m,i4),newLAMBDA(m,i4))
     .     -funct(S_SET(m,i4),LAMBDA(m,i4))
      END DO
      DPL = sum/C
      END IF 

      IF (i1.EQ.i4) THEN
      do m = 1,C
            sum = sum + 
     .      funct(S_SET(m,i1),newLAMBDA(m,i1))
     .    - funct(S_SET(m,i1),LAMBDA(m,i1)) +
     .      funct(S_SET(m,i2),newLAMBDA(m,i2))
     .     -funct(S_SET(m,i2),LAMBDA(m,i2)) + 
     .      funct(S_SET(m,i3),newLAMBDA(m,i3)) 
     .     -funct(S_SET(m,i3),LAMBDA(m,i3))
      END DO
      DPL = sum/C
      END IF

      IF ((i2.EQ.i3).OR.(i4.EQ.i3)) THEN
      DO m = 1,C
            sum = sum + 
     .      funct(S_SET(m,i1),newLAMBDA(m,i1))
     .    - funct(S_SET(m,i1),LAMBDA(m,i1)) +
     .      funct(S_SET(m,i2),newLAMBDA(m,i2))
     .     -funct(S_SET(m,i2),LAMBDA(m,i2)) + 
     .      funct(S_SET(m,i4),newLAMBDA(m,i4))
     .     -funct(S_SET(m,i4),LAMBDA(m,i4))
      END DO
      DPL = sum/C
      END IF

C     CHECK IF THE CHANGE IS ACCEPTED
      IF (genrand_real2().LT.MIN(1.d0,EXP(DPL/TEMP_F))) THEN
            valid = .TRUE.
            NBR = newNBR
            INBR = newINBR
            JJ = newJJ
            
            LAMBDA = newLAMBDA
      END IF

      RETURN
      END SUBROUTINE PSEUDOLIKELIHOOD

C------------------------------------------------------------------------------
      REAL*8 FUNCTION PSEUDO(N,C,S_SET,TEMP,NBR,JJ)
C     THIS FUNCTION CALCULATES THE PSEUDOLIKELIHOOD FOR A GIVEN TEMP,
C     S_SET AND GRAPH DEFINED BY NBR AND JJ
      IMPLICIT NONE
      INTEGER N, C, i, k, m
      REAL*8 TEMP
      REAL*8 L
      REAL*8 sum

      INTEGER S_SET(1:C,1:N)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)

      PSEUDO = 0.0d0
      L = 0.0d0
      sum = 0.0d0

      DO i = 1,N
      DO m = 1,C
            DO k = 1, SIZE(JJ(i)%V)
                  sum = sum + JJ(i)%V(k)*S_SET(m,NBR(i)%V(k))
            ENDDO
            
            L = L + LOG(0.5d0*(1 + S_SET(m,i)*EXP(TEMP*sum)))
            sum = 0.0d0
      ENDDO
      PSEUDO = PSEUDO + L
      L = 0.0d0
      ENDDO

      RETURN
      END FUNCTION PSEUDO
C------------------------------------------------------------------------------

C------------------------------------------------------------------------------
c     MAKE A JJ EXCHANGE OF TYPE 1 OR 2 (WHEN THE GRAPH STRUCTURE IS NOT FIXED)
      SUBROUTINE JJ_CHANGE(N,NBR,INBR,JJ,newNBR,newINBR,newJJ,change,
     .                     i1,i2,i3,i4)

C     THIS SUBROUTINE CHANGES PAIRWISE COUPLING OF i1-i3 FOR THE
C     PAIRWISE COUPLING OF i2-i4, AND GENERATES A THE NEW NBR, INBR
C     AND JJ ARRAYS.

      IMPLICIT NONE
      INTEGER N
      INTEGER i, j, k
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newNBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newINBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newJJ(:)
C     FUNCTION TO GENERATE A U(0,1) RANDOM NUMBER
      REAL*8 GENRAND_REAL2 
      INTEGER i1,k1, i2,k2, i3,k3, i4,k4

      INTEGER ip
      
      LOGICAL CHANGE, VALID, VALID2, VALID3

      CHANGE = .FALSE.
      VALID = .FALSE.
      VALID2 = .FALSE.
      VALID3 = .FALSE.

C     SELECTING A RANDOM NODE
      i1 = INT(genrand_real2()*N) + 1 

C     SELECTING A RANDOM EDGE OF THAT NODE, THE EXTRA +1
C     IS BECAUSE WE CAN SELECT A NODE WITH VALUE 0
      k1 = int(genrand_real2()*size(NBR(i1)%v)) + 1 + 1 

C     IN CASE SIZE(NBR(i,k))=0 K SHOULD BE EQUAL TO 1
      IF (SIZE(NBR(i1)%V) == 0) THEN
            k1 = 1
      ENDIF

C     IF THE MAX EDGES PER NODE IS EXCEDED WE DENY THE CHANGE
      IF (k1 > N-1) THEN
            CHANGE = .FALSE.
            RETURN
      ENDIF

C-------------------------------------------------------------
C     IF JJ(i1,k1) = 0 WE WANT A 0 <----> 1,-1 EXCHANGE
C-------------------------------------------------------------
      IF (k1 > SIZE(NBR(i1)%v)) THEN

C     DIFERENT RANDOM NODE SELECTION i2
      DO WHILE (valid.EQV..FALSE.)
      i2 = INT(genrand_real2()*N) + 1
      IF ((i1.NE.i2).AND.(SIZE(NBR(i2)%v).GT.1)) THEN 
            valid =  .true.
      END IF
      END DO

C     i2 CAN'T BE THE SAME AS i1
C     TO AVOID PROBLEMS I2 NEEDS TO HAVE MORE THAN 1 NEIGHBOR
C     OTHERWISE A NODE COULD END UP WITH 0 EDGES

C     SELECTING A RANDOM NEIGHBOR OF THE NODE i2
      k2 = INT(genrand_real2()*SIZE(NBR(i2)%v)) + 1

C     SELECTING A NODE i3 THAT IS NOT A NEIGBOR  OF i1 ALREADY
      DO WHILE (valid2.EQV..FALSE.)
      i3 = INT(genrand_real2()*N) + 1
C     THIS ALLOW TO SEE IF I3 IS ALREADY A NEIGBOUR OF I1
      ip = floc(NBR(i1)%v,i3)
C     ip IS 0 IF i3 IS NOT IN NBR(i1)%v
      IF ((ip.EQ.0).AND.(i1.NE.i3).AND.(SIZE(NBR(i3)%v).LT.N-1)) THEN
            valid2 = .TRUE.
      END IF
      END DO 

      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

C     ADDING THE NEW EDGE
      CALL ADDTOLIST(newNBR(i1)%v,i3)
      CALL ADDTOLIST(newNBR(i3)%v,i1)
      CALL ADDTOLIST(newJJ(i1)%v,JJ(i2)%v(k2))
      CALL ADDTOLIST(newJJ(i3)%v,JJ(i4)%v(k4))

C     REMOVING THE OLD EDGE
      CALL RMVOFLIST(newNBR(i2)%v,k2)
      CALL RMVOFLIST(newNBR(i4)%v,k4)
      CALL RMVOFLIST(newJJ(i2)%v,k2)
      CALL RMVOFLIST(newJJ(i4)%v,k4)

C     UPDATING THE INBR     
      newINBR = newNBR
      DO i = 1,N 
            DO j = 1,SIZE(newNBR(i)%v)
                  DO k = 1,SIZE(newNBR(newNBR(i)%v(j))%v)
                        IF ((newNBR(newNBR(i)%v(j))%v(k).EQ.i)) THEN
                              newINBR(i)%v(j) = k
                        END IF
                  END DO
            END DO
      END DO

      change = .true.
      RETURN
      END IF

C-------------------------------------------------------------
C     IF JJ(i1,k1) = +1,-1 WE WANT a +-1 <----> -+1 EXCHANGE
C-------------------------------------------------------------
      IF (k1 <= SIZE(NBR(i1)%v)) THEN
      
C     DIFERENT RANDOM NODE SELECTION i2
      DO WHILE (valid.EQV..false.)
      i2 = int(genrand_real2()*N) + 1

      IF ((i1.NE.i2).AND.(size(NBR(i2)%v).gt.1)) THEN
            valid = .TRUE.
      END IF
      END DO

C     i2 CAN'T BE THE SAME AS i1
C     SELECTING A RANDOM NEIGHBOR OF THE NODE i2
      DO WHILE (valid2.EQV..false.)
      k2 = INT(genrand_real2()*SIZE(NBR(i2)%v)) + 1
      IF (NBR(i2)%v(k2).NE.i1) THEN
            valid2 = .TRUE.
      END IF
      END DO
C     IF BOTH EDGES HAVE THE SAME VALUE THERE IS NO CHANGE
      IF (JJ(i1)%v(k1).EQ.JJ(i2)%v(k2)) THEN
            change = .false.
            RETURN
      END IF
      
C     CHANGING THE EDGE VALUES
      i3 = NBR(i1)%v(k1)
      k3 = INBR(i1)%v(k1)
      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

      newJJ(i1)%v(k1) = JJ(i2)%v(k2)
      newJJ(i3)%v(k3) = JJ(i2)%v(k2)
      newJJ(i2)%v(k2) = JJ(i1)%v(k1)
      newJJ(i4)%v(k4) = JJ(i1)%v(k1)

      change = .TRUE.
      RETURN
      END IF 

      END SUBROUTINE JJ_CHANGE
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     MAKE A JJ EXCHANGE OF TYPE 2 (WHEN THE GRAPH STRUCTURE IS FIXED)
      SUBROUTINE JJ_CHANGE2(N,NBR,INBR,JJ,newNBR,newINBR,newJJ,change,
     .                     i1,i2,i3,i4)
C     THIS SUBROUTINE CHANGES PAIRWISE COUPLING OF i1-i3 FOR THE
C     PAIRWISE COUPLING OF i2-i4, AND GENERATES A THE NEW NBR, INBR
C     AND JJ ARRAYS.

      INTEGER N
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newNBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newINBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: newJJ(:)
C     FUNCTION TO GENERATE A U(0,1) RANDOM NUMBER
      REAL*8 genrand_real2 
      INTEGER i1,k1, i2,k2, i3,k3, i4,k4
      
      LOGICAL change, valid, valid2, valid3

      change = .FALSE.
      valid = .FALSE.
      valid2 = .FALSE.
      valid3 = .FALSE.

C     SELECTING A RANDOM NODE
      i1 = INT(genrand_real2()*N) + 1 

C     SELECTING A RANDOM EDGE OF THAT NODE
      k1 = INT(genrand_real2()*SIZE(NBR(i1)%v)) + 1

C     IN CASE SIZE(NBR(i,k))= 0 NO CHANGE POSSIBLE
      IF (SIZE(NBR(i1)%v).EQ.0) THEN
            change = .FALSE.
            RETURN
      END IF

C-------------------------------------------------------------
C     IF JJ(i1,k1) = -1,1 WE WANT A +-1 <----> -+1 EXCHANGE
C-------------------------------------------------------------

      newNBR = NBR
      newINBR = INBR

      IF (k1 <= SIZE(NBR(i1)%v)) THEN
      
C     DIFERENT RANDOM NODE SELECTION i2
      DO WHILE (valid.EQV..false.)
      i2 = int(genrand_real2()*N) + 1

      IF ((i1.NE.i2).AND.(size(NBR(i2)%v).gt.1)) THEN
            valid = .TRUE.
      END IF
      END DO

C     i2 CAN'T BE THE SAME AS i1
C     SELECTING A RANDOM NEIGHBOR OF THE NODE i2
      DO WHILE (valid2.EQV..false.)
      k2 = INT(genrand_real2()*SIZE(NBR(i2)%v)) + 1
      IF (NBR(i2)%v(k2).NE.i1) THEN
            valid2 = .TRUE.
      END IF
      END DO
C     IF BOTH EDGES HAVE THE SAME VALUE THERE IS NO CHANGE
      IF (JJ(i1)%v(k1).EQ.JJ(i2)%v(k2)) THEN
            change = .false.
            RETURN
      END IF
      
C     CHANGING THE EDGE VALUES
      i3 = NBR(i1)%v(k1)
      k3 = INBR(i1)%v(k1)
      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

      newJJ(i1)%v(k1) = JJ(i2)%v(k2)
      newJJ(i3)%v(k3) = JJ(i2)%v(k2)
      newJJ(i2)%v(k2) = JJ(i1)%v(k1)
      newJJ(i4)%v(k4) = JJ(i1)%v(k1)

      change = .TRUE.
      RETURN
      END IF 

      END SUBROUTINE JJ_CHANGE2
C-----------------------------------------------------------------------

C------------------------------------------------------------------
      REAL*8 FUNCTION GAMMAA(N,M,NBR,JJ,NBR_0,JJ_0)
C     THIS FUNCTION CALCULATES THE GAMMA BETWEEN JJ AND JJ_0

      INTEGER N, M
      REAL*8 sum
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR_0(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ_0(:)
      INTEGER ip

      sum = 0.d0
      DO i = 1,N
            DO k = 1,SIZE(NBR_0(i)%v)
                  IF (i.LT.NBR_0(i)%v(k)) THEN
                  ip = floc(NBR(i)%v,NBR_0(i)%v(k))
                  IF (ip.EQ.0) THEN
                        sum = sum + (JJ_0(i)%v(k))**2
                  END IF
                  IF (ip.ne.0) THEN
                        sum = sum + (JJ_0(i)%v(k)-JJ(i)%v(ip))**2
                  END IF
                  END IF
            END DO
      END DO

      GAMMAA = sum/M
      RETURN
      END FUNCTION GAMMAA
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     SAVE FINAL COUPLINGS
      SUBROUTINE SAVE_COUPLINGS(N,str1,str2,NBR,JJ)
C     THIS SUBROTINE SAVES THE FINAL COUPLINGS IN J_str1_str2.dat

      INTEGER i,N
      CHARACTER(3) str1,str2
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      INTEGER LENG(1:N)

      DO i = 1,N
            LENG(i)= SIZE(NBR(i)%V)
      ENDDO
      OPEN(UNIT=2,FILE='J_'//str1//'_'//str2//'.dat')
      WRITE(2,*) LENG
      DO i = 1,N
            WRITE(2,*) NBR(i)%V
      ENDDO
      DO i = 1,N
            WRITE(2,*) JJ(i)%V
      ENDDO
      CLOSE(2)

      END SUBROUTINE SAVE_COUPLINGS
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     READ FINAL COUPLINGS
      SUBROUTINE READ_COUPLINGS(N,str1,str2,NBR,JJ)
C     THIS SUBROTINE READS THE FINAL COUPLINGS IN J_str1_str2.dat

      INTEGER i,N
      CHARACTER(3) str1,str2
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
      INTEGER LENG(1:N)

      OPEN(UNIT=2,FILE='J_'//str1//'_'//str2//'.dat')
      READ(2,*) LENG
      DO i = 1,N
            ALLOCATE(NBR(i)%V(LENG(i)))
            ALLOCATE(JJ(i)%V(LENG(i)))
      ENDDO
      DO i = 1,N
            READ(2,*) NBR(i)%V
      ENDDO
      DO i = 1,N
            READ(2,*) JJ(i)%V
      ENDDO
      CLOSE(2)

      END SUBROUTINE READ_COUPLINGS
C-----------------------------------------------------------------------

      END MODULE SPIN_MODEL