c     EQUILIBRIUM RECONSTRUCTION
c     LLUÍS TORRES HUGAS

      implicit none
c------------------------------------------------------------------
c------------------------------------------------------------------
c     nodes, edges, connectivity
      integer N, M,z
      parameter (N = 100)
      parameter (z = 4)
      parameter (M = z*N/2)
      integer MCTOT, MCINI, SC
c     total MonteCarlo steps (MCS)
      parameter (MCTOT =  30 000)
c     MCS till we consider equilibrium
      parameter (MCINI = MCTOT/3)
c     save spin configuration every SC (MCS)
      parameter (SC = 100)
c     temperature (TEMP = k_B·T)
      real*8 TEMP
      parameter(TEMP = 1.8d0) 
c     +1 -1 edges ratio (1 => all +1), (0 => all -1)
      real*8 p
c------------------------------------------------------------------
c------------------------------------------------------------------
c     number of graphs in simulated for every p
      integer NSEEDS
      parameter(NSEEDS = 100)
c------------------------------------------------------------------
c------------------------------------------------------------------
      integer C
c     sample size (number of spin config)
      parameter(C = 2*(MCTOT-MCINI)/SC)
c     maximum number of neighbours that a node can have
      integer zmax
      parameter (zmax = N-1) 
      integer MCTOT2, MCINI2, MCD2
c     montecarlo total steps
      parameter (MCTOT2 =  2 000)
      parameter (MCINI2 = 1)
c     rate of montecarlo averages
      parameter (MCD2 = 1)
c------------------------------------------------------------------
c------------------------------------------------------------------
      integer i, j, k, IMC, IPAS
c     seed number
      integer SEED, SEEDini
c     multidimensional array type
      type multi_array
      integer,allocatable::v(:)
      end type multi_array
c     principal arrays
      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: INBR(:)
      type(multi_array),allocatable:: JJ(:)

      type(multi_array),allocatable:: NBR_0(:)
      type(multi_array),allocatable:: INBR_0(:)
      type(multi_array),allocatable:: JJ_0(:)

      integer, allocatable:: S_SET(:,:)
      
      logical valid, change
      real*8 TF_STEP
      real*8 PL, DPL, SUMPL, SUMg
      real*4 TIME1, TIME2, time
c     fictitious temperature (TEMP_F = k_B·TF)
      real*8 TEMP_F

c------------------------------------------------------------------
c     number of p values to sumulate
      integer npvalues, np
      real*8 pvalues(1:7)
c------------------------------------------------------------------

      character(4) str5
      character(3) str6

      integer LAMBDA(1:C,1:N)
      real*8 funct(-1:1,-zmax:zmax)

c------------------------------------------------------------------
c     START
c------------------------------------------------------------------

      SEEDini = 100
      
c     all p values
      pvalues = (/0.d0,0.15d0,0.3d0,0.5d0,0.7d0,0.85d0,1.d0/)
      npvalues = size(pvalues)

      call cpu_time(TIME1)

c¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
c     do for all p values
      DO np=1,npvalues

      p = pvalues(np)
      
      print*,' '
      print*,'p# =', p
      print*,' '

      allocate(S_SET(1:C,1:N))

      write( str5,'(f4.2)')   p 

c     gamma(var(q)) data file
      open(unit = 10,file = 'g_value_'//str5//'.dat')

c¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
c     do for all p values
      DO SEED = SEEDini,SEEDini+NSEEDS-1

      write( str6,'(i3)')     SEED

c..................................................................
c	Generate the Original Graph
      print*,'SEED# =', SEED
      print*,' '

      call init_genrand(SEED) 
      call IRG(N,M,p,NBR_0,INBR_0,JJ_0)
      print*,'Original Graph Generated'
      print*,' '
c..................................................................
c     Read the set of spin configurations
      open(unit = 1,file = 'spin_config_'//str5//'_'//str6//'.dat')
      do i = 1, C
            read(1,*) S_SET(i,:)
      enddo
      close(1)

c..................................................................
      call init_genrand(999)
c..................................................................
c     ficticious tempereture 
      TEMP_F = -log(0.5d0*(1+tanh(1./TEMP)))/3
      TF_STEP = TEMP_F/(MCTOT2)
      
c..................................................................
c     Initial NBR, INBR and JJ random arrays generation
      call IRG(N,M,p,NBR,INBR,JJ)

c..................................................................
      call class_function(zmax,TEMP,funct)
      call class_lambda(N,C,S_SET,NBR,JJ,LAMBDA)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     MC SIMULATION
      PL = PEUDOL(N,C,S_SET,TEMP,NBR,JJ)

      SUMPL = 0
      SUMg = 0

      print*,'MC Simulation Started'
      print*,' '
c     MC simulation
      DO IMC = 1, MCTOT2

            do IPAS = 1, M
            
            call pseudolikelihood(N,C,S_SET,valid,TEMP_F,
     .                        DPL,NBR,INBR,JJ,funct,zmax,LAMBDA)
            if (valid) then
                  PL = PL + DPL
            endif
            enddo

            SUMPL = SUMPL + PL

            SUMg = SUMg + g(N,M,NBR,JJ,NBR_0,JJ_0)

c           average values every MCD2 montecarlo steps
            if ((IMC.gt.MCINI2).and.(MCD2*(IMC/MCD2).eq.IMC)) then
                  
                  SUMPL = 0
                  SUMg = 0

            end if

            TEMP_F = TEMP_F - TF_STEP

      ENDDO

      write (10,*) p, SEED, g(N,M,NBR,JJ,NBR_0,JJ_0)
      print*,'γ =',g(N,M,NBR,JJ,NBR_0,JJ_0)
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do i = 1,N
            DEALLOCATE(NBR(i)%v)
            DEALLOCATE(INBR(i)%v)
            DEALLOCATE(JJ(i)%v)
      enddo
      DEALLOCATE(NBR)
      DEALLOCATE(INBR)
      DEALLOCATE(JJ)

c------------------------------------------------------------------

      do i = 1,N
            DEALLOCATE(NBR_0(i)%v)
            DEALLOCATE(INBR_0(i)%v)
            DEALLOCATE(JJ_0(i)%v)
      enddo
      DEALLOCATE(NBR_0)
      DEALLOCATE(INBR_0)
      DEALLOCATE(JJ_0)


      if ((SEED.eq.SEEDini).and.(p.eq.pvalues(1)))then
      call cpu_time(TIME2)
      time = (TIME2-TIME1)*NSEEDS*npvalues
      print*, "ESTIMATED TIME: ", int(time/3600), 'h',
     . int((time/3600-int(time/3600))*60), 'min', 
     . int((time/60-int(time/60))*60), 's'
      print*, ' '
      endif


      ENDDO

      DEALLOCATE(S_SET)

      close(10)

      ENDDO

      call cpu_time(TIME2)
      time = (TIME2-TIME1)
      print*, "CPU TIME: ", int(time/3600), 'h',
     . int((time/3600-int(time/3600))*60), 'min', 
     . int((time/60-int(time/60))*60), 's'
      print*, ' '

      close(10)

      contains

c------------------------------------------------------------------
c     RANDOM ERDÖS-RÉNYI GRAPH GENERATOR
c------------------------------------------------------------------
      subroutine IRG(N,M,p,NBR,INBR,JJ)
c     This subroutine generates a random graph with p*M edges with 
c     a weight of 1 and (1-p)*M edges with a weight of 1
c     and saves it in the NBR, INBR and JJ arrays.

      implicit none
      integer i, j, k
c     nodes, edges
      integer N, M
c     function to generate a U(0,1) random number
      real*8 genrand_real2 
c     fraction of edges with value 1
      real*8 p
      integer edges_p, edges_n

      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: INBR(:)
      type(multi_array),allocatable:: JJ(:)
       
      edges_p = int(p*M)
      edges_n = int((1-p)*M)

      if (edges_p+edges_n.ne.M) then
            write(*,*)
            print*, 'GRAPH PROBLEM: p value not possible'
            write(*,*)
            stop
      endif

      allocate(NBR(N))
      allocate(INBR(N))
      allocate(JJ(N))

      do i = 1,N
            allocate(NBR(i)%v(0))
            allocate(JJ(i)%v(0))
      enddo

c     generate M/2 edges of weight 1
      k = 0
      do while(k<edges_p)
            i = int(genrand_real2()*N) + 1 
            j = int(genrand_real2()*N) + 1 
            if ((ANY(NBR(i)%v == j).eqv..false.).and.(i.ne.j)) then

                  call AddToList(NBR(i)%v,j)
                  call AddToList(NBR(j)%v,i)
                  call AddToList(JJ(i)%v,1)
                  call AddToList(JJ(j)%v,1)

                  k = k + 1 
            endif
      enddo
       
c     generate M/2 edges of weight -1 
      k = 0
      do while(k<edges_n)
            i = int(genrand_real2()*N) + 1 
            j = int(genrand_real2()*N) + 1 
            if ((ANY(NBR(i)%v == j).eqv..false.).and.(i.ne.j)) then
                  
                  call AddToList(NBR(i)%v,j)
                  call AddToList(NBR(j)%v,i)
                  call AddToList(JJ(i)%v,-1)
                  call AddToList(JJ(j)%v,-1)
            
                  k = k + 1 
            endif
      enddo

c	INBR generation
      do i = 1,N
            allocate(INBR(i)%v(size(NBR(i)%v)))
      enddo
      do i=1,N 
            do j=1,size(NBR(i)%v)
                  do k=1,size(NBR(NBR(i)%v(j))%v)

                        if ((NBR(NBR(i)%v(j))%v(k).eq.i)) then
                              INBR(i)%v(j) = k
                        endif

                  enddo
            enddo
      enddo


      return
      end subroutine

c------------------------------------------------------------------
      real*8 function g(N,M,NBR,JJ,NBR_0,JJ_0)
c     This function calculates the pseudolikelihood for a given TEMP,
c     S_SET and graph defined by NBR,INBR and JJ
      implicit none
      integer N, M
      real*8 sum
      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: JJ(:)
      type(multi_array),allocatable:: NBR_0(:)
      type(multi_array),allocatable:: JJ_0(:)
      integer ip


      sum = 0.d0
      do i = 1,N
            do k = 1,size(NBR_0(i)%v)

                  if (i.lt.NBR_0(i)%v(k)) then
                  ip = floc(NBR(i)%v,NBR_0(i)%v(k))

                  if (ip.eq.0) then
                        sum = sum + (JJ_0(i)%v(k))**2
                  endif
                  if (ip.ne.0) then
                        sum = sum + (JJ_0(i)%v(k)-JJ(i)%v(ip))**2
                  endif
                  endif
            enddo

      enddo

      g = sum/M

      return
      end function
      
c------------------------------------------------------------------
c     CALCULATE LAMBDA 
c------------------------------------------------------------------           
      subroutine class_lambda(N,C,S_SET,NBR,JJ,LAMBDA)
c     This subroutine tabulates the values of the sum inside the 
c     tanh of the pseudolikelihood for a given S_SET and graph 
c     defined by NBR,INBR and JJ.
      implicit none
      integer i, N, k, m
      integer C
      integer sum
      integer S_SET(1:C,1:N), LAMBDA(1:C,1:N)
      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: JJ(:)

      sum = 0

      do m = 1,C
            do i = 1,N

            do k = 1, size(NBR(i)%v)
                  sum = sum + JJ(i)%v(k)*S_SET(m,NBR(i)%v(k))
            enddo

            LAMBDA(m,i) = sum

            sum = 0
            enddo
      enddo
      
      return
      end subroutine

c------------------------------------------------------------------
c     CALCULATE function 
c------------------------------------------------------------------           
      subroutine class_function(zmax,TEMP,funct)
c     This subroutine tabulates the function in the sumatories of 
c     the pseudolikelihood.
      implicit none
      integer i, m, zmax
      real*8 TEMP
      real*8 funct(-1:1,-zmax:zmax)

      do i = -1,1,2
            do m = -zmax,zmax
                  funct(i,m) = log(0.5d0*(1 + i*tanh(m/TEMP)))
            enddo
      enddo
      
      return
      end subroutine

c------------------------------------------------------------------
c     PSEUDOLIKELIHOOD ALGORITHM 
c------------------------------------------------------------------
      subroutine pseudolikelihood(N,C,S_SET,
     . valid,TEMP_F,DPL,NBR,INBR,JJ,funct,zmax,LAMBDA)
c     This subroutine proposes a change of pairwise coupling, 
c     calculates the pseudolikelihood variation of the system (ΔPL) due to it,
c     if ΔPL > 0 then the change is accpeted, else if ΔPL < 0 then the
c     change is accepted with a probability of exp(-ΔPL/k_BT'), where
c     T' is the ficticious temperature
      implicit none
      integer N
      logical valid
      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: INBR(:)
      type(multi_array),allocatable:: JJ(:)
      integer C, S_SET(1:C,1:N)
      real*8 TEMP_F
      logical change
      type(multi_array),allocatable:: newNBR(:)
      type(multi_array),allocatable:: newINBR(:)
      type(multi_array),allocatable:: newJJ(:)
c     function to generate a U(0,1) random number
      real*8 genrand_real2 
      integer zmax, m
      
      real*8 sum
      integer i1,i2,i3,i4
      real*8  DPL
      integer ip,ip2

      real*8 funct(-1:1,-zmax:zmax)
      integer LAMBDA(1:C,1:N), newLAMBDA(1:C,1:N)
      
c..................................................................
      valid = .false.
      change = .false. 

      allocate(newNBR(N))
      allocate(newINBR(N))
      allocate(newJJ(N))

      do i = 1,N
            allocate(newNBR(i)%v(size(NBR(i)%v)))
            allocate(newINBR(i)%v(size(NBR(i)%v)))
            allocate(newJJ(i)%v(size(NBR(i)%v)))
      enddo

      newNBR = NBR
      newINBR = INBR
      newJJ = JJ 

      newLAMBDA = LAMBDA

c..................................................................
c     Random pairwise coupling change
      do while (change.eqv..false.)
      call JJ_change(N,NBR,INBR,JJ,newNBR,newINBR,newJJ,change,
     .              i1,i2,i3,i4)
      enddo
c..................................................................
c     Calculate newLAMBDA

      do m=1,C
      ip = floc(NBR(i1)%v,i3)
      ip2 = floc(newNBR(i1)%v,i3)
      if(ip.eq.0)then
      newLAMBDA(m,i1) = newLAMBDA(m,i1) 
     . + newJJ(i1)%v(ip2)*S_SET(m,i3)
      elseif(ip.ne.0)then
      newLAMBDA(m,i1) = newLAMBDA(m,i1) 
     . + 2*newJJ(i1)%v(ip2)*S_SET(m,i3)
      endif
      ip = floc(NBR(i3)%v,i1)
      ip2 = floc(newNBR(i3)%v,i1)
      if(ip.eq.0)then
      newLAMBDA(m,i3) = newLAMBDA(m,i3) 
     . + newJJ(i3)%v(ip2)*S_SET(m,i1)
      elseif(ip.ne.0)then
      newLAMBDA(m,i3) = newLAMBDA(m,i3) 
     . + 2*newJJ(i3)%v(ip2)*S_SET(m,i1)
      endif
      ip = floc(NBR(i2)%v,i4)
      ip2 = floc(newNBR(i2)%v,i4)
      if(ip2.eq.0)then
      newLAMBDA(m,i2) = newLAMBDA(m,i2) 
     . - JJ(i2)%v(ip)*S_SET(m,i4)
      elseif(ip2.ne.0)then
      newLAMBDA(m,i2) = newLAMBDA(m,i2) 
     . + 2*newJJ(i2)%v(ip2)*S_SET(m,i4)
      endif
      ip = floc(NBR(i4)%v,i2)
      ip2 = floc(newNBR(i4)%v,i2)
      if(ip2.eq.0)then
      newLAMBDA(m,i4) = newLAMBDA(m,i4) 
     . - JJ(i4)%v(ip)*S_SET(m,i2)
      elseif(ip2.ne.0)then
      newLAMBDA(m,i4) = newLAMBDA(m,i4) 
     . + 2*newJJ(i4)%v(ip2)*S_SET(m,i2)
      endif
      enddo
      
c..................................................................
c     calculate DPL
      sum = 0

      if ((i1.ne.i4).and.(i2.ne.i3).and.(i3.ne.i4)) then
      do m = 1,C
            sum = sum + 
     .      funct(S_SET(m,i1),newLAMBDA(m,i1))
     .    - funct(S_SET(m,i1),LAMBDA(m,i1)) +
     .      funct(S_SET(m,i2),newLAMBDA(m,i2))
     .     -funct(S_SET(m,i2),LAMBDA(m,i2)) + 
     .      funct(S_SET(m,i3),newLAMBDA(m,i3)) 
     .     -funct(S_SET(m,i3),LAMBDA(m,i3)) +
     .      funct(S_SET(m,i4),newLAMBDA(m,i4))
     .     -funct(S_SET(m,i4),LAMBDA(m,i4))
      enddo
      DPL = sum/C
      endif


      if (i1.eq.i4) then
      do m = 1,C
            sum = sum + 
     .      funct(S_SET(m,i1),newLAMBDA(m,i1))
     .    - funct(S_SET(m,i1),LAMBDA(m,i1)) +
     .      funct(S_SET(m,i2),newLAMBDA(m,i2))
     .     -funct(S_SET(m,i2),LAMBDA(m,i2)) + 
     .      funct(S_SET(m,i3),newLAMBDA(m,i3)) 
     .     -funct(S_SET(m,i3),LAMBDA(m,i3))
      enddo
      DPL = sum/C
      endif

      if ((i2.eq.i3).or.(i4.eq.i3)) then
      do m = 1,C
            sum = sum + 
     .      funct(S_SET(m,i1),newLAMBDA(m,i1))
     .    - funct(S_SET(m,i1),LAMBDA(m,i1)) +
     .      funct(S_SET(m,i2),newLAMBDA(m,i2))
     .     -funct(S_SET(m,i2),LAMBDA(m,i2)) + 
     .      funct(S_SET(m,i4),newLAMBDA(m,i4))
     .     -funct(S_SET(m,i4),LAMBDA(m,i4))
      enddo
      DPL = sum/C
      endif
c..................................................................
c     check if the change is accpeted
      if (genrand_real2().lt.min(1.d0,exp(DPL/TEMP_F))) then
            valid = .true.
            NBR = newNBR
            INBR = newINBR
            JJ = newJJ
            
            LAMBDA = newLAMBDA

      endif
  
      return
      end subroutine

c------------------------------------------------------------------
      real*8 function PEUDOL(N,C,S_SET,TEMP,NBR,JJ)
c     This function calculates the pseudolikelihood for a given TEMP,
c     S_SET and graph defined by NBR,INBR and JJ
      implicit none
      integer i, N, k, m
      integer C, S_SET(1:C,1:N)
      real*8 TEMP
      real*8 PL, L
      real*8 sum

      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: JJ(:)

      PL = 0.0d0
      L = 0.0d0
      sum = 0.0d0

      do i = 1,N
      do m = 1,C
            do k = 1, size(JJ(i)%v)
                  sum = sum + JJ(i)%v(k)*S_SET(m,NBR(i)%v(k))
            enddo
            
            L = L + log(0.5d0*(1 + S_SET(m,i)*tanh(sum/TEMP)))

            sum = 0.0d0
      enddo
            PL = PL + L/C

            L = 0.0d0
      enddo
      
      PEUDOL = PL 
      
      return
      end function

c------------------------------------------------------------------
c     MAKE A JJ CHANGE 
c------------------------------------------------------------------ 
      subroutine JJ_change(N,NBR,INBR,JJ,newNBR,newINBR,newJJ,change,
     .                     i1,i2,i3,i4)
c     This subroutine changes pairwise coupling of i1-i3 for the
c     pairwise coupling of i2-i4, and generates a the new NBR, INBR
c     and JJ arrays.

      implicit none
      integer N
      integer i, k
      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: INBR(:)
      type(multi_array),allocatable:: JJ(:)
      type(multi_array),allocatable:: newNBR(:)
      type(multi_array),allocatable:: newINBR(:)
      type(multi_array),allocatable:: newJJ(:)
c     function to generate a U(0,1) random number
      real*8 genrand_real2 
      integer i1,k1, i2,k2, i3,k3, i4,k4

      integer ip
      
      logical change, valid, valid2, valid3

      change = .false.
      valid = .false.
      valid2 = .false.
      valid3 = .false.

c     selecting a random node
55    i1 = int(genrand_real2()*N) + 1 
      if (i1.eq.(N+1)) go to 55

c     selecting a random edge of that node, the extra +1 
c     is because we can select a node with value 0
      k1 = int(genrand_real2()*size(NBR(i1)%v)) + 1 + 1 

c     in case size(NBR(i,k))= 0 k should be equal to 1
      if (size(NBR(i1)%v) == 0) then
            k1 = 1
      endif

c     if the max edges per node is exceded we deny the change
      if (k1 > N-1) then
            change = .false.
            return
      endif

c------------------------------------------------------------------
c     if JJ(i1,k1) = 0 we want a 0 <----> 1,-1 change
c------------------------------------------------------------------
      if (k1 > size(NBR(i1)%v)) then

c     diferent random node selection i2
c..................................................................
      do while (valid.eqv..false.)

56    i2 = int(genrand_real2()*N) + 1
      if (i2.eq.(N+1)) go to 56

      if ((i1.ne.i2).and.(size(NBR(i2)%v).gt.1)) then 
            valid =  .true.
      endif
      enddo

c     i2 can't be the same as i
c     to avoid problems i2 needs to have more than 1 neighbour
c     otherwise a node could end up with 0 edges
c..................................................................
c     selecting a random neigbour of the node i2 
60    k2 = int(genrand_real2()*size(NBR(i2)%v)) + 1
      if (k2.eq.(size(NBR(i2)%v)+1)) go to 60

c..................................................................
c     selecting a node i3 that is not a neigbour of i already
      do while (valid2.eqv..false.)

57    i3 = int(genrand_real2()*N) + 1
      if (i3.eq.(N+1)) go to 57

c     this allows to see if i3 is already a neigbour of i
      ip = floc(NBR(i1)%v,i3)
c     ip is 0 if i3 is not in NBR(i)%v

      if ((ip.eq.0).and.(i1.ne.i3)
     . .and.(size(NBR(i3)%v).lt.N-1)) then
            valid2 =  .true.
      endif
      enddo

      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

c..................................................................
c     adding new edge

      call AddToList(newNBR(i1)%v,i3)
      call AddToList(newNBR(i3)%v,i1)

      call AddToList(newJJ(i1)%v,JJ(i2)%v(k2))
      call AddToList(newJJ(i3)%v,JJ(i4)%v(k4))
c..................................................................
c     removing the old edge

      call RmvOffList(newNBR(i2)%v,k2)
      call RmvOffList(newNBR(i4)%v,k4)

      call RmvOffList(newJJ(i2)%v,k2)
      call RmvOffList(newJJ(i4)%v,k4)
c..................................................................
c     updating INBR
      newINBR = newNBR
      do i=1,N 
            do j=1,size(newNBR(i)%v)
                  do k=1,size(newNBR(newNBR(i)%v(j))%v)
                        if ((newNBR(newNBR(i)%v(j))%v(k).eq.i)) then
                              newINBR(i)%v(j) = k
                        endif
                  enddo
            enddo
      enddo
c..................................................................
      change = .true.
      return
      endif

c------------------------------------------------------------------
c     if JJ(i1,k1) = -1,1 we want a +-1 <----> +-1 change
c------------------------------------------------------------------
      if (k1 <= size(NBR(i1)%v)) then
      
c     diferent random node selection i2
c..................................................................
      do while (valid.eqv..false.)

59    i2 = int(genrand_real2()*N) + 1
      if (i2.eq.(N+1)) go to 59

      if ((i1.ne.i2).and.(size(NBR(i2)%v).gt.1)) then
      !if (i1.ne.i2) then
            valid =  .true.
      endif
      enddo
c     i2 can't be the same as i
c..................................................................
c     selecting a random neigbour of the node i2 
      do while (valid2.eqv..false.)

61    k2 = int(genrand_real2()*size(NBR(i2)%v)) + 1
      if (k2.eq.(size(NBR(i2)%v)+1)) go to 61

      if (NBR(i2)%v(k2).ne.i1) then
            valid2 =  .true.
      endif
      enddo
c     if both edges have the same value there is no change
      if (JJ(i1)%v(k1).eq.JJ(i2)%v(k2)) then
            change = .false.
            return
      endif
c..................................................................
c     changing the edge
      i3 = NBR(i1)%v(k1)
      k3 = INBR(i1)%v(k1)

      i4 = NBR(i2)%v(k2)
      k4 = INBR(i2)%v(k2)

      newJJ(i1)%v(k1) = JJ(i2)%v(k2)
      newJJ(i3)%v(k3) = JJ(i2)%v(k2)

      newJJ(i2)%v(k2) = JJ(i1)%v(k1)
      newJJ(i4)%v(k4) = JJ(i1)%v(k1)

      change = .true.
      return
      endif

      return
      end subroutine

c------------------------------------------------------------------
c     Find location of a certain value in a array 
c------------------------------------------------------------------

      integer function floc(array,value)

      implicit none
      integer,allocatable:: array(:)
      integer value, r
      logical out

      floc = 0
      do r = 1, size(array)
            if (array(r) .eq. value) then
                  floc = r
                  out = .true.
            endif
      end do

      end function

c------------------------------------------------------------------
c     Add element into array
c------------------------------------------------------------------

      subroutine AddToList(list, element)

      IMPLICIT NONE

      integer :: i, isize
      integer, intent(in) :: element
      integer, dimension(:), allocatable, intent(inout) :: list
      integer, dimension(:), allocatable :: clist


      if(allocated(list)) then

            isize = size(list)
            allocate(clist(isize+1))
            do i=1,isize          
            clist(i) = list(i)
            end do
            clist(isize+1) = element

            deallocate(list)
            call move_alloc(clist, list)

      else
            allocate(list(1))
            list(1) = element
      end if


      return
      end subroutine

c------------------------------------------------------------------
c     Remove specific index off array
c------------------------------------------------------------------

      subroutine RmvOffList(list, index)

      IMPLICIT NONE

      integer :: i, isize
      integer, intent(in) :: index
      integer, dimension(:), allocatable, intent(inout) :: list
      integer, dimension(:), allocatable :: clist


      if(allocated(list)) then

            isize = size(list)
            allocate(clist(isize-1))

            do i=1,index-1
                  clist(i) = list(i)
            end do
            do i=index,isize-1
                  clist(i) = list(i+1)
            end do

            deallocate(list)
            call move_alloc(clist, list)

      end if

      return
      end subroutine
      
      end

c------------------------------------------------------------------

c
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
