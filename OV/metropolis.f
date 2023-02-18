c     METROPOLIS SAMPLE GENERATOR
c     LLUÍS
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
c     number of graphs simulated for every p value
      integer NSEEDS
      parameter(NSEEDS = 100)
c------------------------------------------------------------------
c------------------------------------------------------------------
c     seed number
      integer SEED, SEEDini
      integer i, j, k, IMC, IPAS
c     multidimensional array type
      type multi_array
      integer,allocatable::v(:)
      end type multi_array
c     principal arrays
      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: INBR(:)
      type(multi_array),allocatable:: JJ(:)
c     simulation variables
      integer S1(1:N), S2(1:N)
      logical valid1, valid2
      real*8 DE1, DE2
      real*8 ENE1, ENE2
      real*8 q
c     estimate time variables
      real*4 TIME1, TIME2, time
c     function to generate a U(0,1) random number
      real*8 genrand_real2
c------------------------------------------------------------------
c------------------------------------------------------------------
      integer C
c     sample size (number of spin config)
      parameter(C = 2*(MCTOT-MCINI)/SC)
c     size of the histogram boxes
      real*8 hstep
      parameter(hstep = 0.01d0)
c     set of sambles
      integer,allocatable:: S_SET(:,:)
      real*8 a,b, g
      real*8 ini, fin
      parameter(ini = 0.d0)
      parameter(fin = 1.d0)

      integer sze,l
      parameter(sze = (fin-ini)/hstep)
      real*8 histo_list(0:sze)
      real*8 histo_tot(0:sze)

      character(4) str5
      character(3) str6
c------------------------------------------------------------------
c     number of p values to sumulate
      integer npvalues, np
      real*8 pvalues(1:7)
c------------------------------------------------------------------

c------------------------------------------------------------------
c     START
c------------------------------------------------------------------

      SEEDini = 100
      
c     All p values
      pvalues = (/0.d0,0.15d0,0.3d0,0.5d0,0.7d0,0.85d0,1.d0/)
      npvalues = size(pvalues)

      call cpu_time(TIME1)

      print*,' '
      print*,'N=',N,'z=',Z,'C=',C,'T=',TEMP

c¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
c     Do for all p values
      DO np=1,npvalues

      p = pvalues(np)

      write( str5,'(f4.2)')   p
      
      print*,' '
      print*,'p# =', p
      print*,' '

      histo_tot = 0

      allocate(S_SET(1:C,1:N))

c¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
c     do for all seeds
      DO SEED = SEEDini,SEEDini+NSEEDS-1
      
      write( str6,'(i3)')     SEED

      print*,'SEED# =', SEED
      print*,' '
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call init_genrand(SEED)

c     Initial NBR, INBR and JJ random arrays generation
      call IRG(N,M,p,NBR,INBR,JJ)

      print*,'Graph Generated'
      print*,' '

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c     Spin config data file
      open(unit = 3,file = 'spin_config_'//str5//'_'//str6//'.dat')

c     Genaration of two random initial spin configuration
      do i = 1,N
            S1(i) = int(2*mod(int(2*genrand_real2()),2) - 1)
            S2(i) = int(2*mod(int(2*genrand_real2()),2) - 1)
      end do

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Initial energy calculation
      ENE1 = ENERG(S1,N,NBR,JJ)
      ENE2 = ENERG(S2,N,NBR,JJ)

      print*,'MC starts'
      print*,' '

      k = 1
c     MC simulation
      do IMC = 1, MCTOT
      
            do IPAS = 1, N
                  call metropolis(S1,N,valid1,TEMP,DE1,NBR,JJ)
                  call metropolis(S2,N,valid2,TEMP,DE2,NBR,JJ)
                  if (valid1) then
                        ENE1 = ENE1 + DE1
                  end if
                  if (valid2) then
                        ENE2 = ENE2 + DE2
                  end if
            end do

c           Extract the spin configuratinons every SC montecarlo steps
            if ((IMC.gt.MCINI).and.(SC*(IMC/SC).eq.IMC)) then

                  write(3,*) S1
                  write(3,*) S2

            endif

      enddo

      close(3)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      print*, 'Samples Generated'
      print*,' '

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      open(unit = 1,file = 'spin_config_'//str5//'_'//str6//'.dat')
      do i = 1, C
            read(1,*) S_SET(i,:)
      enddo
      close(1)

      histo_list = 0

      open(unit = 2,file = 'histogram_'//str5//'_'//str6//'.dat')
      open(unit = 5,file = 'histogram_tot_'//str5//'.dat')

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Generating 1-seed and nseed-seed histograms data
      do i = 1,C
      do j = 1,C
            if (j.gt.i) then

      q = 0
      do k = 1,N
            q = q + S_SET(i,k)*S_SET(j,k)
      enddo

      q = abs(q)/N

      a = ini 
      b = ini + hstep

      do l = 0,sze
            if((q.ge.a).and.(q.le.b)) then
                  histo_list(l) = histo_list(l) + 1
                  histo_tot(l) = histo_tot(l) + 1
            endif
            a = a + hstep
            b = b + hstep
      enddo
      
            endif      
      enddo
      enddo

      histo_list = histo_list/sum(histo_list)

      g = ini
      do l = 0,sze
      write(2,"(F6.3,2X,F15.8)") g, histo_list(l)
      g = g + hstep
      enddo

      close(2)

      print*,'Histogram Generated'
      print*,' '
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do i = 1,N
            DEALLOCATE(NBR(i)%v)
            DEALLOCATE(INBR(i)%v)
            DEALLOCATE(JJ(i)%v)
      enddo
      DEALLOCATE(NBR)
      DEALLOCATE(INBR)
      DEALLOCATE(JJ)


      if ((SEED.eq.SEEDini).and.(p.eq.pvalues(1))) then
      call cpu_time(TIME2)
      time = (TIME2-TIME1)*NSEEDS*npvalues
      print*, "ESTIMATED TIME: ", int(time/3600), 'h',
     . int((time/3600-int(time/3600))*60), 'min', 
     . int((time/60-int(time/60))*60), 's'
      print*, ' '
      endif

c__________________________________________________________________
      ENDDO

      close(8)

      histo_tot = histo_tot/(sum(histo_tot))
      
      g = ini
      do l = 0,sze
      if (histo_tot(l).gt.0) then
      write(5,"(F6.3,2X,F15.8)") g, histo_tot(l)
      endif
      g = g + hstep
      enddo

      DEALLOCATE(S_SET)

c__________________________________________________________________
      ENDDO

      call cpu_time(TIME2)
      time = (TIME2-TIME1)
      print*, "CPU TIME: ", int(time/3600), 'h',
     . int((time/3600-int(time/3600))*60), 'min', 
     . int((time/60-int(time/60))*60), 's'
      print*, ' '

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
c     METROPOLIS ALGORITHM 
c------------------------------------------------------------------
      subroutine metropolis(S,N,valid,TEMP,DE,NBR,JJ)
c     This subroutine proposes a change of spin in a random node, 
c     calculates the energy variation of the system (ΔH) due to it,
c     if ΔH < 0 then the change is accpeted, else if ΔH > 0 then the
c     change is accepted with a probability of exp(-ΔH/k_BT).

      implicit none
      integer i, k, N
      integer S(1:N)
      real*8 genrand_real2
      logical valid

      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: JJ(:)
      real*8 suma, TEMP, DE
      
      valid = .false.

c     Random node selection

55      i = int(genrand_real2()*N) + 1

      if (i.eq.(N+1)) go to 55

c     calcul of ΔH
      suma = 0
      do k=1,size(NBR(i)%v)
            suma = suma + JJ(i)%v(k)*S(NBR(i)%v(k))
      enddo
      DE = 2*S(i)*suma

c     check if the change is accpeted
      if (genrand_real2().lt.min(1.d0,exp(-DE/TEMP))) then
            S(i) = -S(i)
            valid = .true.
      endif


      return
      end subroutine

c------------------------------------------------------------------
      real*8 function ENERG(S,N,NBR,JJ)
c     This function calculates the energy of the system given an 
c     S configuration

      implicit none
      integer S(1:N)
      integer i, N
      real*8 ENE

      type(multi_array),allocatable:: NBR(:)
      type(multi_array),allocatable:: JJ(:)

      ENE = 0.0d0
      
      do i = 1,N
            do k=1,size(NBR(i)%v)
                  ENE = ENE - S(i)*S(NBR(i)%v(k))*JJ(i)%v(k)
            enddo
      end do

      ENERG =  ENE/2
      
      return
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
