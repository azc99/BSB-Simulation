      program neighbor_modes
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER config_num,data_interval,num_points
      CHARACTER*30 ARG
      DOUBLE PRECISION LD(4,4)
      DOUBLE PRECISION tau,del_t,take_data
      INTEGER start,end,count_rate
      DOUBLE PRECISION elapsed_time
      
      
c---- Get configuration number from command line
      call getarg(1, ARG)
      read(ARG, *) config_num
      
      write(6,*) 'Config', config_num+1, 'starting...'

c     MARCH 18      
c      LD(1,1) = 0.05652497
c      LD(2,1) = 0.04384723
c      LD(3,1) = 0.04012076

c      LD(1,2) = 0.05475684
c      LD(2,2) = 0.01244233
c      LD(3,2) = 0.05948196

c      LD(1,3) = 0.03306456
c      LD(2,3) = 0.06618305
c      LD(3,3) = 0.03797424

c      LD = 0.d0
c      LD(3,3) = LD(1,2)


C     2X2
c      LD(1,1) = 0.05625
c      LD(1,2) = 0.08536
      
c      LD(2,1) = 0.06919
c      LD(2,2) = -0.06072

c     3x3
c      LD(1,1) = 0.0625
c      LD(1,2) = 0.0776
c      LD(1,3) = -0.0457

c      LD(2,1) = 0.0629
c      LD(2,2) = -2.77d-6
c      LD(2,3) = 0.0909

c      LD(3,1) = 0.0625
c      LD(3,2) = -0.0776
c      LD(3,3) = -0.0457 

c     4x4

      LD(1,1) = 0.0542
      LD(1,2) = 0.0735
      LD(1,3) = 0.0551
      LD(1,4) = 0.0239

      LD(2,1) = 0.0541
      LD(2,2) = 0.0234
      LD(2,3) = -0.0552
      LD(2,4) = -0.0753

      LD(3,1) = 0.0541
      LD(3,2) = -0.023
      LD(3,3) = -0.0552
      LD(3,4) = 0.0753

      LD(4,1) = 0.0542
      LD(4,2) = -0.0735
      LD(4,3) = 0.0551
      LD(4,4) = -0.0239

c     5x5

c      LD(1,1) = 0.0486
c      LD(1,2) = 0.0694
c      LD(1,3) = -0.0586
c      LD(1,4) = 0.0335
c      LD(1,5) = 0.0119

c      LD(2,1) = 0.0482
c      LD(2,2) = 0.0330
c      LD(2,3) = 0.0307
c      LD(2,4) = -0.0705
c      LD(2,5) = -0.0526

c      LD(3,1) = 0.0481
c      LD(3,2) = 1.12d-5
c      LD(3,3) = 0.0569
c      LD(3,4) = 1.66d-5
c      LD(3,5) = 0.0814

c      LD(4,1) = 0.0483
c      LD(4,2) = -0.0330
c      LD(4,3) = 0.0307
c      LD(4,4) = 0.0705
c      LD(4,5) = -0.0526

c      LD(5,1) = 0.0487
c      LD(5,2) = -0.0694
c      LD(5,3) = -0.0586
c      LD(5,4) = -0.0335
c      LD(5,5) = 0.0119

      
c      LD(1,1) = 0.05538204
c      LD(1,2) = 0.05541157
c      LD(1,3) = 0.03148604

c      LD(2,1) = 0.0459863
c      LD(2,2) = 0.00905555
c      LD(2,3) = 0.06678407

c      LD(3,1) = 0.04008227
c      LD(3,2) = 0.06083569
c      LD(3,3) = 0.03813998


c      LD(1,1) = 0.0537462
c      LD(1,2) = 0.0628706
c      LD(1,3) = 0.02973594

c      LD(2,1) = 0.04490783
c      LD(2,2) = 0.01107542
c      LD(2,3) = 0.0631403

c      LD(3,1) = 0.06677761
c      LD(3,2) = 0.06539414
c      LD(3,3) = 0.06504527


c      LD(1,1) = 0.05690
c      LD(1,2) = 0.05408
c      LD(1,3) = 0.03135

c      LD(2,1) = 0.04706
c      LD(2,2) = 0.01146
c      LD(2,3) = 0.07078

c      LD(3,1) = 0.04187
c      LD(3,2) = 0.06359
c      LD(3,3) = 0.03792 

      tau = 800.d-6
      del_t = 1.d-6
      take_data = 25.d-6
      data_interval = NINT(take_data/del_t)
      num_points = NINT(tau/take_data) + 1

      call system_clock(start,count_rate)
      call specific_time(4,4  ,tau,config_num,LD,del_t,data_interval,
     .                   num_points)
      call system_clock(end)
      write (6,*) "time: ", real(end - start) / real(count_rate), " s"
      end 
     
      subroutine specific_time(j,k,tau,m,LD,del_t,data_interval,
     .                         num_points)
      IMPLICIT NONE
      INTEGER j,k,x,y,z,a,b,c,d,ii,iii,unit_num,m,ts,total_ts
      COMPLEX*16 PSI_IN(0:2**j*4**k-1),H(0:7,0:7),H_TROT(0:7,0:7)
      DOUBLE PRECISION LD(j,k),pi,MODE(k,4),
     .                 SPIN(j,2),tau,P(7),PROB
      COMPLEX*16 i,plus,pluscre,plusann,minus,minann,mincre,H_C(0:1,0:1)
      INTEGER rate,low,high,driving_mode,cnt,num_points,
     .        data_interval, mod_result
      DOUBLE PRECISION del_t,t1,t2,t3,sumsq,detune(7,7)
      DOUBLE PRECISION RABI(7),w_k(7),w(7),wqbt(7),w_access(7)
      CHARACTER*50 F1,F2,F3,F4
      
      pi = 4.d0*ATAN(1.d0)
      i = (0.d0,1.d0)

c---- generate input state vector

      MODE = 0.d0; SPIN = 0.d0

c     SPIN(x,y) -> ion j+1-x, spin y      
      do ii = 1, j
c          SPIN(ii,1) = sqrt(995.d-3)
c          SPIN(ii,2) = sqrt(5.d-3)
           SPIN(ii,1) = 1.d0
      end do
c     MODE(x,y) = mo9de k+1-x, phonon y-1
      do ii = 1, k
c         MODE(ii,2) = sqrt(1.d-1)
c         MODE(ii,1) = sqrt(9.d-1)
          MODE(ii,1) = 1.d0
      end do

      PSI_IN = 1.d0
      CALL GEN_PSI(j,k,SPIN,MODE,PSI_IN)       

c---- input values

c      detune(1,1) = 8.79038883e-06
c      detune(1,2) = 9.40519099e-06
c      detune(1,3) = 1.03265252e-05

c      detune(2,1) = 1.09229815e-05
c      detune(2,2) = 8.61293440e-06
c      detune(2,3) = 6.63851244e-06

c      detune(3,1) = 9.19494751e-06
c      detune(3,2) = 8.27223757e-06
c      detune(3,3) = 1.00397020e-05

      detune = 1.d-5

      RABI(1) = 1/(11.32)
      RABI(2) = 1/(9.38)
      RABI(3) = 1/(9.59)
      RABI(4) = 1/(9.88)
      RABI(5) = 1/(8.77)
      RABI(6) = 0.02
      RABI(7) = 0.02

      RABI = RABI*(pi/2)
 
      w_access(1) = (208.499094 - 3.1407)*2*pi
      w_access(2) = (208.499094 - 3.1115)*2*pi
      w_access(3) = (208.499094 - 3.0687)*2*pi
      w_access(4) = (208.499094 - 3.0155)*2*pi
      w_access(5) = (208.499094 - 2.9526)*2*pi
      

      w(1) = (w_access(MOD(k-m,k)+1)) - detune(m+1,1)
      w(2) = (w_access(MOD(k-m+1,k)+1)) - detune(m+1,2)
      w(3) = (w_access(MOD(k-m+2,k)+1)) - detune(m+1,3)
      w(4) = (w_access(MOD(k-m+3,k)+1)) - detune(m+1,4)
      w(5) = (w_access(MOD(k-m+4,k)+1)) - detune(m+1,5)

      do ii = 0, j-1
         write(6,*) MOD(k-m+ii,k)+1
      end do

      wqbt(1) = 208.499094*2*pi
      wqbt(2) = 208.499094*2*pi
      wqbt(3) = 208.499094*2*pi
      wqbt(4) = 208.499094*2*pi
      wqbt(5) = 208.499094*2*pi

      w_k(1) = abs(w_access(1)-wqbt(1))
      w_k(2) = abs(w_access(2)-wqbt(2))
      w_k(3) = abs(w_access(3)-wqbt(3))
      w_k(4) = abs(w_access(4)-wqbt(4)) 
      w_k(5) = abs(w_access(5)-wqbt(5))


c--- set to proper units

      RABI(1) = RABI(1) * (10**6)
      RABI(2) = RABI(2) * (10**6)
      RABI(3) = RABI(3) * (10**6)
      RABI(4) = RABI(4) * (10**6)
      RABI(5) = RABI(5) * (10**6)
      RABI(6) = RABI(6) * (10**6)
      RABI(7) = RABI(7) * (10**6)

      w_k(1) = w_k(1) * (10**6)
      w_k(2) = w_k(2) * (10**6)
      w_k(3) = w_k(3) * (10**6)
      w_k(4) = w_k(4) * (10**6)
      w_k(5) = w_k(5) * (10**6)
      w_k(6) = w_k(6) * (10**6)
      w_k(7) = w_k(7) * (10**6)

      w(1) = w(1) * (10**6)
      w(2) = w(2) * (10**6)
      w(3) = w(3) * (10**6)
      w(4) = w(4) * (10**6)
      w(5) = w(5) * (10**6)
      w(6) = w(6) * (10**6)
      w(7) = w(7) * (10**6)

      wqbt(1) = wqbt(1) * (10**6)
      wqbt(2) = wqbt(2) * (10**6)
      wqbt(3) = wqbt(3) * (10**6)
      wqbt(4) = wqbt(4) * (10**6)
      wqbt(5) = wqbt(5) * (10**6)
      wqbt(6) = wqbt(6) * (10**6)
      wqbt(7) = wqbt(7) * (10**6)
      
c--------------------------------------
       
      total_ts = NINT(tau/del_t)
      write (6,*) "total timesteps: ", total_ts

      cnt = 0
      t1 = 0.d0
      sumsq = 0.d0

      
      write(F4, '(A,I0,A)') 'config',m,'_total.txt'

      open(unit=11, file=F4, status='unknown')

      do x=1,j
         CALL find_prob(PSI_IN,j,k,x,1,PROB)
         P(x) = PROB
      end do
      
      if (m .EQ. 0) then
         write(11,*) t1,(P(x), x = 1, j)
      else
         write(11,*) (P(x), x = 1, j)
      end if
      write(50,*) t1,P(1)

      do ts = 1,total_ts

         t2 = t1 + del_t
c         write(6,*) "spin only terms"
        
c         do x = 1,j
c            CALL sigma_z(j,k,x,t2,t1,PSI_IN,RABI,w,wqbt,w_k)
c         end do

         do x = 1,j
           CALL spin_only(j,k,PSI_IN,x,m,t2,t1,LD,RABI,w,wqbt,w_k)
         enddo

         do x = 1,j
            CALL sigma_z(j,k,x,t2,t1,PSI_IN,RABI,w,wqbt,w_k)
         end do

c         write(6,*) "spin mode terms"
         do x = 1,j
            
            driving_mode = MOD((k-m+x-1),k)+1
c            write(6,*) "dm: ", driving_mode
            low = MAX(driving_mode-1,1)
            high = MIN(driving_mode+1,k)
            do y = 1,k
              CALL spin_mode(j,k,x,y,t2,t1,PSI_IN,LD,RABI,w,wqbt,w_k)
            end do
         end do

         do x = 1,j
            do y = 1,k
               do z = 1,k

               if (y .EQ. z) then
c                  CALL spin_single(j,k,x,y,t2,t1,PSI_IN,
c     .                             LD,RABI,w,wqbt,w_k)
               else
c                 CALL spin_mode_mode(j,k,x,y,z,t2,t1,PSI_IN,
c     .                               LD,RABI,w,wqbt,w_k)
               end if
 
               end do   
            end do      
         end do

         do x = 1,j
            do y = 1,j
               if (x .NE. y) then
                 do z = 1,k
c                    call spin_spin_mode(j,k,x,y,z,t2,t2,t1,PSI_IN,
c     .                                  LD,RABI,w,wqbt,w_k)
                 end do
                else
                  do z = 1,k
c                     call identity_mode(j,k,x,z,t2,t2,t1,PSI_IN,
c     .                                   LD,RABI,w,wqbt,w_k)
                  end do
               end if 
            end do
         end do


         do x = 1,j
            do y = 1,j
               if (x .EQ. y) then
C                  CALL sigma_z(j,k,x,t2,t1,PSI_IN,RABI,w,wqbt,w_k)
               end if
            end do
         end do   

         do a = 1,j
            do b = 1,j
               do c = 1,k
                  do d = 1,k
                  if (a .EQ. b .AND. c .EQ. d) then
c                    call identity_single(j,k,a,c,t3,t2,t1,PSI_IN,
c     .                                   LD,RABI,w,wqbt,w_k)

                  else if (c .EQ. d) then
c                    call spin_spin_single(j,k,a,b,c,t3,t2,t1,PSI_IN,
c     .                                    LD,RABI,w,wqbt,w_k)   
                  else if (a .EQ. b) then
c                    call identity_mode_mode(j,k,a,c,d,t3,t2,t1,PSI_IN,
c     .                                      LD,RABI,w,wqbt,w_k)
                  else
c                    call spin_spin_mode_mode(j,k,a,b,c,d,t3,t2,t1,
c     .                                       PSI_IN,LD,RABI,w,wqbt,w_k) 
                  end if   
                  
                  end do
               end do
            end do      
         end do


c         do a = 1,j
c           do b = 1,j
c             do c = 1,k
c               do d = 1,k
c                 do e = 1,k
c                 end do
c               end do
c             end do
c           end do
c         end do     

c        do x = 1,j
c           do y = 1,j
c              call spin_spin(j,k,x,y,t2,t2,t1,PSI_IN,RABI,w,wqbt,w_k)
c           end do
c        end do

         if (MOD(ts,data_interval) .EQ. 0) then
            write (6,*) "config: ", m, "t: ", t2
            do ii = 1, j
                CALL find_prob(PSI_IN,j,k,ii,1,PROB)
                P(ii) = PROB
            end do
            if (m .EQ. 0) then 
               write(11,*) t2*1e6,(P(ii), ii = 1, j)
            else
               write(11,*) (P(ii), ii = 1, j)
            end if

            write(50,*) t2*1e6, P(1)
            cnt = cnt + 1
         end if
         t1 = t1 + del_t  
      end do

      close(11)

      return
      end

c===================================================================
c                 APPLY SPIN ONLY TERMS TO PSI
C===================================================================

      subroutine spin_only(j,k,PSI_IN,x,mode,t2,t1,
     .                                LD,RABI,w,wqbt,w_k)
      INTEGER j,k,x,mode
      DOUBLE PRECISION LD(j,k),RABI(7),wj,w_k(7),w(7),wqbt(7),t2,t1
      COMPLEX*16 plus,minus
      DOUBLE PRECISION plusx,minusy
      COMPLEX*16 i,H(0:1,0:1),SIGX(0:1,0:1),SIGY(0:1,0:1),
     .           H_C(0:1,0:1),PSI_IN(0:2**j*4**k-1),SIG(0:1,0:1)
      
      i = (0.d0,1.d0)
      wj = wqbt(x)-w(x)

      plus = RABI(x)*(i/wj)
     .            *(exp(-(i)*wj*t2)
     .            -exp(-(i)*wj*t1))
      minus = RABI(x)*((-i)/wj)
     .              *(exp(i*wj*t2)-
     .                exp(i*wj*t1))

c      write(6,*) "RABI*(i/wj)", abs(RABI(x)*(i/wj))
c      write(6,*) "exponential term", abs((exp(-(i)*wj*t2)
c     .            -exp(-(i)*wj*t1)))
c      write(6,*) "full term", abs(plus)

      plusx = RABI(x)*(1/wj)
     .        *(sin(wj*t2)-sin(wj*t1))
      minusy = RABI(x)*(1/wj)
     .        *(cos(wj*t1)-cos(wj*t2)) 

c      H(0,0) = 0.d0
c      H(1,1) = 0.d0
c      H(0,1) = minus
c      H(1,0) = plus

c      CALL TAYLOR_APPROXIMATION(1,H,50,-(0.d0,1.d0),H_C)

c      SIGX(0,0) = cos(plusx/2.d0)
c      SIGX(0,1) = i*sin(plusx/2.d0)
c      SIGX(1,0) = i*sin(plusx/2.d0)
c      SIGX(1,1) = cos(plus/2.d0)

c      SIGY(0,0) = cos(minusy)
c      SIGY(0,1) = -sin(minusy)
c      SIGY(1,0) = sin(minusy)
c      SIGY(1,1) = cos(minusy)

      SIG(0,0) = cos(sqrt(plusx**2+minusy**2))
      SIG(1,1) = cos(sqrt(plusx**2+minusy**2))
      SIG(1,0) =
     .(sin(sqrt(plusx**2+minusy**2))/sqrt(plusx**2+minusy**2))*
     .(-i*plusx+minusy)

      SIG(0,1) = 
     .(sin(sqrt(plusx**2+minusy**2))/sqrt(plusx**2+minusy**2))*
     .(-i*plusx-minusy)
      

c      FIN = MATMUL(SIGX,SIGY)
      
      CALL carrier_operator(SIG,PSI_IN,j,k,x)


      return
      end

      subroutine carrier_operator(SIG,PSI,j,k,ion)
      IMPLICIT NONE
      INTEGER a,b,c,j,k,ion
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER ip,i_leading,i_trailing,itarg,i
      INTEGER IND(0:1)
      COMPLEX*16 INP(0:1),OUTP(0:1),SIGX(0:1,0:1),PSI(0:2**j*4**k-1)
      complex*16 SIGY(0:1,0:1),H(0:1,0:1),SIG(0:1,0:1) 

c      itarg = 1
c      do i=0,2**(j+2*k-1)-1
c         AA(0) = AA(2*i); AA(1) = AA(2*i+1)
c          write(6,*) 2*i,2*i+1
         ! Do your thing with AA
c         i_leading = i/2**itarg; i_trailing = mod(i,2**itarg)
c         ii = i_leading*2**(itarg+1) + i_trailing
c         write(6,*) ii,ii+2**itarg   
c         B(ii) = AA(0); B(ii+2**itarg) = AA(1);
c      enddo
c      stop


      ip = j+k*2-ion
      do ii =0, 2**(j+2*k-1)-1
            call gen_other_bits_carrier(j+2*k-1,ip,ii,base_bits)
            idx = 0
         do a = 0,1
            bit_combo = base_bits
           
            if (a .EQ. 1) bit_combo = IBSET(bit_combo,ip)
            IND(idx) = bit_combo
            idx = idx + 1
         end do

         do a = 0,1
            INP(a) = PSI(IND(a))
         end do
      
c         OUTP = MATMUL(SIGX,INP)
c         OUTP = MATMUL(SIGY,OUTP)
c         OUTP = MATMUL(SIGX,OUTP)
 
          OUTP = MATMUL(SIG,INP)
c         OUTP = MATMUL(H,INP)
           
         do a = 0,1
            PSI(IND(a)) = OUTP(a)
         end do

      end do
      
      return
      end


      subroutine gen_other_bits_carrier(nbits,ip,num,result)
      IMPLICIT NONE
      INTEGER nbits,ip,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
        if (bit_pos .NE. ip) then
            if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
         end if
      end do

      return
      end 

c===============================================================

c===============================================================
c                 APPLY SPIN x MODE TERMS
c===============================================================

      subroutine spin_mode(j,k,x,y,t2,t1,PSI_IN,LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y,a,b,c
      DOUBLE PRECISION wj,w(7),w_k(7),RABI(7),LD(j,k),wqbt(7),t2,t1
      DOUBLE PRECISION plusannx,plusanny,pluscrex,pluscrey,minannx,
     .                 minanny,mincrex,mincrey
      COMPLEX*16 plusann,pluscre,minann,mincre
      COMPLEX*16 i,H(0:7,0:7),PSI_IN(0:2**j*4**k-1),EVENX(0:7,0:7)
      COMPLEX*16 ODDX(0:7,0:7),EVENY(0:7,0:7),ODDY(0:7,0:7),H_C(0:7,0:7)
      COMPLEX*16 ODD(0:7,0:7),EVEN(0:7,0:7)
      

c      if (x .EQ. 1 .AND. y .EQ. 1) then
c          write(6,*) RABI(x), LD(x,y)
c          stop
c      end if
       
      
      i = (0.d0,1.d0)
      wj = wqbt(x)-w(x)

      plusann = i*RABI(x)*LD(x,y)*(i/(w_k(y)+wj))
     .            *(exp((-i)*(w_k(y)+wj)*t2)-
     .              exp((-i)*(w_k(y)+wj)*t1))
      pluscre = i*RABI(x)*LD(x,y)*((-i)/(w_k(y)-wj))
     .           *(exp(i*((w_k(y)-wj)*t2))-
     .             exp(i*((w_k(y)-wj)*t1)))
      minann = (-i)*RABI(x)*LD(x,y)*((-i)/(wj-w_k(y)))*
     .            (exp(i*(wj-w_k(y))*t2)-
     .             exp(i*(wj-w_k(y))*t1))
      mincre = (-i)*RABI(x)*LD(x,y)*((-i)/(w_k(y)+wj))*
     .            (exp(i*(w_k(y)+wj)*t2)-
     .             exp(i*(w_k(y)+wj)*t1))

      
c      write(6,*) "RABI*LD*(i/wj-wk)", 
c     .            abs((-i)*RABI(x)*LD(x,y)*((-i)/(wj-w_k(y))))
c      write(6,*) "exponential term:", abs(exp(i*(wj-w_k(y))*t2)
c     .            -exp(i*(wj-w_k(y))*t1))
c      write(6,*) "full term", abs(minann)

c      plusann = 1.d-10
c      mincre = 1.d-10

      plusannx = -REAL(plusann)
      plusanny = -AIMAG(plusann)
      pluscrex = -REAL(pluscre)
      pluscrey = -AIMAG(pluscre)
      minannx = -REAL(minann)
      minanny = -AIMAG(minann)
      mincrex = -REAL(mincre)
      mincrey = -AIMAG(mincre)

c      H = 0.d0

c      H(0,5) = pluscre
c      H(2,7) = sqrt(3.d0)*pluscre
c      H(5,0) = minann
c      H(7,2) = sqrt(3.d0)*minann
c      H(1,4) = plusann
c      H(3,6) = sqrt(3.d0)*plusann
c      H(4,1) = mincre
c      H(6,3) = sqrt(3.d0)*mincre

c      H(1,6) = sqrt(2.d0)*pluscre
c      H(6,1) = sqrt(2.d0)*minann
      
c      H(2,5) = sqrt(2.d0)*plusann
c      H(5,2) = sqrt(2.d0)*mincre

      EVEN = 0.d0
      ODD = 0.d0

      do a = 0,7
c         EVENX(a,a) = 1.d0
c         EVENY(a,a) = 1.d0
c         ODDX(a,a) = 1.d0
c         ODDY(a,a) = 1.d0
         EVEN(a,a) = 1.d0
         ODD(a,a) = 1.d0
      end do

      EVEN(0,0) = cos(sqrt((pluscrex/2)**2+(minanny/2)**2))
      EVEN(5,5) = cos(sqrt((pluscrex/2)**2+(minanny/2)**2))
      EVEN(5,0) = 
     .(sin(sqrt((pluscrex/2)**2+(minanny/2)**2))/
     .sqrt((pluscrex/2)**2+(minanny/2)**2))*
     .(i*(pluscrex/2)-(minanny/2))
      EVEN(0,5) = 
     .(sin(sqrt((pluscrex/2)**2+(minanny/2)**2))/
     .sqrt((pluscrex/2)**2+(minanny/2)**2))*
     .(i*(pluscrex/2)+(minanny/2))

      EVEN(1,1) = cos(sqrt((plusannx/2)**2+(mincrey/2)**2))
      EVEN(4,4) = cos(sqrt((plusannx/2)**2+(mincrey/2)**2))
      EVEN(4,1) =
     .(sin(sqrt((plusannx/2)**2+(mincrey/2)**2))/
     .sqrt((plusannx/2)**2+(mincrey/2)**2))*
     .(i*(plusannx/2)-(mincrey/2))
      EVEN(1,4) =
     .(sin(sqrt((plusannx/2)**2+(mincrey/2)**2))/
     .sqrt((plusannx/2)**2+(mincrey/2)**2))*
     .(i*(plusannx/2)+(mincrey/2))

      EVEN(2,2) = cos(sqrt((sqrt(3.d0)*pluscrex/2)**2+
     .            (sqrt(3.d0)*minanny/2)**2))
      EVEN(7,7) = cos(sqrt((sqrt(3.d0)*pluscrex/2)**2+
     .            (sqrt(3.d0)*minanny/2)**2))
      EVEN(7,2) =
     .(sin(sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))/
     .sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))*
     .(i*(sqrt(3.d0)*pluscrex/2)-(sqrt(3.d0)*minanny/2))
      EVEN(2,7) =
     .(sin(sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))/
     .sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))*
     .(i*(sqrt(3.d0)*pluscrex/2)+(sqrt(3.d0)*minanny/2))

      EVEN(3,3) = cos(sqrt((sqrt(3.d0)*plusannx/2)**2+
     .            (sqrt(3.d0)*mincrey/2)**2))
      EVEN(6,6) = cos(sqrt((sqrt(3.d0)*plusannx/2)**2+
     .            (sqrt(3.d0)*mincrey/2)**2))
      EVEN(6,3) =
     .(sin(sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))/
     .sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))*
     .(i*(sqrt(3.d0)*plusannx/2)-(sqrt(3.d0)*mincrey/2))
      EVEN(3,6) =
     .(sin(sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))/
     .sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))*
     .(i*(sqrt(3.d0)*plusannx/2)+(sqrt(3.d0)*mincrey/2))
    

      ODD(1,1) = cos(sqrt((sqrt(2.d0)*pluscrex)**2+
     .            (sqrt(2.d0)*minanny)**2))
      ODD(6,6) = cos(sqrt((sqrt(2.d0)*pluscrex)**2+
     .            (sqrt(2.d0)*minanny)**2))
      ODD(6,1) =
     .(sin(sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))/
     .sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))*
     .(i*(sqrt(2.d0)*pluscrex)-(sqrt(2.d0)*minanny))
      ODD(1,6) =
     .(sin(sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))/
     .sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))*
     .(i*(sqrt(2.d0)*pluscrex)+(sqrt(2.d0)*minanny))


      ODD(2,2) = cos(sqrt((sqrt(2.d0)*plusannx)**2+
     .            (sqrt(2.d0)*mincrey)**2))
      ODD(5,5) = cos(sqrt((sqrt(2.d0)*plusannx)**2+
     .            (sqrt(2.d0)*mincrey)**2))
      ODD(5,2) =
     .(sin(sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))/
     .sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))*
     .(i*(sqrt(2.d0)*plusannx)-(sqrt(2.d0)*mincrey))
      ODD(2,5) =
     .(sin(sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))/
     .sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))*
     .(i*(sqrt(2.d0)*plusannx)+(sqrt(2.d0)*mincrey))

c      CALL EXACT_TAYLOR(7,H,H_C)

      CALL spinmode_operators(EVEN,ODD,PSI_IN,
     .                        H_C,j,k,x,y)
      return
      end

      subroutine spinmode_operators(EVEN,ODD,PSI,
     .                              H,j,k,ion,mode)
      IMPLICIT NONE
      INTEGER a,b,c,j,k,ion,mode,x,y
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,ip
      INTEGER IND(0:7)
      COMPLEX*16 INP(0:7),OUTP(0:7),EVENX(0:7,0:7),PSI(0:2**j*4**k-1)
      COMPLEX*16 EVENY(0:7,0:7),ODDX(0:7,0:7),ODDY(0:7,0:7),H(0:7,0:7)
      COMPLEX*16 NOUTP(0:7),EVEN(0:7,0:7),ODD(0:7,0:7)
      
      mp1 = k*2-mode*2
      mp2 = k*2-mode*2+1
      ip = k*2+j-ion

c---- generate order of bitstrings

c---- generate unique bitstrings using binary
      do ii = 0, 2**(j+2*k-3)-1
      
c------- Create base bitstring (all bits except mp1,mp2,ip)
         call gen_other_bits_fast(j+2*k-1,mp1,mp2,ip,ii,base_bits)
c------- Loop through all 8 possibilities of the 3 bits
c         CALL system_clock(t_tmp1)
         idx = 0
         do a = 0,1
            do b = 0,1
               do c = 0,1
c---------------- Build complete bitstring using bitwise
                  bit_combo = base_bits
                  if (a .EQ. 1) bit_combo = IBSET(bit_combo, ip)
                  if (b .EQ. 1) bit_combo = IBSET(bit_combo, mp2)
                  if (c .EQ. 1) bit_combo = IBSET(bit_combo, mp1)
                  
                  IND(idx) = bit_combo
c                  write (6,*) (ii*8+idx),bit_combo
                  idx = idx + 1
               end do
            end do
         end do
         


c         stop   
         do a = 0,7
            OUTP(a) = PSI(IND(a))
         end do

c------- WORKING HERE: implement the spin x mode matrix
c           Manually implement matrix multiplication
            
 
         NOUTP(0) = EVEN(0,0)*OUTP(0) + EVEN(0,5)*OUTP(5)
         NOUTP(1) = EVEN(1,1)*OUTP(1) + EVEN(1,4)*OUTP(4)
         NOUTP(2) = EVEN(2,2)*OUTP(2) + EVEN(2,7)*OUTP(7)
         NOUTP(3) = EVEN(3,3)*OUTP(3) + EVEN(3,6)*OUTP(6)
         NOUTP(4) = EVEN(4,1)*OUTP(1) + EVEN(4,4)*OUTP(4)
         NOUTP(5) = EVEN(5,0)*OUTP(0) + EVEN(5,5)*OUTP(5)
         NOUTP(6) = EVEN(6,3)*OUTP(3) + EVEN(6,6)*OUTP(6)
         NOUTP(7) = EVEN(7,2)*OUTP(2) + EVEN(7,7)*OUTP(7)
         
         OUTP = NOUTP   

c         NOUTP(0) = EVENY(0,0)*OUTP(0) + EVENY(0,5)*OUTP(5)
c         NOUTP(1) = EVENY(1,1)*OUTP(1) + EVENY(1,4)*OUTP(4)
c         NOUTP(2) = EVENY(2,2)*OUTP(2) + EVENY(2,7)*OUTP(7)
c         NOUTP(3) = EVENY(3,3)*OUTP(3) + EVENY(3,6)*OUTP(6)
c         NOUTP(4) = EVENY(4,1)*OUTP(1) + EVENY(4,4)*OUTP(4)
c         NOUTP(5) = EVENY(5,0)*OUTP(0) + EVENY(5,5)*OUTP(5)
c         NOUTP(6) = EVENY(6,3)*OUTP(3) + EVENY(6,6)*OUTP(6)
c         NOUTP(7) = EVENY(7,2)*OUTP(2) + EVENY(7,7)*OUTP(7)
   
c         OUTP = NOUTP
   
         NOUTP(1) = ODD(1,1)*OUTP(1) + ODD(1,6)*OUTP(6)
         NOUTP(2) = ODD(2,2)*OUTP(2) + ODD(2,5)*OUTP(5)
         NOUTP(5) = ODD(5,2)*OUTP(2) + ODD(5,5)*OUTP(5)
         NOUTP(6) = ODD(6,1)*OUTP(1) + ODD(6,6)*OUTP(6)

         OUTP = NOUTP

         NOUTP(0) = EVEN(0,0)*OUTP(0) + EVEN(0,5)*OUTP(5)
         NOUTP(1) = EVEN(1,1)*OUTP(1) + EVEN(1,4)*OUTP(4)
         NOUTP(2) = EVEN(2,2)*OUTP(2) + EVEN(2,7)*OUTP(7)
         NOUTP(3) = EVEN(3,3)*OUTP(3) + EVEN(3,6)*OUTP(6)
         NOUTP(4) = EVEN(4,1)*OUTP(1) + EVEN(4,4)*OUTP(4)
         NOUTP(5) = EVEN(5,0)*OUTP(0) + EVEN(5,5)*OUTP(5)
         NOUTP(6) = EVEN(6,3)*OUTP(3) + EVEN(6,6)*OUTP(6)
         NOUTP(7) = EVEN(7,2)*OUTP(2) + EVEN(7,7)*OUTP(7)
         
         OUTP = NOUTP 
         
c         NOUTP(1) = ODDY(1,1)*OUTP(1) + ODDY(1,6)*OUTP(6)
c         NOUTP(2) = ODDY(2,2)*OUTP(2) + ODDY(2,5)*OUTP(5)
c         NOUTP(5) = ODDY(5,2)*OUTP(2) + ODDY(5,5)*OUTP(5)
c         NOUTP(6) = ODDY(6,1)*OUTP(1) + ODDY(6,6)*OUTP(6)
         
c         OUTP = NOUTP 
            
c         do a = 0,7
c            OUTP(a) = PSI(IND(a))
c         end do
         
c         OUTP = MATMUL(H,OUTP)

c         OUTP = MATMUL(EVENX,OUTP)   
      
c         OUTP = MATMUL(EVENY,OUTP) 

c         OUTP = MATMUL(ODDX,OUTP)

c         OUTP = MATMUL(ODDY,OUTP)

         do a = 0,7
            PSI(IND(a)) = OUTP(a)
         end do
 
      end do
c      stop 
       
      return
      end

c==== generate other bits in bitstring
      subroutine gen_other_bits_fast(nbits,mp1,mp2,ip,num,result)
      IMPLICIT NONE
      INTEGER nbits,mp1,mp2,ip,num,result
      INTEGER bit_pos,source_bit,temp
      
      result = 0
      source_bit = 0
      temp = num
      
      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2 
     .       .AND. bit_pos .NE. ip) then
            if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
         end if
      end do
      
      return
      end

C====================================================================
C           FIRST ORDER MAGNUS, HIGHER ORDER TAYLOR SERIES
C====================================================================

      subroutine spin_mode_mode(j,k,x,y,z,t2,t1,PSI_IN,
     .                          LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y,z
      DOUBLE PRECISION t2,t1,LD(j,k),RABI(7),w(7),w_k(7),wqbt(7),wj
      COMPLEX*16 pac,pca,paa,pcc
      COMPLEX*16 i,H(0:31,0:31),H_C(0:31,0:31),PSI_IN(0:2**j*4**k-1)

      i = (0.d0,1.d0)
      wj = w(x) - wqbt(x)

      pac = RABI(x)*LD(x,y)*LD(x,z)*
     .     (sin((wj - w_k(z) + w_k(y))*t2)
     .     + i*cos((wj - w_k(z) + w_k(y))*t2)
     .     - sin((wj - w_k(z) + w_k(y))*t1)
     .     - i*cos((wj - w_k(z) + w_k(y))*t1))
     .     / (wj - w_k(z) + w_k(y))

      pca = RABI(x)*LD(x,y)*LD(x,z)*
     .     (sin((wj + w_k(z) - w_k(y))*t2)
     .     + i*cos((wj + w_k(z) - w_k(y))*t2)
     .     - sin((wj + w_k(z) - w_k(y))*t1)
     .     - i*cos((wj + w_k(z) - w_k(y))*t1))
     .     / (wj + w_k(z) - w_k(y))

      paa = RABI(x)*LD(x,y)*LD(x,z)*
     .     (sin((wj + w_k(z) + w_k(y))*t2)
     .     + i*cos((wj + w_k(z) + w_k(y))*t2)
     .     - sin((wj + w_k(z) + w_k(y))*t1)
     .     - i*cos((wj + w_k(z) + w_k(y))*t1))
     .     / (wj + w_k(z) + w_k(y))

      pcc = RABI(x)*LD(x,y)*LD(x,z)*
     .     (sin((-wj + w_k(z) + w_k(y))*t2)
     .     - i*cos((-wj + w_k(z) + w_k(y))*t2)
     .     - sin((-wj + w_k(z) + w_k(y))*t1)
     .     + i*cos((-wj + w_k(z) + w_k(y))*t1))
     .     / (-wj + w_k(z) + w_k(y))


      H = 0.d0
C=======================================================================
C     top-right block: rows 0:15, cols 16:31
C=======================================================================

      H(0,21) = pcc

      H(1,20) = pca
      H(1,22) = sqrt(2.d0)*pcc

      H(2,21) = sqrt(2.d0)*pca
      H(2,23) = sqrt(3.d0)*pcc

      H(3,22) = sqrt(3.d0)*pca

      H(4,17) = pac
      H(4,25) = sqrt(2.d0)*pcc

      H(5,16) = paa
      H(5,18) = sqrt(2.d0)*pac
      H(5,24) = sqrt(2.d0)*pca
      H(5,26) = sqrt(4.d0)*pcc

      H(6,17) = sqrt(2.d0)*paa
      H(6,19) = sqrt(3.d0)*pac
      H(6,25) = sqrt(4.d0)*pca
      H(6,27) = sqrt(6.d0)*pcc

      H(7,18) = sqrt(3.d0)*paa
      H(7,26) = sqrt(6.d0)*pca

      H(8,21) = sqrt(2.d0)*pac
      H(8,29) = sqrt(3.d0)*pcc

      H(9,20) = sqrt(2.d0)*paa
      H(9,22) = sqrt(4.d0)*pac
      H(9,28) = sqrt(3.d0)*pca
      H(9,30) = sqrt(6.d0)*pcc

      H(10,21) = sqrt(4.d0)*paa
      H(10,23) = sqrt(6.d0)*pac
      H(10,29) = sqrt(6.d0)*pca
      H(10,31) = sqrt(9.d0)*pcc

      H(11,22) = sqrt(6.d0)*paa
      H(11,30) = sqrt(9.d0)*pca

      H(12,25) = sqrt(3.d0)*pac

      H(13,24) = sqrt(3.d0)*paa
      H(13,26) = sqrt(6.d0)*pac

      H(14,25) = sqrt(6.d0)*paa
      H(14,27) = sqrt(9.d0)*pac

      H(15,26) = sqrt(9.d0)*paa


C=======================================================================
C     bottom-left block: rows 16:31, cols 0:15
C     Hermitian conjugate of top-right block
C=======================================================================

      H(21,0) = conjg(pcc)

      H(20,1) = conjg(pca)
      H(22,1) = sqrt(2.d0)*conjg(pcc)

      H(21,2) = sqrt(2.d0)*conjg(pca)
      H(23,2) = sqrt(3.d0)*conjg(pcc)

      H(22,3) = sqrt(3.d0)*conjg(pca)

      H(17,4) = conjg(pac)
      H(25,4) = sqrt(2.d0)*conjg(pcc)

      H(16,5) = conjg(paa)
      H(18,5) = sqrt(2.d0)*conjg(pac)
      H(24,5) = sqrt(2.d0)*conjg(pca)
      H(26,5) = sqrt(4.d0)*conjg(pcc)

      H(17,6) = sqrt(2.d0)*conjg(paa)
      H(19,6) = sqrt(3.d0)*conjg(pac)
      H(25,6) = sqrt(4.d0)*conjg(pca)
      H(27,6) = sqrt(6.d0)*conjg(pcc)

      H(18,7) = sqrt(3.d0)*conjg(paa)
      H(26,7) = sqrt(6.d0)*conjg(pca)

      H(21,8) = sqrt(2.d0)*conjg(pac)
      H(29,8) = sqrt(3.d0)*conjg(pcc)

      H(20,9) = sqrt(2.d0)*conjg(paa)
      H(22,9) = sqrt(4.d0)*conjg(pac)
      H(28,9) = sqrt(3.d0)*conjg(pca)
      H(30,9) = sqrt(6.d0)*conjg(pcc)

      H(21,10) = sqrt(4.d0)*conjg(paa)
      H(23,10) = sqrt(6.d0)*conjg(pac)
      H(29,10) = sqrt(6.d0)*conjg(pca)
      H(31,10) = sqrt(9.d0)*conjg(pcc)

      H(22,11) = sqrt(6.d0)*conjg(paa)
      H(30,11) = sqrt(9.d0)*conjg(pca)

      H(25,12) = sqrt(3.d0)*conjg(pac)

      H(24,13) = sqrt(3.d0)*conjg(paa)
      H(26,13) = sqrt(6.d0)*conjg(pac)

      H(25,14) = sqrt(6.d0)*conjg(paa)
      H(27,14) = sqrt(9.d0)*conjg(pac)

      H(26,15) = sqrt(9.d0)*conjg(paa)

      CALL EXACT_TAYLOR(31,H,H_C)
      CALL smm_operators(H_C,PSI_IN,j,k,x,y,z)

      return
      end

      subroutine smm_operators(H,PSI,j,k,ion,m1,m2)
      IMPLICIT NONE
      INTEGER a,b,c,d,e,f,j,k,ion,m1,m2,x,y
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,mp3,mp4,ip
      INTEGER IND(0:31)
      COMPLEX*16 INP(0:31),OUTP(0:31),H(0:31,0:31),PSI(0:2**j*4**k-1)

      mp4 = k*2-m2*2+1
      mp3 = k*2-m2*2
      mp1 = k*2-m1*2
      mp2 = k*2-m1*2+1
      ip = k*2+j-ion

c      write(6,*) mp1,mp2,mp3,mp4,ip
c      stop

c---- generate unique bitstrings using binary
      do ii = 0, 2**(j+2*k-5)-1

c------- Create base bitstring (all bits except mp1,mp2,ip)
        call gen_other_bits_smm(j+2*k-1,mp1,mp2,mp3,mp4,ip,ii,base_bits)
c------- Loop through all 8 possibilities of the 3 bits
c         CALL system_clock(t_tmp1)
         idx = 0
         do a = 0,1
            do b = 0,1
               do c = 0,1
                 do d = 0,1
                   do e = 0,1
c------------------ Build complete bitstring using bitwise
                      bit_combo = base_bits
                      if (a .EQ. 0) bit_combo = IBSET(bit_combo, ip)
                      if (b .EQ. 0) bit_combo = IBSET(bit_combo, mp4)
                      if (c .EQ. 0) bit_combo = IBSET(bit_combo, mp3)
                      if (d .EQ. 0) bit_combo = IBSET(bit_combo, mp2)
                      if (e .EQ. 0) bit_combo = IBSET(bit_combo, mp1)
                      IND(idx) = bit_combo
                      idx = idx + 1
                    end do
                  end do
               end do
            end do
         end do

         OUTP = 0.d0

         do a = 0,31
            INP(a) = PSI(IND(a))
         end do

         OUTP = MATMUL(H,INP)

         do a = 0,31
            PSI(IND(a)) = OUTP(a)
         end do

      end do

      return
      end

      subroutine gen_other_bits_smm(nbits,mp1,mp2,mp3,mp4,
     .                              ip,num,result)
      IMPLICIT NONE
      INTEGER nbits,mp1,mp2,mp3,mp4,ip,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2
     .       .AND. bit_pos .NE. ip .AND. bit_pos .NE. mp3
     .       .AND. bit_pos .NE. mp4) then
            if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
         end if
      end do

      return
      end


C====================================================================
C                      SPIN X SINGLE
C====================================================================

      subroutine spin_single(j,k,x,y,t2,t1,PSI_IN,
     .                       LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y
      DOUBLE PRECISION LD(j,k),RABI(7),w(7),w_k(7),wqbt(7),t2,t1,wj
      COMPLEX*16 paa,pac,pca,pcc
      COMPLEX*16 i,H(0:7,0:7),H_C(0:7,0:7),PSI_IN(0:2**j*4**k-1)

      i = (0.d0,1.d0)
      wj = w(x) - wqbt(x)

       pac = RABI(x)*LD(x,y)*LD(x,y)*
     .     (sin((wj - w_k(y) + w_k(y))*t2)
     .     + i*cos((wj - w_k(y) + w_k(y))*t2)
     .     - sin((wj - w_k(y) + w_k(y))*t1)
     .     - i*cos((wj - w_k(y) + w_k(y))*t1))
     .     / (wj - w_k(y) + w_k(y))

      pca = RABI(x)*LD(x,y)*LD(x,y)*
     .     (sin((wj + w_k(y) - w_k(y))*t2)
     .     + i*cos((wj + w_k(y) - w_k(y))*t2)
     .     - sin((wj + w_k(y) - w_k(y))*t1)
     .     - i*cos((wj + w_k(y) - w_k(y))*t1))
     .     / (wj + w_k(y) - w_k(y))

      paa = RABI(x)*LD(x,y)*LD(x,y)*
     .     (sin((wj + w_k(y) + w_k(y))*t2)
     .     + i*cos((wj + w_k(y) + w_k(y))*t2)
     .     - sin((wj + w_k(y) + w_k(y))*t1)
     .     - i*cos((wj + w_k(y) + w_k(y))*t1))
     .     / (wj + w_k(y) + w_k(y))

      pcc = RABI(x)*LD(x,y)*LD(x,y)*
     .     (sin((-wj + w_k(y) + w_k(y))*t2)
     .     - i*cos((-wj + w_k(y) + w_k(y))*t2)
     .     - sin((-wj + w_k(y) + w_k(y))*t1)
     .     + i*cos((-wj + w_k(y) + w_k(y))*t1))
     .     / (-wj + w_k(y) + w_k(y))

      H = (0.d0,0.d0)

C=======================================================================
C     top-right block: rows 0:3, cols 4:7
C=======================================================================

      H(0,4) = pac
      H(0,6) = sqrt(2.d0)*paa

      H(1,5) = pca + 2.d0*pac
      H(1,7) = sqrt(6.d0)*paa

      H(2,4) = sqrt(2.d0)*pcc
      H(2,6) = 2.d0*pca + 3.d0*pac

      H(3,5) = sqrt(6.d0)*pcc
      H(3,7) = 3.d0*pca


C=======================================================================
C     bottom-left block: rows 4:7, cols 0:3
C     Hermitian conjugate
C=======================================================================

      H(4,0) = conjg(pac)
      H(6,0) = sqrt(2.d0)*conjg(paa)

      H(5,1) = conjg(pca + 2.d0*pac)
      H(7,1) = sqrt(6.d0)*conjg(paa)

      H(4,2) = sqrt(2.d0)*conjg(pcc)
      H(6,2) = conjg(2.d0*pca + 3.d0*pac)

      H(5,3) = sqrt(6.d0)*conjg(pcc)
      H(7,3) = 3.d0*conjg(pca)
      
      call EXACT_TAYLOR(7,H,H_C)
      call ssingle_operators(H_C,PSI_IN,j,k,x,y)
      
      return
      end


      subroutine ssingle_operators(H,PSI,
     .                              j,k,ion,mode)
      IMPLICIT NONE
      INTEGER a,b,c,j,k,ion,mode,x,y
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,ip
      INTEGER IND(0:7)
      COMPLEX*16 INP(0:7),OUTP(0:7),PSI(0:2**j*4**k-1)
      COMPLEX*16 H(0:7,0:7)
      
      mp1 = k*2-mode*2
      mp2 = k*2-mode*2+1
      ip = k*2+j-ion

c---- generate order of bitstrings

c---- generate unique bitstrings using binary
      do ii = 0, 2**(j+2*k-3)-1
      
c------- Create base bitstring (all bits except mp1,mp2,ip)
         call gen_other_bits_fast(j+2*k-1,mp1,mp2,ip,ii,base_bits)
c------- Loop through all 8 possibilities of the 3 bits
c         CALL system_clock(t_tmp1)
         idx = 0
         do a = 0,1
            do b = 0,1
               do c = 0,1
c---------------- Build complete bitstring using bitwise
                  bit_combo = base_bits
                  if (a .EQ. 1) bit_combo = IBSET(bit_combo, ip)
                  if (b .EQ. 1) bit_combo = IBSET(bit_combo, mp2)
                  if (c .EQ. 1) bit_combo = IBSET(bit_combo, mp1)
                  
                  IND(idx) = bit_combo
c                  write (6,*) (ii*8+idx),bit_combo
                  idx = idx + 1
               end do
            end do
         end do
         


c         stop   
         do a = 0,7
            OUTP(a) = PSI(IND(a))
         end do

         OUTP = MATMUL(H,OUTP)

         do a = 0,7
            PSI(IND(a)) = OUTP(a)
         end do
 
      end do
c      stop 
       
      return
      end



C====================================================================
C                 SECOND ORDER MAGNUS IN LD REGIME
C====================================================================           

C====================================================================
C                       APPLY SPIN X SPIN
C====================================================================
      subroutine spin_spin(j,k,x,y,t3,t2,t1,PSI_IN,RABI,
     .                     w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y
      DOUBLE PRECISION t1,t2,t3
      DOUBLE PRECISION RABI(7),w(7),wqbt(7),w_k(7),wj,wjp
      DOUBLE PRECISION ppx,ppy,pmx,pmy
      COMPLEX*16 PSI_IN(0:2**j*4**k-1),pp,pm
      COMPLEX*16 H(0:3,0:3),H_C(0:3,0:3),i
      
      i = (0.d0,1.d0)
      wj = w(x) - wqbt(x)
      wjp = w(y) - wqbt(y)
      
      ppx = 5.d-1*(RABI(x)*RABI(y)/(wj*wjp))*
     . (-cos(wjp*t3 + wj*t2)
     .  + cos(wjp*t3 + wj*t1)
     .  + cos(wj*t3 + wjp*t2)
     .  - cos(wj*t3 + wjp*t1)
     .  - cos(wjp*t2 + wj*t1)
     .  + cos(wj*t2 + wjp*t1))

      ppy = 5.d-1*(RABI(x)*RABI(y)/(wj*wjp))*
     .  (sin(wjp*t3 + wj*t2)
     .  - sin(wjp*t3 + wj*t1)
     .  - sin(wj*t3 + wjp*t2)
     .  + sin(wj*t3 + wjp*t1)
     .  + sin(wjp*t2 + wj*t1)
     .  - sin(wj*t2 + wjp*t1))

      pmx = 5.d-1*(RABI(x)*RABI(y)/(wj*wjp))*
     . (-cos(wjp*t3 - wj*t2)
     .  + cos(wjp*t3 - wj*t1)
     .  + cos(wj*t3 - wjp*t2)
     .  - cos(wj*t3 - wjp*t1)
     .  - cos(wjp*t2 - wj*t1)
     .  + cos(wj*t2 - wjp*t1))

      pmy = 5.d-1*(RABI(x)*RABI(y)/(wj*wjp))*
     . (-sin(wjp*t3 - wj*t2)
     .  + sin(wjp*t3 - wj*t1)
     .  - sin(wj*t3 - wjp*t2)
     .  + sin(wj*t3 - wjp*t1)
     .  - sin(wjp*t2 - wj*t1)
     .  - sin(wj*t2 - wjp*t1))
      
      pp = ppx + i*ppy
      pm = pmx + i*pmy

      H = 0.d0
      H(0,3) = pp
      H(3,0) = conjg(pp)
      H(1,2) = pm
      H(2,1) = conjg(pm)

      call exact_taylor(3,H,H_C)

      call ss_operators(H_C,PSI_IN,j,k,x,y)

      return
      end 


      subroutine ss_operators(H,PSI,j,k,ion1,ion2)
      IMPLICIT NONE
      INTEGER a,b,c,j,k,ion1,ion2
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER ip1,ip2
      INTEGER IND(0:3)
      COMPLEX*16 OUTP(0:3),PSI(0:2**j*4**k-1)
      COMPLEX*16 H(0:3,0:3)


      ip1 = j+k*2-ion1
      ip2 = j+k*2-ion2

      do ii = 0,2**(j+2*k-2)-1
            call gen_other_bits_ss(j+2*k-1,ip1,ip2,ii,base_bits)
            idx = 0
          do a = 0,1
            do b = 0,1
               bit_combo = base_bits

               if (a .eq. 1) bit_combo = IBSET(bit_combo, ip1)
               if (b .eq. 1) bit_combo = IBSET(bit_combo, ip2)

               IND(idx) = bit_combo
               idx = idx + 1
            end do
         end do
      
      
         do a = 0,3
            OUTP(a) = PSI(IND(a))
         end do

         OUTP = MATMUL(H,OUTP)

         do a = 0,3
            PSI(IND(a)) = OUTP(a)   
         end do

      end do
            

      return
      end


      subroutine gen_other_bits_ss(nbits,ip1,ip2,num,result)
      IMPLICIT NONE
      INTEGER nbits,ip1,ip2,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
         if (bit_pos .NE. ip1 .AND. bit_pos .NE. ip2)
     .      then
            if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
          end if
      end do

      return
      end

C=====================================================================
C                             SIGMA Z
C=====================================================================


      subroutine sigma_z(j,k,x,t2,t1,PSI_IN,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x
      DOUBLE PRECISION t2,t1,RABI(7),w(7),wqbt(7),w_k(7),wj
      COMPLEX*16 i,pm,H(0:1,0:1),H_C(0:1,0:1),PSI_IN(0:2**j*4**k-1)

      i = (0.d0,1.d0)
      wj = wqbt(x) - w(x)

      pm = -5.d-1*i*RABI(x)*RABI(x)*
     .     (2.d0*i*(sin(wj*t2 - wj*t1) - wj*t2 + wj*t1))
     .     / (wj**2)


      H_C = 0.d0
      H_C(0,0) = exp(-i*pm)
      H_C(1,1) = exp(i*pm)
     
      CALL carrier_operator(H_C,PSI_IN,j,k,x)

      return
      end


C=====================================================================
c                       APPLY SPIN X SPIN X MODE
C=====================================================================
      subroutine spin_spin_mode(j,k,x,y,z,t3,t2,t1,PSI_IN,
     .                          LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y,z,a
      DOUBLE PRECISION t1,t2,t3
      DOUBLE PRECISION LD(j,k),RABI(7),w(7),wqbt(7),w_k(7),wj,wjp
      DOUBLE PRECISION ppannx,ppanny,ppcrex,ppcrey
      DOUBLE PRECISION pmannx,pmanny,pmcrex,pmcrey
      DOUBLE PRECISION pannpx,pannpy,pcrepx,pcrepy
      DOUBLE PRECISION pannmx,pannmy,pcremx,pcremy
      COMPLEX*16 PSI_IN(0:2**j*4**k-1)
      COMPLEX*16 EVEN(0:15,0:15),ODD(0:15,0:15)
      COMPLEX*16 i
      
      i = (0.d0,1.d0)
      
      wj = w(x) - wqbt(x)
      wjp = w(y) - wqbt(y)

      ppanny = 5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(wjp+w_k(z))))*
     . (-cos((w_k(z)+wjp)*t3 + wj*t2)
     .  + cos((w_k(z)+wjp)*t3 + wj*t1)
     .  + cos(wj*t3 + (w_k(z)+wjp)*t2)
     .  - cos(wj*t3 + (w_k(z)+wjp)*t1)
     .  - cos((w_k(z)+wjp)*t2 + wj*t1)
     .  + cos(wj*t2 + (w_k(z)+wjp)*t1))

      ppannx = 5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(wjp+w_k(z))))*
     .  (sin((w_k(z)+wjp)*t3 + wj*t2)
     .  - sin((w_k(z)+wjp)*t3 + wj*t1)
     .  - sin(wj*t3 + (w_k(z)+wjp)*t2)
     .  + sin(wj*t3 + (w_k(z)+wjp)*t1)
     .  + sin((w_k(z)+wjp)*t2 + wj*t1)
     .  - sin(wj*t2 + (w_k(z)+wjp)*t1))

      ppcrey = -5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(w_k(z)-wjp)))*
     . (-cos((wjp-w_k(z))*t3 + wj*t2)
     .  + cos((wjp-w_k(z))*t3 + wj*t1)
     .  + cos(wj*t3 + (wjp-w_k(z))*t2)
     .  - cos(wj*t3 + (wjp-w_k(z))*t1)
     .  - cos((wjp-w_k(z))*t2 + wj*t1)
     .  + cos(wj*t2 + (wjp-w_k(z))*t1))

      ppcrex = -5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(w_k(z)-wjp)))*
     .  (sin((wjp-w_k(z))*t3 + wj*t2)
     .  - sin((wjp-w_k(z))*t3 + wj*t1)
     .  - sin(wj*t3 + (wjp-w_k(z))*t2)
     .  + sin(wj*t3 + (wjp-w_k(z))*t1)
     .  + sin((wjp-w_k(z))*t2 + wj*t1)
     .  - sin(wj*t2 + (wjp-w_k(z))*t1))

      pmanny = 5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(w_k(z)-wjp)))*
     .  (cos((wjp-w_k(z))*t3 - wj*t2)
     .  - cos((wjp-w_k(z))*t3 - wj*t1)
     .  - cos(wj*t3 + (w_k(z)-wjp)*t2)
     .  + cos(wj*t3 + (w_k(z)-wjp)*t1)
     .  + cos((wjp-w_k(z))*t2 - wj*t1)
     .  - cos(wj*t2 + (w_k(z)-wjp)*t1)
     .  + cos((w_k(z)-wjp+wj)*t2)
     .  - cos((-w_k(z)+wjp-wj)*t2))

      pmannx = 5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(w_k(z)-wjp)))*
     .  (sin((wjp-w_k(z))*t3 - wj*t2)
     .  - sin((wjp-w_k(z))*t3 - wj*t1)
     .  + sin(wj*t3 + (w_k(z)-wjp)*t2)
     .  - sin(wj*t3 + (w_k(z)-wjp)*t1)
     .  + sin((wjp-w_k(z))*t2 - wj*t1)
     .  + sin(wj*t2 + (w_k(z)-wjp)*t1)
     .  - sin((w_k(z)-wjp+wj)*t2)
     .  - sin((-w_k(z)+wjp-wj)*t2))

      pmcrey = 5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(wjp+w_k(z))))*
     . (-cos((w_k(z)+wjp)*t3 - wj*t2)
     .  + cos((w_k(z)+wjp)*t3 - wj*t1)
     .  + cos(wj*t3 + (-w_k(z)-wjp)*t2)
     .  - cos(wj*t3 + (-w_k(z)-wjp)*t1)
     .  - cos((w_k(z)+wjp)*t2 - wj*t1)
     .  + cos(wj*t2 + (-w_k(z)-wjp)*t1)
     .  + cos((w_k(z)+wjp-wj)*t2)
     .  - cos((-w_k(z)-wjp+wj)*t2))

      pmcrex = 5.d-1*(RABI(x)*RABI(y)*LD(y,z)/
     .        (wj*(wjp+w_k(z))))*
     . (-sin((w_k(z)+wjp)*t3 - wj*t2)
     .  + sin((w_k(z)+wjp)*t3 - wj*t1)
     .  - sin(wj*t3 + (-w_k(z)-wjp)*t2)
     .  + sin(wj*t3 + (-w_k(z)-wjp)*t1)
     .  - sin((w_k(z)+wjp)*t2 - wj*t1)
     .  - sin(wj*t2 + (-w_k(z)-wjp)*t1)
     .  + sin((w_k(z)+wjp-wj)*t2)
     .  + sin((-w_k(z)-wjp+wj)*t2))

      pannpx = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(sin(wjp*t3 + (w_k(y)+wj)*t2)
     .       - sin(wjp*t3 + (w_k(y)+wj)*t1)
     .       - sin((w_k(y)+wj)*t3 + wjp*t2)
     .       + sin((w_k(y)+wj)*t3 + wjp*t1)
     .       + sin(wjp*t2 + (w_k(y)+wj)*t1)
     .       - sin((w_k(y)+wj)*t2 + wjp*t1))
     .       / ((w_k(y)+wj)*wjp)

      pannpy = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(cos(wjp*t3 + (w_k(y)+wj)*t2)
     .       - cos(wjp*t3 + (w_k(y)+wj)*t1)
     .       - cos((w_k(y)+wj)*t3 + wjp*t2)
     .       + cos((w_k(y)+wj)*t3 + wjp*t1)
     .       + cos(wjp*t2 + (w_k(y)+wj)*t1)
     .       - cos((w_k(y)+wj)*t2 + wjp*t1))
     .       / ((w_k(y)+wj)*wjp)

      pannmx = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(sin(wjp*t3 - (w_k(y)+wj)*t2)
     .       - sin(wjp*t3 - (w_k(y)+wj)*t1)
     .       + sin((w_k(y)+wj)*t3 - wjp*t2)
     .       - sin((w_k(y)+wj)*t3 - wjp*t1)
     .       + sin(wjp*t2 - (w_k(y)+wj)*t1)
     .       + sin((w_k(y)+wj)*t2 - wjp*t1))
     .       / ((w_k(y)+wj)*wjp)

      pannmy = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(-cos(wjp*t3 - (w_k(y)+wj)*t2)
     .       + cos(wjp*t3 - (w_k(y)+wj)*t1)
     .       + cos((w_k(y)+wj)*t3 - wjp*t2)
     .       - cos((w_k(y)+wj)*t3 - wjp*t1)
     .       - cos(wjp*t2 - (w_k(y)+wj)*t1)
     .       + cos((w_k(y)+wj)*t2 - wjp*t1))
     .       / ((w_k(y)+wj)*wjp)

      pcrepx = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(-sin(wjp*t3 + (wj-w_k(y))*t2)
     .       + sin(wjp*t3 + (wj-w_k(y))*t1)
     .       + sin((wj-w_k(y))*t3 + wjp*t2)
     .       - sin((wj-w_k(y))*t3 + wjp*t1)
     .       - sin(wjp*t2 + (wj-w_k(y))*t1)
     .       + sin((wj-w_k(y))*t2 + wjp*t1))
     .       / ((w_k(y)-wj)*wjp)

      pcrepy = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(-cos(wjp*t3 + (wj-w_k(y))*t2)
     .       + cos(wjp*t3 + (wj-w_k(y))*t1)
     .       + cos((wj-w_k(y))*t3 + wjp*t2)
     .       - cos((wj-w_k(y))*t3 + wjp*t1)
     .       - cos(wjp*t2 + (wj-w_k(y))*t1)
     .       + cos((wj-w_k(y))*t2 + wjp*t1))
     .       / ((w_k(y)-wj)*wjp)

      pcremx = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(-sin(wjp*t3 + (w_k(y)-wj)*t2)
     .       + sin(wjp*t3 + (w_k(y)-wj)*t1)
     .       - sin((wj-w_k(y))*t3 - wjp*t2)
     .       + sin((wj-w_k(y))*t3 - wjp*t1)
     .       - sin(wjp*t2 + (w_k(y)-wj)*t1)
     .       - sin((wj-w_k(y))*t2 - wjp*t1))
     .       / ((w_k(y)-wj)*wjp)

      pcremy = 5.d-1*RABI(x)*RABI(y)*LD(x,z)
     .       *(cos(wjp*t3 + (w_k(y)-wj)*t2)
     .       - cos(wjp*t3 + (w_k(y)-wj)*t1)
     .       - cos((wj-w_k(y))*t3 - wjp*t2)
     .       + cos((wj-w_k(y))*t3 - wjp*t1)
     .       + cos(wjp*t2 + (w_k(y)-wj)*t1)
     .       - cos((wj-w_k(y))*t2 - wjp*t1))
     .       / ((w_k(y)-wj)*wjp)

      ppcrex = ppcrex + pcrepx
      ppcrey = ppcrey + pcrepy

      ppannx = ppannx + pannpx
      ppanny = ppanny + pannpy

      pmcrex = pmcrex + pcremx
      pmcrey = pmcrey + pcremy

      pmannx = pmannx + pannmx
      pmanny = pmanny + pannmy 

      EVEN = 0.d0
      ODD = 0.d0

      do a = 0,15
         EVEN(a,a) = 1.d0
         ODD(a,a) = 1.d0
      end do
            
      EVEN(0,0) = cos(sqrt((ppcrex/2)**2+(ppcrey/2)**2))
      EVEN(13,13) = cos(sqrt((ppcrex/2)**2+(ppcrey/2)**2))
      EVEN(13,0) =
     .(sin(sqrt((ppcrex/2)**2+(ppcrey/2)**2))/
     .sqrt((ppcrex/2)**2+(ppcrey/2)**2))*
     .(i*(ppcrex/2)-(ppcrey/2))
      EVEN(0,13) =
     .(sin(sqrt((ppcrex/2)**2+(ppcrey/2)**2))/
     .sqrt((ppcrex/2)**2+(ppcrey/2)**2))*
     .(i*(ppcrex/2)+(ppcrey/2))
      
      EVEN(4,4) = cos(sqrt((pmcrex/2)**2+(pmcrey/2)**2))
      EVEN(9,9) = cos(sqrt((pmcrex/2)**2+(pmcrey/2)**2))
      EVEN(9,4) =
     .(sin(sqrt((pmcrex/2)**2+(pmcrey/2)**2))/
     .sqrt((pmcrex/2)**2+(pmcrey/2)**2))*
     .(i*(pmcrex/2)-(pmcrey/2))
      EVEN(4,9) =
     .(sin(sqrt((pmcrex/2)**2+(pmcrey/2)**2))/
     .sqrt((pmcrex/2)**2+(pmcrey/2)**2))*
     .(i*(pmcrex/2)+(pmcrey/2))

      EVEN(1,1) = cos(sqrt((ppannx/2)**2+(ppanny/2)**2))
      EVEN(12,12) = cos(sqrt((ppannx/2)**2+(ppanny/2)**2))
      EVEN(12,1) =
     .(sin(sqrt((ppannx/2)**2+(ppanny/2)**2))/
     .sqrt((ppannx/2)**2+(ppanny/2)**2))*
     .(i*(ppannx/2)-(ppanny/2))
      EVEN(1,12) =
     .(sin(sqrt((ppannx/2)**2+(ppanny/2)**2))/
     .sqrt((ppannx/2)**2+(ppanny/2)**2))*
     .(i*(ppannx/2)+(ppanny/2))

      EVEN(5,5) = cos(sqrt((pmannx/2)**2+(pmanny/2)**2))
      EVEN(8,8) = cos(sqrt((pmannx/2)**2+(pmanny/2)**2))
      EVEN(8,5) =
     .(sin(sqrt((pmannx/2)**2+(pmanny/2)**2))/
     .sqrt((pmannx/2)**2+(pmanny/2)**2))*
     .(i*(pmannx/2)-(pmanny/2))
      EVEN(5,8) =
     .(sin(sqrt((pmannx/2)**2+(pmanny/2)**2))/
     .sqrt((pmannx/2)**2+(pmanny/2)**2))*
     .(i*(pmannx/2)+(pmanny/2))

      EVEN(2,2) = cos(sqrt((sqrt(3.d0)*ppcrex/2)**2+
     .            (sqrt(3.d0)*ppcrey/2)**2))
      EVEN(15,15) = cos(sqrt((sqrt(3.d0)*ppcrex/2)**2+
     .            (sqrt(3.d0)*ppcrey/2)**2))
      EVEN(15,2) =
     .(sin(sqrt((sqrt(3.d0)*ppcrex/2)**2+(sqrt(3.d0)*ppcrey/2)**2))/
     .sqrt((sqrt(3.d0)*ppcrex/2)**2+(sqrt(3.d0)*ppcrey/2)**2))*
     .(i*(sqrt(3.d0)*ppcrex/2)-(sqrt(3.d0)*ppcrey/2))
      EVEN(2,15) =
     .(sin(sqrt((sqrt(3.d0)*ppcrex/2)**2+(sqrt(3.d0)*ppcrey/2)**2))/
     .sqrt((sqrt(3.d0)*ppcrex/2)**2+(sqrt(3.d0)*ppcrey/2)**2))*
     .(i*(sqrt(3.d0)*ppcrex/2)+(sqrt(3.d0)*ppcrey/2))

      EVEN(6,6) = cos(sqrt((sqrt(3.d0)*pmcrex/2)**2+
     .            (sqrt(3.d0)*pmcrey/2)**2))
      EVEN(11,11) = cos(sqrt((sqrt(3.d0)*pmcrex/2)**2+
     .            (sqrt(3.d0)*pmcrey/2)**2))
      EVEN(11,6) =
     .(sin(sqrt((sqrt(3.d0)*pmcrex/2)**2+(sqrt(3.d0)*pmcrey/2)**2))/
     .sqrt((sqrt(3.d0)*pmcrex/2)**2+(sqrt(3.d0)*pmcrey/2)**2))*
     .(i*(sqrt(3.d0)*pmcrex/2)-(sqrt(3.d0)*pmcrey/2))
      EVEN(6,11) =
     .(sin(sqrt((sqrt(3.d0)*pmcrex/2)**2+(sqrt(3.d0)*pmcrey/2)**2))/
     .sqrt((sqrt(3.d0)*pmcrex/2)**2+(sqrt(3.d0)*pmcrey/2)**2))*
     .(i*(sqrt(3.d0)*pmcrex/2)+(sqrt(3.d0)*pmcrey/2))

      EVEN(3,3) = cos(sqrt((sqrt(3.d0)*ppannx/2)**2+
     .            (sqrt(3.d0)*ppanny/2)**2))
      EVEN(14,14) = cos(sqrt((sqrt(3.d0)*ppannx/2)**2+
     .            (sqrt(3.d0)*ppanny/2)**2))
      EVEN(14,3) =
     .(sin(sqrt((sqrt(3.d0)*ppannx/2)**2+(sqrt(3.d0)*ppanny/2)**2))/
     .sqrt((sqrt(3.d0)*ppannx/2)**2+(sqrt(3.d0)*ppanny/2)**2))*
     .(i*(sqrt(3.d0)*ppannx/2)-(sqrt(3.d0)*ppanny/2))
      EVEN(3,14) =
     .(sin(sqrt((sqrt(3.d0)*ppannx/2)**2+(sqrt(3.d0)*ppanny/2)**2))/
     .sqrt((sqrt(3.d0)*ppannx/2)**2+(sqrt(3.d0)*ppanny/2)**2))*
     .(i*(sqrt(3.d0)*ppannx/2)+(sqrt(3.d0)*ppanny/2))

      EVEN(7,7) = cos(sqrt((sqrt(3.d0)*pmannx/2)**2+
     .            (sqrt(3.d0)*pmanny/2)**2))
      EVEN(10,10) = cos(sqrt((sqrt(3.d0)*pmannx/2)**2+
     .            (sqrt(3.d0)*pmanny/2)**2))
      EVEN(10,7) =
     .(sin(sqrt((sqrt(3.d0)*pmannx/2)**2+(sqrt(3.d0)*pmanny/2)**2))/
     .sqrt((sqrt(3.d0)*pmannx/2)**2+(sqrt(3.d0)*pmanny/2)**2))*
     .(i*(sqrt(3.d0)*pmannx/2)-(sqrt(3.d0)*pmanny/2))
      EVEN(7,10) =
     .(sin(sqrt((sqrt(3.d0)*pmannx/2)**2+(sqrt(3.d0)*pmanny/2)**2))/
     .sqrt((sqrt(3.d0)*pmannx/2)**2+(sqrt(3.d0)*pmanny/2)**2))*
     .(i*(sqrt(3.d0)*pmannx/2)+(sqrt(3.d0)*pmanny/2))

      


      ODD(1,1) = cos(sqrt((sqrt(2.d0)*ppcrex)**2+
     .            (sqrt(2.d0)*ppcrey)**2))
      ODD(14,14) = cos(sqrt((sqrt(2.d0)*ppcrex)**2+
     .            (sqrt(2.d0)*ppcrey)**2))
      ODD(14,1) =
     .(sin(sqrt((sqrt(2.d0)*ppcrex)**2+(sqrt(2.d0)*ppcrey)**2))/
     .sqrt((sqrt(2.d0)*ppcrex)**2+(sqrt(2.d0)*ppcrey)**2))*
     .(i*(sqrt(2.d0)*ppcrex)-(sqrt(2.d0)*ppcrey))
      ODD(1,14) =
     .(sin(sqrt((sqrt(2.d0)*ppcrex)**2+(sqrt(2.d0)*ppcrey)**2))/
     .sqrt((sqrt(2.d0)*ppcrex)**2+(sqrt(2.d0)*ppcrey)**2))*
     .(i*(sqrt(2.d0)*ppcrex)+(sqrt(2.d0)*ppcrey))


      ODD(5,5) = cos(sqrt((sqrt(2.d0)*pmcrex)**2+
     .            (sqrt(2.d0)*pmcrey)**2))
      ODD(10,10) = cos(sqrt((sqrt(2.d0)*pmcrex)**2+
     .            (sqrt(2.d0)*pmcrey)**2))
      ODD(10,5) =
     .(sin(sqrt((sqrt(2.d0)*pmcrex)**2+(sqrt(2.d0)*pmcrey)**2))/
     .sqrt((sqrt(2.d0)*pmcrex)**2+(sqrt(2.d0)*pmcrey)**2))* 
     .(i*(sqrt(2.d0)*pmcrex)-(sqrt(2.d0)*pmcrey))
      ODD(5,10) =
     .(sin(sqrt((sqrt(2.d0)*pmcrex)**2+(sqrt(2.d0)*pmcrey)**2))/
     .sqrt((sqrt(2.d0)*pmcrex)**2+(sqrt(2.d0)*pmcrey)**2))*
     .(i*(sqrt(2.d0)*pmcrex)+(sqrt(2.d0)*pmcrey))


      ODD(2,2) = cos(sqrt((sqrt(2.d0)*ppannx)**2+
     .            (sqrt(2.d0)*ppanny)**2))
      ODD(13,13) = cos(sqrt((sqrt(2.d0)*ppannx)**2+
     .            (sqrt(2.d0)*ppanny)**2))
      ODD(13,2) =
     .(sin(sqrt((sqrt(2.d0)*ppannx)**2+(sqrt(2.d0)*ppanny)**2))/
     .sqrt((sqrt(2.d0)*ppannx)**2+(sqrt(2.d0)*ppanny)**2))*
     .(i*(sqrt(2.d0)*ppannx)-(sqrt(2.d0)*ppanny))
      ODD(2,13) =
     .(sin(sqrt((sqrt(2.d0)*ppannx)**2+(sqrt(2.d0)*ppanny)**2))/
     .sqrt((sqrt(2.d0)*ppannx)**2+(sqrt(2.d0)*ppanny)**2))*
     .(i*(sqrt(2.d0)*ppannx)+(sqrt(2.d0)*ppanny))

      ODD(6,6) = cos(sqrt((sqrt(2.d0)*pmannx)**2+
     .            (sqrt(2.d0)*pmanny)**2))
      ODD(9,9) = cos(sqrt((sqrt(2.d0)*pmannx)**2+
     .            (sqrt(2.d0)*pmanny)**2))
      ODD(9,6) =
     .(sin(sqrt((sqrt(2.d0)*pmannx)**2+(sqrt(2.d0)*pmanny)**2))/
     .sqrt((sqrt(2.d0)*pmannx)**2+(sqrt(2.d0)*pmanny)**2))*
     .(i*(sqrt(2.d0)*pmannx)-(sqrt(2.d0)*pmanny))
      ODD(6,9) =
     .(sin(sqrt((sqrt(2.d0)*pmannx)**2+(sqrt(2.d0)*pmanny)**2))/
     .sqrt((sqrt(2.d0)*pmannx)**2+(sqrt(2.d0)*pmanny)**2))*
     .(i*(sqrt(2.d0)*pmannx)+(sqrt(2.d0)*pmanny))

     
      call ssm_operators(EVEN,ODD,PSI_IN,j,k,x,y,z)

      return
      end

      subroutine ssm_operators(EVEN,ODD,PSI,j,k,ion1,ion2,mode)
      IMPLICIT NONE
      INTEGER a,b,c,d,j,k,ion1,ion2,mode,x,y
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,ip1,ip2
      INTEGER IND(0:15)
      COMPLEX*16 OUTP(0:15),EVEN(0:15,0:15),ODD(0:15,0:15),NOUTP(0:15)
      COMPLEX*16 PSI(0:2**j*4**k-1)

      mp1 = k*2-mode*2
      mp2 = k*2-mode*2+1
      ip1 = k*2+j-ion1
      ip2 = k*2+j-ion2

      do ii = 0, 2**(j+2*k-4)-1
         call gen_other_bits_ssm(j+2*k-1,mp1,mp2,ip1,ip2,ii,base_bits)
         idx = 0
         
         do a = 0,1
            do b = 0,1
              do c = 0,1
                do d = 0,1
                bit_combo = base_bits
                if (a .eq. 1) bit_combo = IBSET(bit_combo, ip1)
                if (b .eq. 1) bit_combo = IBSET(bit_combo, ip2)
                if (c .eq. 1) bit_combo = IBSET(bit_combo, mp2)
                if (d .eq. 1) bit_combo = IBSET(bit_combo, mp1)
            
                IND(idx) = bit_combo
                idx = idx + 1
                end do
              end do
            end do
         end do

        do a = 0,15
          OUTP(a) = PSI(IND(a))
        end do
        
        NOUTP(0) = EVEN(0,0)*OUTP(0) + EVEN(0,13)*OUTP(13)
        NOUTP(1) = EVEN(1,1)*OUTP(1) + EVEN(1,12)*OUTP(12)
        NOUTP(2) = EVEN(2,2)*OUTP(2) + EVEN(2,15)*OUTP(15)
        NOUTP(3) = EVEN(3,3)*OUTP(3) + EVEN(3,14)*OUTP(14)

        NOUTP(4) = EVEN(4,4)*OUTP(4) + EVEN(4,9)*OUTP(9)
        NOUTP(5) = EVEN(5,5)*OUTP(5) + EVEN(5,8)*OUTP(8)
        NOUTP(6) = EVEN(6,6)*OUTP(6) + EVEN(6,11)*OUTP(11)
        NOUTP(7) = EVEN(7,7)*OUTP(7) + EVEN(7,10)*OUTP(10)

        NOUTP(8) = EVEN(8,5)*OUTP(5) + EVEN(8,8)*OUTP(8)
        NOUTP(9) = EVEN(9,4)*OUTP(4) + EVEN(9,9)*OUTP(9)
        NOUTP(10) = EVEN(10,7)*OUTP(7) + EVEN(10,10)*OUTP(10)
        NOUTP(11) = EVEN(11,6)*OUTP(6) + EVEN(11,11)*OUTP(11)

        NOUTP(12) = EVEN(12,1)*OUTP(1) + EVEN(12,12)*OUTP(12)
        NOUTP(13) = EVEN(13,0)*OUTP(0) + EVEN(13,13)*OUTP(13)
        NOUTP(14) = EVEN(14,3)*OUTP(3) + EVEN(14,14)*OUTP(14)
        NOUTP(15) = EVEN(15,2)*OUTP(2) + EVEN(15,15)*OUTP(15) 

        OUTP = NOUTP

        NOUTP(1) = ODD(1,1)*OUTP(1) + ODD(1,14)*OUTP(14)
        NOUTP(2) = ODD(2,2)*OUTP(2) + ODD(2,13)*OUTP(13)
        NOUTP(5) = ODD(5,5)*OUTP(5) + ODD(5,10)*OUTP(10)
        NOUTP(6) = ODD(6,6)*OUTP(6) + ODD(6,9)*OUTP(9)
        NOUTP(9) = ODD(9,6)*OUTP(6) + ODD(9,9)*OUTP(9)
        NOUTP(10) = ODD(10,5)*OUTP(5) + ODD(10,10)*OUTP(10)
        NOUTP(13) = ODD(13,2)*OUTP(2) + ODD(13,13)*OUTP(13)
        NOUTP(14) = ODD(14,1)*OUTP(1) + ODD(14,14)*OUTP(14)
     
        OUTP = NOUTP

        NOUTP(0) = EVEN(0,0)*OUTP(0) + EVEN(0,13)*OUTP(13)
        NOUTP(1) = EVEN(1,1)*OUTP(1) + EVEN(1,12)*OUTP(12)
        NOUTP(2) = EVEN(2,2)*OUTP(2) + EVEN(2,15)*OUTP(15)
        NOUTP(3) = EVEN(3,3)*OUTP(3) + EVEN(3,14)*OUTP(14)

        NOUTP(4) = EVEN(4,4)*OUTP(4) + EVEN(4,9)*OUTP(9)
        NOUTP(5) = EVEN(5,5)*OUTP(5) + EVEN(5,8)*OUTP(8)
        NOUTP(6) = EVEN(6,6)*OUTP(6) + EVEN(6,11)*OUTP(11)
        NOUTP(7) = EVEN(7,7)*OUTP(7) + EVEN(7,10)*OUTP(10)

        NOUTP(8) = EVEN(8,5)*OUTP(5) + EVEN(8,8)*OUTP(8)
        NOUTP(9) = EVEN(9,4)*OUTP(4) + EVEN(9,9)*OUTP(9)
        NOUTP(10) = EVEN(10,7)*OUTP(7) + EVEN(10,10)*OUTP(10)
        NOUTP(11) = EVEN(11,6)*OUTP(6) + EVEN(11,11)*OUTP(11)

        NOUTP(12) = EVEN(12,1)*OUTP(1) + EVEN(12,12)*OUTP(12)
        NOUTP(13) = EVEN(13,0)*OUTP(0) + EVEN(13,13)*OUTP(13)
        NOUTP(14) = EVEN(14,3)*OUTP(3) + EVEN(14,14)*OUTP(14)
        NOUTP(15) = EVEN(15,2)*OUTP(2) + EVEN(15,15)*OUTP(15)

        OUTP = NOUTP

        do a = 0,15
           PSI(IND(a)) = OUTP(a)
        end do
      end do

      return
      end

      subroutine gen_other_bits_ssm(nbits,mp1,mp2,ip1,ip2,num,result)
      IMPLICIT NONE
      INTEGER nbits,mp1,mp2,ip1,ip2,num,result
      INTEGER bit_pos,source_bit,temp
      
      result = 0
      source_bit = 0
      temp = num
      
      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2
     .       .AND. bit_pos .NE. ip1 .AND. bit_pos .NE. ip2)
     .      then
            if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
          end if
      end do
      
      return
      end


C=====================================================================
C                       IDENTITY x MODE
C=====================================================================

      subroutine identity_mode(j,k,x,z,t3,t2,t1,PSI_IN,
     .                         LD,RABI,w,wqbt,w_k)
      INTEGER j,k,x,z
      DOUBLE PRECISION t3,t2,t1,RABI(7),LD(j,k),w(7),wqbt(7),w_k(7)
      DOUBLE PRECISION pmanny,pmannx,pmcrex,pmcrey,wj,wjp
      DOUBLE PRECISION pannmy,pannmx,pcremx,pcremy
      COMPLEX*16 i,cre,ann
      COMPLEX*16 PSI_IN(0:2**j*4**k-1)
      COMPLEX*16 H(0:7,0:7), H_C(0:7,0:7)
      
      i = (0.d0,1.d0)
      wj = w(x) - wqbt(x)
      wjp = w(x) - wqbt(x)
      

      pannmx = 5.d-1*RABI(x)*RABI(x)*LD(x,z)
     .       *(sin(wjp*t3 - (w_k(z)+wj)*t2)
     .       - sin(wjp*t3 - (w_k(z)+wj)*t1)
     .       + sin((w_k(z)+wj)*t3 - wjp*t2)
     .       - sin((w_k(z)+wj)*t3 - wjp*t1)
     .       + sin(wjp*t2 - (w_k(z)+wj)*t1)
     .       + sin((w_k(z)+wj)*t2 - wjp*t1))
     .       / ((w_k(z)+wj)*wjp)

      pannmy = 5.d-1*RABI(x)*RABI(x)*LD(x,z)
     .       *(-cos(wjp*t3 - (w_k(z)+wj)*t2)
     .       + cos(wjp*t3 - (w_k(z)+wj)*t1)
     .       + cos((w_k(z)+wj)*t3 - wjp*t2)
     .       - cos((w_k(z)+wj)*t3 - wjp*t1)
     .       - cos(wjp*t2 - (w_k(z)+wj)*t1)
     .       + cos((w_k(z)+wj)*t2 - wjp*t1))
     .       / ((w_k(z)+wj)*wjp)

      pcremx = 5.d-1*RABI(x)*RABI(x)*LD(x,z)
     .       *(-sin(wjp*t3 + (w_k(z)-wj)*t2)
     .       + sin(wjp*t3 + (w_k(z)-wj)*t1)
     .       - sin((wj-w_k(z))*t3 - wjp*t2)
     .       + sin((wj-w_k(z))*t3 - wjp*t1)
     .       - sin(wjp*t2 + (w_k(z)-wj)*t1)
     .       - sin((wj-w_k(z))*t2 - wjp*t1))
     .       / ((w_k(z)-wj)*wjp)

      pcremy = 5.d-1*RABI(x)*RABI(x)*LD(x,z)
     .       *(cos(wjp*t3 + (w_k(z)-wj)*t2)
     .       - cos(wjp*t3 + (w_k(z)-wj)*t1)
     .       - cos((wj-w_k(z))*t3 - wjp*t2)
     .       + cos((wj-w_k(z))*t3 - wjp*t1)
     .       + cos(wjp*t2 + (w_k(z)-wj)*t1)
     .       - cos((wj-w_k(z))*t2 - wjp*t1))
     .       / ((w_k(z)-wj)*wjp)

      pmanny = 5.d-1*(RABI(x)*RABI(x)*LD(x,z)/
     .        (wj*(w_k(z)-wjp)))*
     .  (cos((wjp-w_k(z))*t3 - wj*t2)
     .  - cos((wjp-w_k(z))*t3 - wj*t1)
     .  - cos(wj*t3 + (w_k(z)-wjp)*t2)
     .  + cos(wj*t3 + (w_k(z)-wjp)*t1)
     .  + cos((wjp-w_k(z))*t2 - wj*t1)
     .  - cos(wj*t2 + (w_k(z)-wjp)*t1)
     .  + cos((w_k(z)-wjp+wj)*t2)
     .  - cos((-w_k(z)+wjp-wj)*t2))

      pmannx = 5.d-1*(RABI(x)*RABI(x)*LD(x,z)/
     .        (wj*(w_k(z)-wjp)))*
     .  (sin((wjp-w_k(z))*t3 - wj*t2)
     .  - sin((wjp-w_k(z))*t3 - wj*t1)
     .  + sin(wj*t3 + (w_k(z)-wjp)*t2)
     .  - sin(wj*t3 + (w_k(z)-wjp)*t1)
     .  + sin((wjp-w_k(z))*t2 - wj*t1)
     .  + sin(wj*t2 + (w_k(z)-wjp)*t1)
     .  - sin((w_k(z)-wjp+wj)*t2)
     .  - sin((-w_k(z)+wjp-wj)*t2))

      pmcrey = 5.d-1*(RABI(x)*RABI(x)*LD(x,z)/
     .        (wj*(wjp+w_k(z))))*
     . (-cos((w_k(z)+wjp)*t3 - wj*t2)
     .  + cos((w_k(z)+wjp)*t3 - wj*t1)
     .  + cos(wj*t3 + (-w_k(z)-wjp)*t2)
     .  - cos(wj*t3 + (-w_k(z)-wjp)*t1)
     .  - cos((w_k(z)+wjp)*t2 - wj*t1)
     .  + cos(wj*t2 + (-w_k(z)-wjp)*t1)
     .  + cos((w_k(z)+wjp-wj)*t2)
     .  - cos((-w_k(z)-wjp+wj)*t2))

      pmcrex = 5.d-1*(RABI(x)*RABI(x)*LD(x,z)/
     .        (wj*(wjp+w_k(z))))*
     . (-sin((w_k(z)+wjp)*t3 - wj*t2)
     .  + sin((w_k(z)+wjp)*t3 - wj*t1)
     .  - sin(wj*t3 + (-w_k(z)-wjp)*t2)
     .  + sin(wj*t3 + (-w_k(z)-wjp)*t1)
     .  - sin((w_k(z)+wjp)*t2 - wj*t1)
     .  - sin(wj*t2 + (-w_k(z)-wjp)*t1)
     .  + sin((w_k(z)+wjp-wj)*t2)
     .  + sin((-w_k(z)-wjp+wj)*t2))

      cre = pmcrex + pcremx + i*(pmcrey + pcremy)
      ann = pmannx + pannmx + i*(pmanny + pannmy)

      H = 0.d0
      
      H(0,1) = cre
      H(1,2) = sqrt(2.d0)*cre
      H(2,3) = sqrt(3.d0)*cre

      H(1,0) = ann
      H(2,1) = sqrt(2.d0)*ann
      H(3,2) = sqrt(3.d0)*ann

      H(4,5) = conjg(cre)
      H(5,6) = sqrt(2.d0)*conjg(cre)
      H(6,7) = sqrt(3.d0)*conjg(cre)

      H(5,4) = conjg(ann)
      H(6,5) = sqrt(2.d0)*conjg(ann)
      H(7,6) = sqrt(3.d0)*conjg(ann)

      CALL EXACT_TAYLOR(7,H,H_C)

      CALL idmode_operators(H_C,PSI_IN,j,k,x,z)
      return 
      end

      subroutine idmode_operators(H,PSI,j,k,ion,mode)
      IMPLICIT NONE
      INTEGER a,b,c,j,k,ion,mode
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,ip
      INTEGER IND(0:7)
      COMPLEX*16 OUTP(0:7),PSI(0:2**j*4**k-1),H(0:7,0:7)

      mp1 = k*2 - mode*2
      mp2 = k*2 - mode*2 + 1
      ip  = k*2 + j - ion

      do ii = 0, 2**(j+2*k-3)-1
         call gen_other_bits_idm(j+2*k-1,mp1,mp2,ip,ii,base_bits)
         idx = 0

         do a = 0,1
            do b = 0,1
              do c = 0,1
                 bit_combo = base_bits
                 if (a .eq. 1) bit_combo = IBSET(bit_combo, ip)
                 if (b .eq. 1) bit_combo = IBSET(bit_combo, mp2)
                 if (c .eq. 1) bit_combo = IBSET(bit_combo, mp1)

                 IND(idx) = bit_combo
                 idx = idx + 1
              end do
            end do
        end do

      do a = 0,7
         OUTP(a) = PSI(IND(a))
      end do

      OUTP = MATMUL(H,OUTP)

      do a = 0,7
         PSI(IND(a)) = OUTP(a)
      end do
     
      end do 
      return
      end

      subroutine gen_other_bits_idm(nbits,mp1,mp2,ip,num,result)
      IMPLICIT NONE
      INTEGER nbits,mp1,mp2,ip,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2
     .      .AND. bit_pos .NE. ip) then
            if (BTEST(temp, source_bit)) then
                result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
         end if
      end do        

      return
      end

C=====================================================================
C                       SPIN x SPIN x MODE x MODE
C=====================================================================


      subroutine spin_spin_mode_mode(j,k,a,b,c,d,t3,t2,t1,PSI_IN,
     .                               LD,RABI,w,wqbt,w_k)
      INTEGER j,k,a,b,c,d
      DOUBLE PRECISION t3,t2,t1,RABI(7),LD(j,k),w(7),wqbt(7),w_k(7)
      DOUBLE PRECISION wj,wjp
      COMPLEX*16 ppcc,ppca,ppac,ppaa,pmcc,pmca,pmac,pmaa
      COMPLEX*16 i,PSI_IN(0:2**j*4**k-1)
      COMPLEX*16 H(0:63,0:63),H_C(0:63,0:63)

      i = (0.d0,1.d0)
      wj = w(a) - wqbt(a)
      wjp = w(b) - wqbt(b)
      
      
      ppac = 5.d-1*RABI(a)*RABI(b)*LD(b,c)*LD(b,d)*
     .        (-i*sin((-w_k(d)+w_k(c)+wjp)*t3 + wj*t2)
     .        + cos((-w_k(d)+w_k(c)+wjp)*t3 + wj*t2)
     .        + i*sin((-w_k(d)+w_k(c)+wjp)*t3 + wj*t1)
     .        - cos((-w_k(d)+w_k(c)+wjp)*t3 + wj*t1)
     .        + i*sin(wj*t3 + (-w_k(d)+w_k(c)+wjp)*t2)
     .        - cos(wj*t3 + (-w_k(d)+w_k(c)+wjp)*t2)
     .        - i*sin(wj*t3 + (-w_k(d)+w_k(c)+wjp)*t1)
     .        + cos(wj*t3 + (-w_k(d)+w_k(c)+wjp)*t1)
     .        - i*sin((-w_k(d)+w_k(c)+wjp)*t2 + wj*t1)
     .        + cos((-w_k(d)+w_k(c)+wjp)*t2 + wj*t1)
     .        + i*sin(wj*t2 + (-w_k(d)+w_k(c)+wjp)*t1)
     .        - cos(wj*t2 + (-w_k(d)+w_k(c)+wjp)*t1))
     .        / (wj*(w_k(d)-w_k(c)-wjp))
      
            
      pmca = 5.d-1*RABI(a)*RABI(b)*LD(b,c)*LD(b,d)*
     .        (i*sin((-w_k(d)+w_k(c)+wjp)*t3 - wj*t2)
     .        + cos((-w_k(d)+w_k(c)+wjp)*t3 - wj*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wjp)*t3 - wj*t1)
     .        - cos((-w_k(d)+w_k(c)+wjp)*t3 - wj*t1)
     .        + i*sin(wj*t3 + (w_k(d)-w_k(c)-wjp)*t2)
     .        - cos(wj*t3 + (w_k(d)-w_k(c)-wjp)*t2)
     .        - i*sin(wj*t3 + (w_k(d)-w_k(c)-wjp)*t1)
     .        + cos(wj*t3 + (w_k(d)-w_k(c)-wjp)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wjp)*t2 - wj*t1)
     .        + cos((-w_k(d)+w_k(c)+wjp)*t2 - wj*t1)
     .        + i*sin(wj*t2 + (w_k(d)-w_k(c)-wjp)*t1)
     .        - cos(wj*t2 + (w_k(d)-w_k(c)-wjp)*t1)
     .        - i*sin((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(d)+w_k(c)+wjp-wj)*t2))
     .        / (wj*(w_k(d)-w_k(c)-wjp))
      

      ppaa = 5.d-1*RABI(a)*RABI(b)*LD(b,c)*LD(b,d)*
     .        (i*sin((-w_k(d)+w_k(c)+wj)*t3 + wjp*t2)
     .        - cos((-w_k(d)+w_k(c)+wj)*t3 + wjp*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t3 + wjp*t1)
     .        + cos((-w_k(d)+w_k(c)+wj)*t3 + wjp*t1)
     .        - i*sin(wjp*t3 + (-w_k(d)+w_k(c)+wj)*t2)
     .        + cos(wjp*t3 + (-w_k(d)+w_k(c)+wj)*t2)
     .        + i*sin(wjp*t3 + (-w_k(d)+w_k(c)+wj)*t1)
     .        - cos(wjp*t3 + (-w_k(d)+w_k(c)+wj)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t2 + wjp*t1)
     .        - cos((-w_k(d)+w_k(c)+wj)*t2 + wjp*t1)
     .        - i*sin(wjp*t2 + (-w_k(d)+w_k(c)+wj)*t1)
     .        + cos(wjp*t2 + (-w_k(d)+w_k(c)+wj)*t1))
     .        / (wjp*(w_k(d)-w_k(c)-wj))



      pmaa = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .        (-i*sin((-w_k(d)+w_k(c)+wj)*t3 - wjp*t2)
     .        + cos((-w_k(d)+w_k(c)+wj)*t3 - wjp*t2)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t3 - wjp*t1)
     .        - cos((-w_k(d)+w_k(c)+wj)*t3 - wjp*t1)
     .        - i*sin(wjp*t3 + (w_k(d)-w_k(c)-wj)*t2)
     .        - cos(wjp*t3 + (w_k(d)-w_k(c)-wj)*t2)
     .        + i*sin(wjp*t3 + (w_k(d)-w_k(c)-wj)*t1)
     .        + cos(wjp*t3 + (w_k(d)-w_k(c)-wj)*t1)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t2 - wjp*t1)
     .        + cos((-w_k(d)+w_k(c)+wj)*t2 - wjp*t1)
     .        - i*sin(wjp*t2 + (w_k(d)-w_k(c)-wj)*t1)
     .        - cos(wjp*t2 + (w_k(d)-w_k(c)-wj)*t1)
     .        + i*sin((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(d)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(d)+w_k(c)-wjp+wj)*t2))
     .        / (wjp*(w_k(d)-w_k(c)-wj))

      ppaa = ppaa
     .     + 
     .        5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .        (i*sin((w_k(d)+wjp)*t3 + (w_k(c)+wj)*t2)
     .        - cos((w_k(d)+wjp)*t3 + (w_k(c)+wj)*t2)
     .        - i*sin((w_k(d)+wjp)*t3 + (w_k(c)+wj)*t1)
     .        + cos((w_k(d)+wjp)*t3 + (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t3 + (w_k(d)+wjp)*t2)
     .        + cos((w_k(c)+wj)*t3 + (w_k(d)+wjp)*t2)
     .        + i*sin((w_k(c)+wj)*t3 + (w_k(d)+wjp)*t1)
     .        - cos((w_k(c)+wj)*t3 + (w_k(d)+wjp)*t1)
     .        + i*sin((w_k(d)+wjp)*t2 + (w_k(c)+wj)*t1)
     .        - cos((w_k(d)+wjp)*t2 + (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t2 + (w_k(d)+wjp)*t1)
     .        + cos((w_k(c)+wj)*t2 + (w_k(d)+wjp)*t1))
     .        / ((w_k(c)+wj)*(w_k(d)+wjp))

      ppac = ppac
     .     + 
     .        5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .        (-i*sin((wjp-w_k(d))*t3 + (w_k(c)+wj)*t2)
     .        + cos((wjp-w_k(d))*t3 + (w_k(c)+wj)*t2)
     .        + i*sin((wjp-w_k(d))*t3 + (w_k(c)+wj)*t1)
     .        - cos((wjp-w_k(d))*t3 + (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t3 + (wjp-w_k(d))*t2)
     .        - cos((w_k(c)+wj)*t3 + (wjp-w_k(d))*t2)
     .        - i*sin((w_k(c)+wj)*t3 + (wjp-w_k(d))*t1)
     .        + cos((w_k(c)+wj)*t3 + (wjp-w_k(d))*t1)
     .        - i*sin((wjp-w_k(d))*t2 + (w_k(c)+wj)*t1)
     .        + cos((wjp-w_k(d))*t2 + (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t2 + (wjp-w_k(d))*t1)
     .        - cos((w_k(c)+wj)*t2 + (wjp-w_k(d))*t1))
     .        / ((w_k(c)+wj)*(w_k(d)-wjp))


      pmac = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .        (-i*sin((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        - cos((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        + i*sin((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        + cos((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t2)
     .        + cos((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t2)
     .        + i*sin((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t1)
     .        - cos((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t1)
     .        - i*sin((w_k(d)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - cos((w_k(d)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t2 - (w_k(d)+wjp)*t1)
     .        + cos((w_k(c)+wj)*t2 - (w_k(d)+wjp)*t1)
     .        + i*sin((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(d)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(d)+w_k(c)-wjp+wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(d)+wjp))

      pmaa = pmaa
     .     +  5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .        (i*sin((wjp-w_k(d))*t3 - (w_k(c)+wj)*t2)
     .        + cos((wjp-w_k(d))*t3 - (w_k(c)+wj)*t2)
     .        - i*sin((wjp-w_k(d))*t3 - (w_k(c)+wj)*t1)
     .        - cos((wjp-w_k(d))*t3 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t2)
     .        - cos((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t2)
     .        - i*sin((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t1)
     .        + cos((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t1)
     .        + i*sin((wjp-w_k(d))*t2 - (w_k(c)+wj)*t1)
     .        + cos((wjp-w_k(d))*t2 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t2 + (w_k(d)-wjp)*t1)
     .        - cos((w_k(c)+wj)*t2 + (w_k(d)-wjp)*t1)
     .        - i*sin((w_k(d)+w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(d)+w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(d)-w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(d)-w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(d)-wjp))


      ppca = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .        (i*sin((w_k(d)+wjp)*t3 + (wj-w_k(c))*t2)
     .        - cos((w_k(d)+wjp)*t3 + (wj-w_k(c))*t2)
     .        - i*sin((w_k(d)+wjp)*t3 + (wj-w_k(c))*t1)
     .        + cos((w_k(d)+wjp)*t3 + (wj-w_k(c))*t1)
     .        - i*sin((wj-w_k(c))*t3 + (w_k(d)+wjp)*t2)
     .        + cos((wj-w_k(c))*t3 + (w_k(d)+wjp)*t2)
     .        + i*sin((wj-w_k(c))*t3 + (w_k(d)+wjp)*t1)
     .        - cos((wj-w_k(c))*t3 + (w_k(d)+wjp)*t1)
     .        + i*sin((w_k(d)+wjp)*t2 + (wj-w_k(c))*t1)
     .        - cos((w_k(d)+wjp)*t2 + (wj-w_k(c))*t1)
     .        - i*sin((wj-w_k(c))*t2 + (w_k(d)+wjp)*t1)
     .        + cos((wj-w_k(c))*t2 + (w_k(d)+wjp)*t1))
     .        / ((w_k(c)-wj)*(w_k(d)+wjp))

      pmcc = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .       (-i*sin((wjp-w_k(d))*t3 + (wj-w_k(c))*t2)
     .        + cos((wjp-w_k(d))*t3 + (wj-w_k(c))*t2)
     .        + i*sin((wjp-w_k(d))*t3 + (wj-w_k(c))*t1)
     .        - cos((wjp-w_k(d))*t3 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t3 + (wjp-w_k(d))*t2)
     .        - cos((wj-w_k(c))*t3 + (wjp-w_k(d))*t2)
     .        - i*sin((wj-w_k(c))*t3 + (wjp-w_k(d))*t1)
     .        + cos((wj-w_k(c))*t3 + (wjp-w_k(d))*t1)
     .        - i*sin((wjp-w_k(d))*t2 + (wj-w_k(c))*t1)
     .        + cos((wjp-w_k(d))*t2 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t2 + (wjp-w_k(d))*t1)
     .        - cos((wj-w_k(c))*t2 + (wjp-w_k(d))*t1))
     .        / ((w_k(c)-wj)*(w_k(d)-wjp))


      pmca = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,d)*
     .       -(i*sin((wjp-w_k(d))*t3 + (w_k(c)-wj)*t2)
     .        + cos((wjp-w_k(d))*t3 + (w_k(c)-wj)*t2)
     .        - i*sin((wjp-w_k(d))*t3 + (w_k(c)-wj)*t1)
     .        - cos((wjp-w_k(d))*t3 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t3 + (w_k(d)-wjp)*t2)
     .        - cos((wj-w_k(c))*t3 + (w_k(d)-wjp)*t2)
     .        - i*sin((wj-w_k(c))*t3 + (w_k(d)-wjp)*t1)
     .        + cos((wj-w_k(c))*t3 + (w_k(d)-wjp)*t1)
     .        + i*sin((wjp-w_k(d))*t2 + (w_k(c)-wj)*t1)
     .        + cos((wjp-w_k(d))*t2 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t2 + (w_k(d)-wjp)*t1)
     .        - cos((wj-w_k(c))*t2 + (w_k(d)-wjp)*t1)
     .        - i*sin((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(d)+w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)-wj)*(w_k(d)-wjp))



      H = (0.d0,0.d0)

C=======================================================================
C     pp block: rows 0:15, cols 48:63
C=======================================================================

      H(0,53) = ppcc

      H(1,52) = ppca
      H(1,54) = sqrt(2.d0)*ppcc

      H(2,53) = sqrt(2.d0)*ppca
      H(2,55) = sqrt(3.d0)*ppcc

      H(3,54) = sqrt(3.d0)*ppca

      H(4,49) = ppac
      H(4,57) = sqrt(2.d0)*ppcc

      H(5,48) = ppaa
      H(5,50) = sqrt(2.d0)*ppac
      H(5,56) = sqrt(2.d0)*ppca
      H(5,58) = sqrt(4.d0)*ppcc

      H(6,49) = sqrt(2.d0)*ppaa
      H(6,51) = sqrt(3.d0)*ppac
      H(6,57) = sqrt(4.d0)*ppca
      H(6,59) = sqrt(6.d0)*ppcc

      H(7,50) = sqrt(3.d0)*ppaa
      H(7,58) = sqrt(6.d0)*ppca

      H(8,53) = sqrt(2.d0)*ppac
      H(8,61) = sqrt(3.d0)*ppcc

      H(9,52) = sqrt(2.d0)*ppaa
      H(9,54) = sqrt(4.d0)*ppac
      H(9,60) = sqrt(3.d0)*ppca
      H(9,62) = sqrt(6.d0)*ppcc

      H(10,53) = sqrt(4.d0)*ppaa
      H(10,55) = sqrt(6.d0)*ppac
      H(10,61) = sqrt(6.d0)*ppca
      H(10,63) = sqrt(9.d0)*ppcc

      H(11,54) = sqrt(6.d0)*ppaa
      H(11,62) = sqrt(9.d0)*ppca

      H(12,57) = sqrt(3.d0)*ppac

      H(13,56) = sqrt(3.d0)*ppaa
      H(13,58) = sqrt(6.d0)*ppac

      H(14,57) = sqrt(6.d0)*ppaa
      H(14,59) = sqrt(9.d0)*ppac

      H(15,58) = sqrt(9.d0)*ppaa


C=======================================================================
C     pm block: rows 16:31, cols 32:47
C=======================================================================

      H(16,37) = pmcc

      H(17,36) = pmca
      H(17,38) = sqrt(2.d0)*pmcc

      H(18,37) = sqrt(2.d0)*pmca
      H(18,39) = sqrt(3.d0)*pmcc

      H(19,38) = sqrt(3.d0)*pmca

      H(20,33) = pmac
      H(20,41) = sqrt(2.d0)*pmcc

      H(21,32) = pmaa
      H(21,34) = sqrt(2.d0)*pmac
      H(21,40) = sqrt(2.d0)*pmca
      H(21,42) = sqrt(4.d0)*pmcc

      H(22,33) = sqrt(2.d0)*pmaa
      H(22,35) = sqrt(3.d0)*pmac
      H(22,41) = sqrt(4.d0)*pmca
      H(22,43) = sqrt(6.d0)*pmcc

      H(23,34) = sqrt(3.d0)*pmaa
      H(23,42) = sqrt(6.d0)*pmca

      H(24,37) = sqrt(2.d0)*pmac
      H(24,45) = sqrt(3.d0)*pmcc

      H(25,36) = sqrt(2.d0)*pmaa
      H(25,38) = sqrt(4.d0)*pmac
      H(25,44) = sqrt(3.d0)*pmca
      H(25,46) = sqrt(6.d0)*pmcc

      H(26,37) = sqrt(4.d0)*pmaa
      H(26,39) = sqrt(6.d0)*pmac
      H(26,45) = sqrt(6.d0)*pmca
      H(26,47) = sqrt(9.d0)*pmcc

      H(27,38) = sqrt(6.d0)*pmaa
      H(27,46) = sqrt(9.d0)*pmca

      H(28,41) = sqrt(3.d0)*pmac

      H(29,40) = sqrt(3.d0)*pmaa
      H(29,42) = sqrt(6.d0)*pmac

      H(30,41) = sqrt(6.d0)*pmaa
      H(30,43) = sqrt(9.d0)*pmac

      H(31,42) = sqrt(9.d0)*pmaa


C=======================================================================
C     mp block: rows 32:47, cols 16:31
C     Hermitian conjugate of pm block
C=======================================================================

      H(37,16) = conjg(pmcc)

      H(36,17) = conjg(pmca)
      H(38,17) = sqrt(2.d0)*conjg(pmcc)

      H(37,18) = sqrt(2.d0)*conjg(pmca)
      H(39,18) = sqrt(3.d0)*conjg(pmcc)

      H(38,19) = sqrt(3.d0)*conjg(pmca)

      H(33,20) = conjg(pmac)
      H(41,20) = sqrt(2.d0)*conjg(pmcc)

      H(32,21) = conjg(pmaa)
      H(34,21) = sqrt(2.d0)*conjg(pmac)
      H(40,21) = sqrt(2.d0)*conjg(pmca)
      H(42,21) = sqrt(4.d0)*conjg(pmcc)

      H(33,22) = sqrt(2.d0)*conjg(pmaa)
      H(35,22) = sqrt(3.d0)*conjg(pmac)
      H(41,22) = sqrt(4.d0)*conjg(pmca)
      H(43,22) = sqrt(6.d0)*conjg(pmcc)

      H(34,23) = sqrt(3.d0)*conjg(pmaa)
      H(42,23) = sqrt(6.d0)*conjg(pmca)

      H(37,24) = sqrt(2.d0)*conjg(pmac)
      H(45,24) = sqrt(3.d0)*conjg(pmcc)

      H(36,25) = sqrt(2.d0)*conjg(pmaa)
      H(38,25) = sqrt(4.d0)*conjg(pmac)
      H(44,25) = sqrt(3.d0)*conjg(pmca)
      H(46,25) = sqrt(6.d0)*conjg(pmcc)

      H(37,26) = sqrt(4.d0)*conjg(pmaa)
      H(39,26) = sqrt(6.d0)*conjg(pmac)
      H(45,26) = sqrt(6.d0)*conjg(pmca)
      H(47,26) = sqrt(9.d0)*conjg(pmcc)

      H(38,27) = sqrt(6.d0)*conjg(pmaa)
      H(46,27) = sqrt(9.d0)*conjg(pmca)

      H(41,28) = sqrt(3.d0)*conjg(pmac)

      H(40,29) = sqrt(3.d0)*conjg(pmaa)
      H(42,29) = sqrt(6.d0)*conjg(pmac)

      H(41,30) = sqrt(6.d0)*conjg(pmaa)
      H(43,30) = sqrt(9.d0)*conjg(pmac)

      H(42,31) = sqrt(9.d0)*conjg(pmaa)


C=======================================================================
C     mm block: rows 48:63, cols 0:15
C     Hermitian conjugate of pp block
C=======================================================================

      H(53,0) = conjg(ppcc)

      H(52,1) = conjg(ppca)
      H(54,1) = sqrt(2.d0)*conjg(ppcc)

      H(53,2) = sqrt(2.d0)*conjg(ppca)
      H(55,2) = sqrt(3.d0)*conjg(ppcc)

      H(54,3) = sqrt(3.d0)*conjg(ppca)

      H(49,4) = conjg(ppac)
      H(57,4) = sqrt(2.d0)*conjg(ppcc)

      H(48,5) = conjg(ppaa)
      H(50,5) = sqrt(2.d0)*conjg(ppac)
      H(56,5) = sqrt(2.d0)*conjg(ppca)
      H(58,5) = sqrt(4.d0)*conjg(ppcc)

      H(49,6) = sqrt(2.d0)*conjg(ppaa)
      H(51,6) = sqrt(3.d0)*conjg(ppac)
      H(57,6) = sqrt(4.d0)*conjg(ppca)
      H(59,6) = sqrt(6.d0)*conjg(ppcc)

      H(50,7) = sqrt(3.d0)*conjg(ppaa)
      H(58,7) = sqrt(6.d0)*conjg(ppca)

      H(53,8) = sqrt(2.d0)*conjg(ppac)
      H(61,8) = sqrt(3.d0)*conjg(ppcc)

      H(52,9) = sqrt(2.d0)*conjg(ppaa)
      H(54,9) = sqrt(4.d0)*conjg(ppac)
      H(60,9) = sqrt(3.d0)*conjg(ppca)
      H(62,9) = sqrt(6.d0)*conjg(ppcc)

      H(53,10) = sqrt(4.d0)*conjg(ppaa)
      H(55,10) = sqrt(6.d0)*conjg(ppac)
      H(61,10) = sqrt(6.d0)*conjg(ppca)
      H(63,10) = sqrt(9.d0)*conjg(ppcc)

      H(54,11) = sqrt(6.d0)*conjg(ppaa)
      H(62,11) = sqrt(9.d0)*conjg(ppca)

      H(57,12) = sqrt(3.d0)*conjg(ppac)

      H(56,13) = sqrt(3.d0)*conjg(ppaa)
      H(58,13) = sqrt(6.d0)*conjg(ppac)

      H(57,14) = sqrt(6.d0)*conjg(ppaa)
      H(59,14) = sqrt(9.d0)*conjg(ppac)

      H(58,15) = sqrt(9.d0)*conjg(ppaa)


      CALL EXACT_TAYLOR(63,H,H_C)
      
      CALL ssmm_operators(H_C,PSI_IN,j,k,a,b,c,d) 
      
      return
      end

      subroutine ssmm_operators(H,PSI,j,k,ion1,ion2,mode1,mode2)
      IMPLICIT NONE
      INTEGER a,b,c,d,f,e,j,k,ion1,ion2,mode1,mode2
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,mp3,mp4,ip1,ip2
      INTEGER IND(0:63)
      COMPLEX*16 OUTP(0:63)
      COMPLEX*16 PSI(0:2**j*4**k-1),H(0:63,0:63)

      mp1 = k*2-mode1*2
      mp2 = k*2-mode1*2+1
      mp3 = k*2-mode2*2
      mp4 = k*2-mode2*2+1
      ip1 = k*2+j-ion1
      ip2 = k*2+j-ion2

      do ii = 0, 2**(j+2*k-6)-1
         call gen_other_bits_ssmm(j+2*k-1,mp1,mp2,mp3,mp4,ip1,ip2,
     .                            ii,base_bits)
         idx = 0

         do a = 0,1
            do b = 0,1
              do c = 0,1
                do d = 0,1
                  do e = 0,1
                    do f = 0,1
                    bit_combo = base_bits
                    if (a .eq. 1) bit_combo = IBSET(bit_combo, ip1)
                    if (b .eq. 1) bit_combo = IBSET(bit_combo, ip2)
                    if (c .eq. 1) bit_combo = IBSET(bit_combo, mp4)
                    if (d .eq. 1) bit_combo = IBSET(bit_combo, mp3)
                    if (e .eq. 1) bit_combo = IBSET(bit_combo, mp2)
                    if (f .eq. 1) bit_combo = IBSET(bit_combo, mp1)

            
                    IND(idx) = bit_combo
                    idx = idx + 1
                    end do
                  end do
                 end do
               end do
            end do
         end do

         do a = 0,63
            OUTP(a) = PSI(IND(a))
         end do

         OUTP = MATMUL(H,OUTP)

         do a = 0,63
            PSI(IND(a)) = OUTP(a)
         end do
      end do
                        

      return
      end

      subroutine gen_other_bits_ssmm(nbits,mp1,mp2,mp3,mp4,ip1,ip2,
     .                               num,result)
      IMPLICIT NONE
      INTEGER nbits,mp4,mp3,mp2,mp1,ip1,ip2,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2
     .       .AND. bit_pos .NE. mp3 .AND. bit_pos .NE. mp4
     .       .AND. bit_pos .NE. ip1 .AND. bit_pos .NE. ip2)
     .      then
            if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
          end if
      end do

      return
      end


C=====================================================================
C                       IDENTITY X MODE X MODE
C=====================================================================

      subroutine identity_mode_mode(j,k,a,c,d,t3,t2,t1,PSI_IN,
     .                              LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,a,c,d
      DOUBLE PRECISION t1,t2,t3
      DOUBLE PRECISION LD(j,k),RABI(7),w(7),wqbt(7),w_k(7),wj,wjp
      COMPLEX*16 cc,ca,ac,aa,pmcc,pmca,pmac,pmaa,i
      COMPLEX*16 H(0:31,0:31),PSI_IN(0:2**j*4**k-1),H_C(0:31,0:31)

      i = (0.d0,1.d0)
      wj = w(a) - wqbt(a)
      wjp = wj

      pmca = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,d)*
     .        (i*sin((-w_k(d)+w_k(c)+wjp)*t3 - wj*t2)
     .        + cos((-w_k(d)+w_k(c)+wjp)*t3 - wj*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wjp)*t3 - wj*t1)
     .        - cos((-w_k(d)+w_k(c)+wjp)*t3 - wj*t1)
     .        + i*sin(wj*t3 + (w_k(d)-w_k(c)-wjp)*t2)
     .        - cos(wj*t3 + (w_k(d)-w_k(c)-wjp)*t2)
     .        - i*sin(wj*t3 + (w_k(d)-w_k(c)-wjp)*t1)
     .        + cos(wj*t3 + (w_k(d)-w_k(c)-wjp)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wjp)*t2 - wj*t1)
     .        + cos((-w_k(d)+w_k(c)+wjp)*t2 - wj*t1)
     .        + i*sin(wj*t2 + (w_k(d)-w_k(c)-wjp)*t1)
     .        - cos(wj*t2 + (w_k(d)-w_k(c)-wjp)*t1)
     .        - i*sin((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(d)+w_k(c)+wjp-wj)*t2))
     .        / (wj*(w_k(d)-w_k(c)-wjp))
      


      pmaa = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,d)*
     .        (-i*sin((-w_k(d)+w_k(c)+wj)*t3 - wjp*t2)
     .        + cos((-w_k(d)+w_k(c)+wj)*t3 - wjp*t2)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t3 - wjp*t1)
     .        - cos((-w_k(d)+w_k(c)+wj)*t3 - wjp*t1)
     .        - i*sin(wjp*t3 + (w_k(d)-w_k(c)-wj)*t2)
     .        - cos(wjp*t3 + (w_k(d)-w_k(c)-wj)*t2)
     .        + i*sin(wjp*t3 + (w_k(d)-w_k(c)-wj)*t1)
     .        + cos(wjp*t3 + (w_k(d)-w_k(c)-wj)*t1)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t2 - wjp*t1)
     .        + cos((-w_k(d)+w_k(c)+wj)*t2 - wjp*t1)
     .        - i*sin(wjp*t2 + (w_k(d)-w_k(c)-wj)*t1)
     .        - cos(wjp*t2 + (w_k(d)-w_k(c)-wj)*t1)
     .        + i*sin((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(d)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(d)+w_k(c)-wjp+wj)*t2))
     .        / (wjp*(w_k(d)-w_k(c)-wj))

      pmac = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,d)*
     .        (-i*sin((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        - cos((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        + i*sin((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        + cos((w_k(d)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t2)
     .        + cos((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t2)
     .        + i*sin((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t1)
     .        - cos((w_k(c)+wj)*t3 - (w_k(d)+wjp)*t1)
     .        - i*sin((w_k(d)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - cos((w_k(d)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t2 - (w_k(d)+wjp)*t1)
     .        + cos((w_k(c)+wj)*t2 - (w_k(d)+wjp)*t1)
     .        + i*sin((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(d)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(d)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(d)+w_k(c)-wjp+wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(d)+wjp))

      pmaa = pmaa
     .     +  5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,d)*
     .        (i*sin((wjp-w_k(d))*t3 - (w_k(c)+wj)*t2)
     .        + cos((wjp-w_k(d))*t3 - (w_k(c)+wj)*t2)
     .        - i*sin((wjp-w_k(d))*t3 - (w_k(c)+wj)*t1)
     .        - cos((wjp-w_k(d))*t3 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t2)
     .        - cos((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t2)
     .        - i*sin((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t1)
     .        + cos((w_k(c)+wj)*t3 + (w_k(d)-wjp)*t1)
     .        + i*sin((wjp-w_k(d))*t2 - (w_k(c)+wj)*t1)
     .        + cos((wjp-w_k(d))*t2 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t2 + (w_k(d)-wjp)*t1)
     .        - cos((w_k(c)+wj)*t2 + (w_k(d)-wjp)*t1)
     .        - i*sin((w_k(d)+w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(d)+w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(d)-w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(d)-w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(d)-wjp))

      pmcc = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,d)*
     .       (-i*sin((wjp-w_k(d))*t3 + (wj-w_k(c))*t2)
     .        + cos((wjp-w_k(d))*t3 + (wj-w_k(c))*t2)
     .        + i*sin((wjp-w_k(d))*t3 + (wj-w_k(c))*t1)
     .        - cos((wjp-w_k(d))*t3 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t3 + (wjp-w_k(d))*t2)
     .        - cos((wj-w_k(c))*t3 + (wjp-w_k(d))*t2)
     .        - i*sin((wj-w_k(c))*t3 + (wjp-w_k(d))*t1)
     .        + cos((wj-w_k(c))*t3 + (wjp-w_k(d))*t1)
     .        - i*sin((wjp-w_k(d))*t2 + (wj-w_k(c))*t1)
     .        + cos((wjp-w_k(d))*t2 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t2 + (wjp-w_k(d))*t1)
     .        - cos((wj-w_k(c))*t2 + (wjp-w_k(d))*t1))
     .        / ((w_k(c)-wj)*(w_k(d)-wjp))


      pmca = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,d)*
     .       -(i*sin((wjp-w_k(d))*t3 + (w_k(c)-wj)*t2)
     .        + cos((wjp-w_k(d))*t3 + (w_k(c)-wj)*t2)
     .        - i*sin((wjp-w_k(d))*t3 + (w_k(c)-wj)*t1)
     .        - cos((wjp-w_k(d))*t3 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t3 + (w_k(d)-wjp)*t2)
     .        - cos((wj-w_k(c))*t3 + (w_k(d)-wjp)*t2)
     .        - i*sin((wj-w_k(c))*t3 + (w_k(d)-wjp)*t1)
     .        + cos((wj-w_k(c))*t3 + (w_k(d)-wjp)*t1)
     .        + i*sin((wjp-w_k(d))*t2 + (w_k(c)-wj)*t1)
     .        + cos((wjp-w_k(d))*t2 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t2 + (w_k(d)-wjp)*t1)
     .        - cos((wj-w_k(c))*t2 + (w_k(d)-wjp)*t1)
     .        - i*sin((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(d)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(d)+w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)-wj)*(w_k(d)-wjp))



      cc = pmcc
      ca = pmca
      ac = pmac
      aa = pmaa


      H = (0.d0,0.d0)

C=======================================================================
C     top-left block: rows 0:15, cols 0:15
C=======================================================================

      H(0,5) = cc

      H(1,4) = ca
      H(1,6) = sqrt(2.d0)*cc

      H(2,5) = sqrt(2.d0)*ca
      H(2,7) = sqrt(3.d0)*cc

      H(3,6) = sqrt(3.d0)*ca

      H(4,1) = ac
      H(4,9) = sqrt(2.d0)*cc

      H(5,0) = aa
      H(5,2) = sqrt(2.d0)*ac
      H(5,8) = sqrt(2.d0)*ca
      H(5,10) = sqrt(4.d0)*cc

      H(6,1) = sqrt(2.d0)*aa
      H(6,3) = sqrt(3.d0)*ac
      H(6,9) = sqrt(4.d0)*ca
      H(6,11) = sqrt(6.d0)*cc

      H(7,2) = sqrt(3.d0)*aa
      H(7,10) = sqrt(6.d0)*ca

      H(8,5) = sqrt(2.d0)*ac
      H(8,13) = sqrt(3.d0)*cc

      H(9,4) = sqrt(2.d0)*aa
      H(9,6) = sqrt(4.d0)*ac
      H(9,12) = sqrt(3.d0)*ca
      H(9,14) = sqrt(6.d0)*cc

      H(10,5) = sqrt(4.d0)*aa
      H(10,7) = sqrt(6.d0)*ac
      H(10,13) = sqrt(6.d0)*ca
      H(10,15) = sqrt(9.d0)*cc

      H(11,6) = sqrt(6.d0)*aa
      H(11,14) = sqrt(9.d0)*ca

      H(12,9) = sqrt(3.d0)*ac

      H(13,8) = sqrt(3.d0)*aa
      H(13,10) = sqrt(6.d0)*ac

      H(14,9) = sqrt(6.d0)*aa
      H(14,11) = sqrt(9.d0)*ac

      H(15,10) = sqrt(9.d0)*aa


C=======================================================================
C     bottom-right block: rows 16:31, cols 16:31
C=======================================================================

      H(16,21) = cc

      H(17,20) = ca
      H(17,22) = sqrt(2.d0)*cc

      H(18,21) = sqrt(2.d0)*ca
      H(18,23) = sqrt(3.d0)*cc

      H(19,22) = sqrt(3.d0)*ca

      H(20,17) = ac
      H(20,25) = sqrt(2.d0)*cc

      H(21,16) = aa
      H(21,18) = sqrt(2.d0)*ac
      H(21,24) = sqrt(2.d0)*ca
      H(21,26) = sqrt(4.d0)*cc

      H(22,17) = sqrt(2.d0)*aa
      H(22,19) = sqrt(3.d0)*ac
      H(22,25) = sqrt(4.d0)*ca
      H(22,27) = sqrt(6.d0)*cc

      H(23,18) = sqrt(3.d0)*aa
      H(23,26) = sqrt(6.d0)*ca

      H(24,21) = sqrt(2.d0)*ac
      H(24,29) = sqrt(3.d0)*cc

      H(25,20) = sqrt(2.d0)*aa
      H(25,22) = sqrt(4.d0)*ac
      H(25,28) = sqrt(3.d0)*ca
      H(25,30) = sqrt(6.d0)*cc

      H(26,21) = sqrt(4.d0)*aa
      H(26,23) = sqrt(6.d0)*ac
      H(26,29) = sqrt(6.d0)*ca
      H(26,31) = sqrt(9.d0)*cc

      H(27,22) = sqrt(6.d0)*aa
      H(27,30) = sqrt(9.d0)*ca

      H(28,25) = sqrt(3.d0)*ac

      H(29,24) = sqrt(3.d0)*aa
      H(29,26) = sqrt(6.d0)*ac

      H(30,25) = sqrt(6.d0)*aa
      H(30,27) = sqrt(9.d0)*ac

      H(31,26) = sqrt(9.d0)*aa

      CALL EXACT_TAYLOR(31,H,H_C)
      CALL idmm_operators(H_C,PSI_IN,j,k,a,c,d)

      return
      end

      subroutine idmm_operators(H,PSI,j,k,ion,mode1,mode2)
      IMPLICIT NONE
      INTEGER j,k,a,b,c,d,e,ion,mode1,mode2
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,mp3,mp4,ip,IND(0:31)
      COMPLEX*16 OUTP(0:31),H(0:31,0:31)
      COMPLEX*16 PSI(0:2**j*4**k-1)

      mp1 = k*2-mode2*2
      mp2 = k*2-mode2*2+1
      mp3 = k*2-mode1*2
      mp4 = k*2-mode1*2+1
      ip = k*2+j-ion

      do ii = 0, 2**(j+2*k-5)-1      
         call gen_other_bits_idmm(j+2*k-1,mp1,mp2,mp3,mp4,ip,
     .                             ii,base_bits)
         idx = 0
         do a = 0,1
           do b = 0,1
             do c = 0,1
               do d = 0,1
                 do e = 0,1
                  bit_combo = base_bits
                  if (a .EQ. 1) bit_combo = IBSET(bit_combo, ip)
                  if (b .EQ. 1) bit_combo = IBSET(bit_combo, mp4)
                  if (c .EQ. 1) bit_combo = IBSET(bit_combo, mp3)
                  if (d .EQ. 1) bit_combo = IBSET(bit_combo, mp2)
                  if (e .EQ. 1) bit_combo = IBSET(bit_combo, mp1)
                  
                  IND(idx) = bit_combo
                  idx = idx + 1
                  end do
               end do
             end do
           end do
         end do

         do a = 0,31
            OUTP(a) = PSI(IND(a))
         end do
      
         OUTP = MATMUL(H,OUTP)

        do a = 0,31
           PSI(IND(a)) = OUTP(a)
        end do
      end do   
                  
      
      return
      end

      subroutine gen_other_bits_idmm(nbits,mp1,mp2,mp3,mp4,ip,
     .                               num,result)
      IMPLICIT NONE
      INTEGER nbits,mp4,mp3,mp2,mp1,ip,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2
     .       .AND. bit_pos .NE. mp3 .AND. bit_pos .NE. mp4
     ,       .AND. bit_pos .NE. ip) then
           if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
           end if
           source_bit = source_bit + 1
         end if
      end do
      return
      end


C=====================================================================
C                       SPIN x SPIN x SINGLE OPERATOR
C=====================================================================
      

      subroutine spin_spin_single(j,k,a,b,c,t3,t2,t1,PSI_IN,
     .                              LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,a,b,c
      DOUBLE PRECISION t1,t2,t3
      DOUBLE PRECISION LD(j,k),RABI(7),w(7),wqbt(7),w_k(7),wj,wjp
      COMPLEX*16 ppaa,ppac,ppca,ppcc,pmaa,pmac,pmca,pmcc,i
      COMPLEX*16 H(0:15,0:15), H_C(0:15,0:15), PSI_IN(0:2**j*4**k-1)


       i = (0.d0,1.d0)
      wj = w(a) - wqbt(a)
      wjp = w(b) - wqbt(b)
      
      
      ppac = 5.d-1*RABI(a)*RABI(b)*LD(b,c)*LD(b,c)*
     .        (-i*sin((-w_k(c)+w_k(c)+wjp)*t3 + wj*t2)
     .        + cos((-w_k(c)+w_k(c)+wjp)*t3 + wj*t2)
     .        + i*sin((-w_k(c)+w_k(c)+wjp)*t3 + wj*t1)
     .        - cos((-w_k(c)+w_k(c)+wjp)*t3 + wj*t1)
     .        + i*sin(wj*t3 + (-w_k(c)+w_k(c)+wjp)*t2)
     .        - cos(wj*t3 + (-w_k(c)+w_k(c)+wjp)*t2)
     .        - i*sin(wj*t3 + (-w_k(c)+w_k(c)+wjp)*t1)
     .        + cos(wj*t3 + (-w_k(c)+w_k(c)+wjp)*t1)
     .        - i*sin((-w_k(c)+w_k(c)+wjp)*t2 + wj*t1)
     .        + cos((-w_k(c)+w_k(c)+wjp)*t2 + wj*t1)
     .        + i*sin(wj*t2 + (-w_k(c)+w_k(c)+wjp)*t1)
     .        - cos(wj*t2 + (-w_k(c)+w_k(c)+wjp)*t1))
     .        / (wj*(w_k(c)-w_k(c)-wjp))
      
            
      pmca = 5.d-1*RABI(a)*RABI(b)*LD(b,c)*LD(b,c)*
     .        (i*sin((-w_k(c)+w_k(c)+wjp)*t3 - wj*t2)
     .        + cos((-w_k(c)+w_k(c)+wjp)*t3 - wj*t2)
     .        - i*sin((-w_k(c)+w_k(c)+wjp)*t3 - wj*t1)
     .        - cos((-w_k(c)+w_k(c)+wjp)*t3 - wj*t1)
     .        + i*sin(wj*t3 + (w_k(c)-w_k(c)-wjp)*t2)
     .        - cos(wj*t3 + (w_k(c)-w_k(c)-wjp)*t2)
     .        - i*sin(wj*t3 + (w_k(c)-w_k(c)-wjp)*t1)
     .        + cos(wj*t3 + (w_k(c)-w_k(c)-wjp)*t1)
     .        + i*sin((-w_k(c)+w_k(c)+wjp)*t2 - wj*t1)
     .        + cos((-w_k(c)+w_k(c)+wjp)*t2 - wj*t1)
     .        + i*sin(wj*t2 + (w_k(c)-w_k(c)-wjp)*t1)
     .        - cos(wj*t2 + (w_k(c)-w_k(c)-wjp)*t1)
     .        - i*sin((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(c)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(c)+w_k(c)+wjp-wj)*t2))
     .        / (wj*(w_k(c)-w_k(c)-wjp))
      

      ppaa = 5.d-1*RABI(a)*RABI(b)*LD(b,c)*LD(b,c)*
     .        (i*sin((-w_k(c)+w_k(c)+wj)*t3 + wjp*t2)
     .        - cos((-w_k(c)+w_k(c)+wj)*t3 + wjp*t2)
     .        - i*sin((-w_k(c)+w_k(c)+wj)*t3 + wjp*t1)
     .        + cos((-w_k(c)+w_k(c)+wj)*t3 + wjp*t1)
     .        - i*sin(wjp*t3 + (-w_k(c)+w_k(c)+wj)*t2)
     .        + cos(wjp*t3 + (-w_k(c)+w_k(c)+wj)*t2)
     .        + i*sin(wjp*t3 + (-w_k(c)+w_k(c)+wj)*t1)
     .        - cos(wjp*t3 + (-w_k(c)+w_k(c)+wj)*t1)
     .        + i*sin((-w_k(c)+w_k(c)+wj)*t2 + wjp*t1)
     .        - cos((-w_k(c)+w_k(c)+wj)*t2 + wjp*t1)
     .        - i*sin(wjp*t2 + (-w_k(c)+w_k(c)+wj)*t1)
     .        + cos(wjp*t2 + (-w_k(c)+w_k(c)+wj)*t1))
     .        / (wjp*(w_k(c)-w_k(c)-wj))



      pmaa = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .        (-i*sin((-w_k(c)+w_k(c)+wj)*t3 - wjp*t2)
     .        + cos((-w_k(c)+w_k(c)+wj)*t3 - wjp*t2)
     .        + i*sin((-w_k(c)+w_k(c)+wj)*t3 - wjp*t1)
     .        - cos((-w_k(c)+w_k(c)+wj)*t3 - wjp*t1)
     .        - i*sin(wjp*t3 + (w_k(c)-w_k(c)-wj)*t2)
     .        - cos(wjp*t3 + (w_k(c)-w_k(c)-wj)*t2)
     .        + i*sin(wjp*t3 + (w_k(c)-w_k(c)-wj)*t1)
     .        + cos(wjp*t3 + (w_k(c)-w_k(c)-wj)*t1)
     .        - i*sin((-w_k(c)+w_k(c)+wj)*t2 - wjp*t1)
     .        + cos((-w_k(c)+w_k(c)+wj)*t2 - wjp*t1)
     .        - i*sin(wjp*t2 + (w_k(c)-w_k(c)-wj)*t1)
     .        - cos(wjp*t2 + (w_k(c)-w_k(c)-wj)*t1)
     .        + i*sin((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(c)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(c)+w_k(c)-wjp+wj)*t2))
     .        / (wjp*(w_k(c)-w_k(c)-wj))

      ppaa = ppaa
     .     + 
     .        5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .        (i*sin((w_k(c)+wjp)*t3 + (w_k(c)+wj)*t2)
     .        - cos((w_k(c)+wjp)*t3 + (w_k(c)+wj)*t2)
     .        - i*sin((w_k(c)+wjp)*t3 + (w_k(c)+wj)*t1)
     .        + cos((w_k(c)+wjp)*t3 + (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t3 + (w_k(c)+wjp)*t2)
     .        + cos((w_k(c)+wj)*t3 + (w_k(c)+wjp)*t2)
     .        + i*sin((w_k(c)+wj)*t3 + (w_k(c)+wjp)*t1)
     .        - cos((w_k(c)+wj)*t3 + (w_k(c)+wjp)*t1)
     .        + i*sin((w_k(c)+wjp)*t2 + (w_k(c)+wj)*t1)
     .        - cos((w_k(c)+wjp)*t2 + (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t2 + (w_k(c)+wjp)*t1)
     .        + cos((w_k(c)+wj)*t2 + (w_k(c)+wjp)*t1))
     .        / ((w_k(c)+wj)*(w_k(c)+wjp))

      ppac = ppac
     .     + 
     .        5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .        (-i*sin((wjp-w_k(c))*t3 + (w_k(c)+wj)*t2)
     .        + cos((wjp-w_k(c))*t3 + (w_k(c)+wj)*t2)
     .        + i*sin((wjp-w_k(c))*t3 + (w_k(c)+wj)*t1)
     .        - cos((wjp-w_k(c))*t3 + (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t3 + (wjp-w_k(c))*t2)
     .        - cos((w_k(c)+wj)*t3 + (wjp-w_k(c))*t2)
     .        - i*sin((w_k(c)+wj)*t3 + (wjp-w_k(c))*t1)
     .        + cos((w_k(c)+wj)*t3 + (wjp-w_k(c))*t1)
     .        - i*sin((wjp-w_k(c))*t2 + (w_k(c)+wj)*t1)
     .        + cos((wjp-w_k(c))*t2 + (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t2 + (wjp-w_k(c))*t1)
     .        - cos((w_k(c)+wj)*t2 + (wjp-w_k(c))*t1))
     .        / ((w_k(c)+wj)*(w_k(c)-wjp))


      pmac = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .        (-i*sin((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        - cos((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        + i*sin((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        + cos((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t2)
     .        + cos((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t2)
     .        + i*sin((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t1)
     .        - cos((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t1)
     .        - i*sin((w_k(c)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - cos((w_k(c)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t2 - (w_k(c)+wjp)*t1)
     .        + cos((w_k(c)+wj)*t2 - (w_k(c)+wjp)*t1)
     .        + i*sin((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(c)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(c)+w_k(c)-wjp+wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(c)+wjp))

      pmaa = pmaa
     .     +  5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .        (i*sin((wjp-w_k(c))*t3 - (w_k(c)+wj)*t2)
     .        + cos((wjp-w_k(c))*t3 - (w_k(c)+wj)*t2)
     .        - i*sin((wjp-w_k(c))*t3 - (w_k(c)+wj)*t1)
     .        - cos((wjp-w_k(c))*t3 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t2)
     .        - cos((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t2)
     .        - i*sin((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t1)
     .        + cos((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t1)
     .        + i*sin((wjp-w_k(c))*t2 - (w_k(c)+wj)*t1)
     .        + cos((wjp-w_k(c))*t2 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t2 + (w_k(c)-wjp)*t1)
     .        - cos((w_k(c)+wj)*t2 + (w_k(c)-wjp)*t1)
     .        - i*sin((w_k(c)+w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(c)+w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(c)-w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(c)-w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(c)-wjp))


      ppca = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .        (i*sin((w_k(c)+wjp)*t3 + (wj-w_k(c))*t2)
     .        - cos((w_k(c)+wjp)*t3 + (wj-w_k(c))*t2)
     .        - i*sin((w_k(c)+wjp)*t3 + (wj-w_k(c))*t1)
     .        + cos((w_k(c)+wjp)*t3 + (wj-w_k(c))*t1)
     .        - i*sin((wj-w_k(c))*t3 + (w_k(c)+wjp)*t2)
     .        + cos((wj-w_k(c))*t3 + (w_k(c)+wjp)*t2)
     .        + i*sin((wj-w_k(c))*t3 + (w_k(c)+wjp)*t1)
     .        - cos((wj-w_k(c))*t3 + (w_k(c)+wjp)*t1)
     .        + i*sin((w_k(c)+wjp)*t2 + (wj-w_k(c))*t1)
     .        - cos((w_k(c)+wjp)*t2 + (wj-w_k(c))*t1)
     .        - i*sin((wj-w_k(c))*t2 + (w_k(c)+wjp)*t1)
     .        + cos((wj-w_k(c))*t2 + (w_k(c)+wjp)*t1))
     .        / ((w_k(c)-wj)*(w_k(c)+wjp))

      pmcc = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .       (-i*sin((wjp-w_k(c))*t3 + (wj-w_k(c))*t2)
     .        + cos((wjp-w_k(c))*t3 + (wj-w_k(c))*t2)
     .        + i*sin((wjp-w_k(c))*t3 + (wj-w_k(c))*t1)
     .        - cos((wjp-w_k(c))*t3 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t3 + (wjp-w_k(c))*t2)
     .        - cos((wj-w_k(c))*t3 + (wjp-w_k(c))*t2)
     .        - i*sin((wj-w_k(c))*t3 + (wjp-w_k(c))*t1)
     .        + cos((wj-w_k(c))*t3 + (wjp-w_k(c))*t1)
     .        - i*sin((wjp-w_k(c))*t2 + (wj-w_k(c))*t1)
     .        + cos((wjp-w_k(c))*t2 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t2 + (wjp-w_k(c))*t1)
     .        - cos((wj-w_k(c))*t2 + (wjp-w_k(c))*t1))
     .        / ((w_k(c)-wj)*(w_k(c)-wjp))


      pmca = 5.d-1*RABI(a)*RABI(b)*LD(a,c)*LD(b,c)*
     .       -(i*sin((wjp-w_k(c))*t3 + (w_k(c)-wj)*t2)
     .        + cos((wjp-w_k(c))*t3 + (w_k(c)-wj)*t2)
     .        - i*sin((wjp-w_k(c))*t3 + (w_k(c)-wj)*t1)
     .        - cos((wjp-w_k(c))*t3 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t3 + (w_k(c)-wjp)*t2)
     .        - cos((wj-w_k(c))*t3 + (w_k(c)-wjp)*t2)
     .        - i*sin((wj-w_k(c))*t3 + (w_k(c)-wjp)*t1)
     .        + cos((wj-w_k(c))*t3 + (w_k(c)-wjp)*t1)
     .        + i*sin((wjp-w_k(c))*t2 + (w_k(c)-wj)*t1)
     .        + cos((wjp-w_k(c))*t2 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t2 + (w_k(c)-wjp)*t1)
     .        - cos((wj-w_k(c))*t2 + (w_k(c)-wjp)*t1)
     .        - i*sin((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(c)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(c)+w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)-wj)*(w_k(c)-wjp))

      


      H = (0.d0,0.d0)

C=======================================================================
C     pp block: rows 0:3, cols 12:15
C     qcc=ppcc, qca=ppca, qac=ppac, qaa=ppaa
C=======================================================================

      H(0,12) = ppac
      H(0,14) = sqrt(2.d0)*ppaa

      H(1,13) = ppca + 2.d0*ppac
      H(1,15) = sqrt(6.d0)*ppaa

      H(2,12) = sqrt(2.d0)*ppcc
      H(2,14) = 2.d0*ppca + 3.d0*ppac

      H(3,13) = sqrt(6.d0)*ppcc
      H(3,15) = 3.d0*ppca


C=======================================================================
C     pm block: rows 4:7, cols 8:11
C     qcc=pmcc, qca=pmca, qac=pmac, qaa=pmaa
C=======================================================================

      H(4,8)  = pmac
      H(4,10) = sqrt(2.d0)*pmaa

      H(5,9)  = pmca + 2.d0*pmac
      H(5,11) = sqrt(6.d0)*pmaa

      H(6,8)  = sqrt(2.d0)*pmcc
      H(6,10) = 2.d0*pmca + 3.d0*pmac

      H(7,9)  = sqrt(6.d0)*pmcc
      H(7,11) = 3.d0*pmca


C=======================================================================
C     mp block: rows 8:11, cols 4:7
C     Hermitian conjugate of pm block
C=======================================================================

      H(8,4)   = conjg(pmac)
      H(10,4)  = sqrt(2.d0)*conjg(pmaa)

      H(9,5)   = conjg(pmca + 2.d0*pmac)
      H(11,5)  = sqrt(6.d0)*conjg(pmaa)

      H(8,6)   = sqrt(2.d0)*conjg(pmcc)
      H(10,6)  = conjg(2.d0*pmca + 3.d0*pmac)

      H(9,7)   = sqrt(6.d0)*conjg(pmcc)
      H(11,7)  = 3.d0*conjg(pmca)


C=======================================================================
C     mm block: rows 12:15, cols 0:3
C     Hermitian conjugate of pp block
C=======================================================================

      H(12,0)  = conjg(ppac)
      H(14,0)  = sqrt(2.d0)*conjg(ppaa)

      H(13,1)  = conjg(ppca + 2.d0*ppac)
      H(15,1)  = sqrt(6.d0)*conjg(ppaa)

      H(12,2)  = sqrt(2.d0)*conjg(ppcc)
      H(14,2)  = conjg(2.d0*ppca + 3.d0*ppac)

      H(13,3)  = sqrt(6.d0)*conjg(ppcc)
      H(15,3)  = 3.d0*conjg(ppca)

      CALL EXACT_TAYLOR(15,H,H_C)

      call sss_operators(H_C,PSI_IN,j,k,a,b,c)

      return
      end

      subroutine sss_operators(H,PSI,j,k,ion1,ion2,mode)
      IMPLICIT NONE
      INTEGER j,k,a,b,c,d,ion1,ion2,mode
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,ip1,ip2,IND(0:15)
      COMPLEX*16 OUTP(0:15),H(0:15,0:15)
      COMPLEX*16 PSI(0:2**j*4**k-1)
      
      mp1 = k*2-mode*2
      mp2 = k*2-mode*2+1
      ip1 = k*2+j-ion1
      ip2 = k*2+j-ion2

      do ii = 0, 2**(j+2*k-4)-1
         call gen_other_bits_sss(j+2*k-1,mp1,mp2,ip1,ip2,ii,base_bits)
         
         idx = 0
         do a = 0,1
           do b = 0,1
            do c = 0,1
              do d = 0,1
                bit_combo = base_bits
                if (a .EQ. 1) bit_combo = IBSET(bit_combo, ip1)
                if (b .EQ. 1) bit_combo = IBSET(bit_combo, ip2)
                if (c .EQ. 1) bit_combo = IBSET(bit_combo, mp2)
                if (d .EQ. 1) bit_combo = IBSET(bit_combo, mp1)

                IND(idx) = bit_combo
                idx = idx + 1
              end do
             end do
           end do
         end do 

         do a = 0,15
            OUTP(a) = PSI(IND(a))
         end do
         
         OUTP = MATMUL(H,OUTP)
            

         do a = 0,15
            PSI(IND(a)) = OUTP(a)
         end do 
        
      end do

      return
      end

      subroutine gen_other_bits_sss(nbits,mp1,mp2,ip1,ip2,num,result)
      IMPLICIT NONE 
      INTEGER nbits,mp2,mp1,ip1,ip2,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2
     .      .AND. bit_pos .NE. ip1 .AND. bit_pos .NE. ip2) then
           if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
           end if
           source_bit = source_bit + 1
         end if
      end do
      
      return
      end



C=====================================================================
C                       IDENTITY X SINGLE
C=====================================================================

      subroutine identity_single(j,k,a,c,t3,t2,t1,PSI_IN,
     .                           LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,a,c
      DOUBLE PRECISION t1,t2,t3
      DOUBLE PRECISION LD(j,k),RABI(7),w(7),wqbt(7),w_k(7),wj,wjp
      COMPLEX*16 pmaa,pmac,pmca,pmcc,i
      COMPLEX*16 PSI_IN(0:2**j*4**k-1),H(0:7,0:7),H_C(0:7,0:7)
      
      i = (0.d0,1.d0)
      wj = w(a) - wqbt(a)
      wjp = wj

      pmca = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,c)*
     .        (i*sin((-w_k(c)+w_k(c)+wjp)*t3 - wj*t2)
     .        + cos((-w_k(c)+w_k(c)+wjp)*t3 - wj*t2)
     .        - i*sin((-w_k(c)+w_k(c)+wjp)*t3 - wj*t1)
     .        - cos((-w_k(c)+w_k(c)+wjp)*t3 - wj*t1)
     .        + i*sin(wj*t3 + (w_k(c)-w_k(c)-wjp)*t2)
     .        - cos(wj*t3 + (w_k(c)-w_k(c)-wjp)*t2)
     .        - i*sin(wj*t3 + (w_k(c)-w_k(c)-wjp)*t1)
     .        + cos(wj*t3 + (w_k(c)-w_k(c)-wjp)*t1)
     .        + i*sin((-w_k(c)+w_k(c)+wjp)*t2 - wj*t1)
     .        + cos((-w_k(c)+w_k(c)+wjp)*t2 - wj*t1)
     .        + i*sin(wj*t2 + (w_k(c)-w_k(c)-wjp)*t1)
     .        - cos(wj*t2 + (w_k(c)-w_k(c)-wjp)*t1)
     .        - i*sin((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(c)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(c)+w_k(c)+wjp-wj)*t2))
     .        / (wj*(w_k(c)-w_k(c)-wjp))
      


      pmaa = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,c)*
     .        (-i*sin((-w_k(c)+w_k(c)+wj)*t3 - wjp*t2)
     .        + cos((-w_k(c)+w_k(c)+wj)*t3 - wjp*t2)
     .        + i*sin((-w_k(c)+w_k(c)+wj)*t3 - wjp*t1)
     .        - cos((-w_k(c)+w_k(c)+wj)*t3 - wjp*t1)
     .        - i*sin(wjp*t3 + (w_k(c)-w_k(c)-wj)*t2)
     .        - cos(wjp*t3 + (w_k(c)-w_k(c)-wj)*t2)
     .        + i*sin(wjp*t3 + (w_k(c)-w_k(c)-wj)*t1)
     .        + cos(wjp*t3 + (w_k(c)-w_k(c)-wj)*t1)
     .        - i*sin((-w_k(c)+w_k(c)+wj)*t2 - wjp*t1)
     .        + cos((-w_k(c)+w_k(c)+wj)*t2 - wjp*t1)
     .        - i*sin(wjp*t2 + (w_k(c)-w_k(c)-wj)*t1)
     .        - cos(wjp*t2 + (w_k(c)-w_k(c)-wj)*t1)
     .        + i*sin((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(c)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(c)+w_k(c)-wjp+wj)*t2))
     .        / (wjp*(w_k(c)-w_k(c)-wj))

      pmac = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,c)*
     .        (-i*sin((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        - cos((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t2)
     .        + i*sin((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        + cos((w_k(c)+wjp)*t3 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t2)
     .        + cos((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t2)
     .        + i*sin((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t1)
     .        - cos((w_k(c)+wj)*t3 - (w_k(c)+wjp)*t1)
     .        - i*sin((w_k(c)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - cos((w_k(c)+wjp)*t2 - (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t2 - (w_k(c)+wjp)*t1)
     .        + cos((w_k(c)+wj)*t2 - (w_k(c)+wjp)*t1)
     .        + i*sin((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + cos((w_k(c)-w_k(c)+wjp-wj)*t2)
     .        + i*sin((-w_k(c)+w_k(c)-wjp+wj)*t2)
     .        - cos((-w_k(c)+w_k(c)-wjp+wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(c)+wjp))

      pmaa = pmaa
     .     +  5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,c)*
     .        (i*sin((wjp-w_k(c))*t3 - (w_k(c)+wj)*t2)
     .        + cos((wjp-w_k(c))*t3 - (w_k(c)+wj)*t2)
     .        - i*sin((wjp-w_k(c))*t3 - (w_k(c)+wj)*t1)
     .        - cos((wjp-w_k(c))*t3 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t2)
     .        - cos((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t2)
     .        - i*sin((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t1)
     .        + cos((w_k(c)+wj)*t3 + (w_k(c)-wjp)*t1)
     .        + i*sin((wjp-w_k(c))*t2 - (w_k(c)+wj)*t1)
     .        + cos((wjp-w_k(c))*t2 - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t2 + (w_k(c)-wjp)*t1)
     .        - cos((w_k(c)+wj)*t2 + (w_k(c)-wjp)*t1)
     .        - i*sin((w_k(c)+w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(c)+w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(c)-w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(c)-w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(c)-wjp))

      pmcc = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,c)*
     .       (-i*sin((wjp-w_k(c))*t3 + (wj-w_k(c))*t2)
     .        + cos((wjp-w_k(c))*t3 + (wj-w_k(c))*t2)
     .        + i*sin((wjp-w_k(c))*t3 + (wj-w_k(c))*t1)
     .        - cos((wjp-w_k(c))*t3 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t3 + (wjp-w_k(c))*t2)
     .        - cos((wj-w_k(c))*t3 + (wjp-w_k(c))*t2)
     .        - i*sin((wj-w_k(c))*t3 + (wjp-w_k(c))*t1)
     .        + cos((wj-w_k(c))*t3 + (wjp-w_k(c))*t1)
     .        - i*sin((wjp-w_k(c))*t2 + (wj-w_k(c))*t1)
     .        + cos((wjp-w_k(c))*t2 + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t2 + (wjp-w_k(c))*t1)
     .        - cos((wj-w_k(c))*t2 + (wjp-w_k(c))*t1))
     .        / ((w_k(c)-wj)*(w_k(c)-wjp))


      pmca = 5.d-1*RABI(a)*RABI(a)*LD(a,c)*LD(a,c)*
     .       -(i*sin((wjp-w_k(c))*t3 + (w_k(c)-wj)*t2)
     .        + cos((wjp-w_k(c))*t3 + (w_k(c)-wj)*t2)
     .        - i*sin((wjp-w_k(c))*t3 + (w_k(c)-wj)*t1)
     .        - cos((wjp-w_k(c))*t3 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t3 + (w_k(c)-wjp)*t2)
     .        - cos((wj-w_k(c))*t3 + (w_k(c)-wjp)*t2)
     .        - i*sin((wj-w_k(c))*t3 + (w_k(c)-wjp)*t1)
     .        + cos((wj-w_k(c))*t3 + (w_k(c)-wjp)*t1)
     .        + i*sin((wjp-w_k(c))*t2 + (w_k(c)-wj)*t1)
     .        + cos((wjp-w_k(c))*t2 + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t2 + (w_k(c)-wjp)*t1)
     .        - cos((wj-w_k(c))*t2 + (w_k(c)-wjp)*t1)
     .        - i*sin((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        + cos((w_k(c)-w_k(c)-wjp+wj)*t2)
     .        - i*sin((-w_k(c)+w_k(c)+wjp-wj)*t2)
     .        - cos((-w_k(c)+w_k(c)+wjp-wj)*t2))
     .        / ((w_k(c)-wj)*(w_k(c)-wjp))


      
      H = (0.d0,0.d0)

C=======================================================================
C     top-left block: rows 0:3, cols 0:3
C     qcc=pmcc, qca=pmca, qac=pmac, qaa=pmaa
C=======================================================================

      H(0,0) = pmac
      H(0,2) = sqrt(2.d0)*pmaa

      H(1,1) = pmca + 2.d0*pmac
      H(1,3) = sqrt(6.d0)*pmaa

      H(2,0) = sqrt(2.d0)*pmcc
      H(2,2) = 2.d0*pmca + 3.d0*pmac

      H(3,1) = sqrt(6.d0)*pmcc
      H(3,3) = 3.d0*pmca


C=======================================================================
C     bottom-right block: rows 4:7, cols 4:7
C     conjugate of top-left block
C=======================================================================

      H(4,4) = conjg(pmac)
      H(4,6) = sqrt(2.d0)*conjg(pmaa)

      H(5,5) = conjg(pmca + 2.d0*pmac)
      H(5,7) = sqrt(6.d0)*conjg(pmaa)

      H(6,4) = sqrt(2.d0)*conjg(pmcc)
      H(6,6) = conjg(2.d0*pmca + 3.d0*pmac)

      H(7,5) = sqrt(6.d0)*conjg(pmcc)
      H(7,7) = 3.d0*conjg(pmca)
      
      CALL EXACT_TAYLOR(7,H,H_C)
      CALL idmode_operators(H_C,PSI_IN,j,k,a,c)

      return
      end


C=====================================================================
C                   SPIN x SPIN x MODE x MODE x MODE
C=====================================================================


      subroutine spin_spin_mode_mode_mode(j,k,a,b,c,d,e,t3,t2,t1,
     .                                    PSI_IN,LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,a,b,c,d,e,s1,s2,n1,n2,n3,col,row
      DOUBLE PRECISION LD(j,k),RABI(7),w(7),wqbt(7),w_k(7),amp,wj
      DOUBLE PRECISION t3,t2,t1,wjp
      COMPLEX*16 ppaac,pmaca,ppcac,pmcca,ppaca,ppacc,pmacc,i
      COMPLEX*16 H(0:255,0:255),H_C(0:255,0:255),PSI_IN(0:2**j*4**k-1)


      i = (0.d0,1.d0)
      wj = w(a) - wqbt(a)
      wjp = w(b) - wqbt(b)

      ppaac = -(i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)+wj)*t2)
     .        - cos((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)+wj)*t2)
     .        - i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)+wj)*t1)
     .        + cos((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t2)
     .        + cos((w_k(c)+wj)*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t2)
     .        + i*sin((w_k(c)+wj)*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t1)
     .        - cos((w_k(c)+wj)*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t1)
     .        + i*sin((-w_k(e)+w_k(d)+wjp)*t2
     .        + (w_k(c)+wj)*t1)
     .        - cos((-w_k(e)+w_k(d)+wjp)*t2
     .        + (w_k(c)+wj)*t1)
     .        - i*sin((w_k(c)+wj)*t2
     .        + (-w_k(e)+w_k(d)+wjp)*t1)
     .        + cos((w_k(c)+wj)*t2
     .        + (-w_k(e)+w_k(d)+wjp)*t1))
     .        / ((w_k(c)+wj)*(w_k(e)-w_k(d)-wjp))

      pmaca = (i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        - (w_k(c)+wj)*t2)
     .        + cos((-w_k(e)+w_k(d)+wjp)*t3
     .        - (w_k(c)+wj)*t2)
     .        - i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        - (w_k(c)+wj)*t1)
     .        - cos((-w_k(e)+w_k(d)+wjp)*t3
     .        - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t3
     .        + (w_k(e)-w_k(d)-wjp)*t2)
     .        - cos((w_k(c)+wj)*t3
     .        + (w_k(e)-w_k(d)-wjp)*t2)
     .        - i*sin((w_k(c)+wj)*t3
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        + cos((w_k(c)+wj)*t3
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        + i*sin((-w_k(e)+w_k(d)+wjp)*t2
     .        - (w_k(c)+wj)*t1)
     .        + cos((-w_k(e)+w_k(d)+wjp)*t2
     .        - (w_k(c)+wj)*t1)
     .        + i*sin((w_k(c)+wj)*t2
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        - cos((w_k(c)+wj)*t2
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        - i*sin((w_k(e)-w_k(d)-wjp+w_k(c)+wj)*t2)
     .        + cos((w_k(e)-w_k(d)-wjp+w_k(c)+wj)*t2)
     .        - i*sin((-w_k(e)+w_k(d)+wjp-w_k(c)-wj)*t2)
     .        - cos((-w_k(e)+w_k(d)+wjp-w_k(c)-wj)*t2))
     .        / ((w_k(c)+wj)*(w_k(e)-w_k(d)-wjp))

      ppcac = -(-i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        + (wj-w_k(c))*t2)
     .        + cos((-w_k(e)+w_k(d)+wjp)*t3
     .        + (wj-w_k(c))*t2)
     .        + i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        + (wj-w_k(c))*t1)
     .        - cos((-w_k(e)+w_k(d)+wjp)*t3
     .        + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t2)
     .        - cos((wj-w_k(c))*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t2)
     .        - i*sin((wj-w_k(c))*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t1)
     .        + cos((wj-w_k(c))*t3
     .        + (-w_k(e)+w_k(d)+wjp)*t1)
     .        - i*sin((-w_k(e)+w_k(d)+wjp)*t2
     .        + (wj-w_k(c))*t1)
     .        + cos((-w_k(e)+w_k(d)+wjp)*t2
     .        + (wj-w_k(c))*t1)
     .        + i*sin((wj-w_k(c))*t2
     .        + (-w_k(e)+w_k(d)+wjp)*t1)
     .        - cos((wj-w_k(c))*t2
     .        + (-w_k(e)+w_k(d)+wjp)*t1))
     .        / ((w_k(c)-wj)*(w_k(e)-w_k(d)-wjp))


      
      pmcca = -(i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)-wj)*t2)
     .        + cos((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)-wj)*t2)
     .        - i*sin((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)-wj)*t1)
     .        - cos((-w_k(e)+w_k(d)+wjp)*t3
     .        + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t3
     .        + (w_k(e)-w_k(d)-wjp)*t2)
     .        - cos((wj-w_k(c))*t3
     .        + (w_k(e)-w_k(d)-wjp)*t2)
     .        - i*sin((wj-w_k(c))*t3
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        + cos((wj-w_k(c))*t3
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        + i*sin((-w_k(e)+w_k(d)+wjp)*t2
     .        + (w_k(c)-wj)*t1)
     .        + cos((-w_k(e)+w_k(d)+wjp)*t2
     .        + (w_k(c)-wj)*t1)
     .        + i*sin((wj-w_k(c))*t2
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        - cos((wj-w_k(c))*t2
     .        + (w_k(e)-w_k(d)-wjp)*t1)
     .        - i*sin((w_k(e)-w_k(d)-wjp-w_k(c)+wj)*t2)
     .        + cos((w_k(e)-w_k(d)-wjp-w_k(c)+wj)*t2)
     .        - i*sin((-w_k(e)+w_k(d)+wjp+w_k(c)-wj)*t2)
     .        - cos((-w_k(e)+w_k(d)+wjp+w_k(c)-wj)*t2))
     .        / ((w_k(c)-wj)*(w_k(e)-w_k(d)-wjp))

      ppaca = -(i*sin((w_k(e)+wjp)*t3
     .        + (-w_k(d)+w_k(c)+wj)*t2)
     .        - cos((w_k(e)+wjp)*t3
     .        + (-w_k(d)+w_k(c)+wj)*t2)
     .        - i*sin((w_k(e)+wjp)*t3
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        + cos((w_k(e)+wjp)*t3
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)+wjp)*t2)
     .        + cos((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)+wjp)*t2)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)+wjp)*t1)
     .        - cos((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)+wjp)*t1)
     .        + i*sin((w_k(e)+wjp)*t2
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        - cos((w_k(e)+wjp)*t2
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t2
     .        + (w_k(e)+wjp)*t1)
     .        + cos((-w_k(d)+w_k(c)+wj)*t2
     .        + (w_k(e)+wjp)*t1))
     .        / ((w_k(d)-w_k(c)-wj)*(w_k(e)+wjp))

      ppacc = -(-i*sin((wjp-w_k(e))*t3
     .        + (-w_k(d)+w_k(c)+wj)*t2)
     .        + cos((wjp-w_k(e))*t3
     .        + (-w_k(d)+w_k(c)+wj)*t2)
     .        + i*sin((wjp-w_k(e))*t3
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        - cos((wjp-w_k(e))*t3
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        + (wjp-w_k(e))*t2)
     .        - cos((-w_k(d)+w_k(c)+wj)*t3
     .        + (wjp-w_k(e))*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        + (wjp-w_k(e))*t1)
     .        + cos((-w_k(d)+w_k(c)+wj)*t3
     .        + (wjp-w_k(e))*t1)
     .        - i*sin((wjp-w_k(e))*t2
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        + cos((wjp-w_k(e))*t2
     .        + (-w_k(d)+w_k(c)+wj)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t2
     .        + (wjp-w_k(e))*t1)
     .        - cos((-w_k(d)+w_k(c)+wj)*t2
     .        + (wjp-w_k(e))*t1))
     .        / ((w_k(d)-w_k(c)-wj)*(w_k(e)-wjp))

      pmacc = (i*sin((w_k(e)+wjp)*t3
     .        + (w_k(d)-w_k(c)-wj)*t2)
     .        + cos((w_k(e)+wjp)*t3
     .        + (w_k(d)-w_k(c)-wj)*t2)
     .        - i*sin((w_k(e)+wjp)*t3
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        - cos((w_k(e)+wjp)*t3
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        - (w_k(e)+wjp)*t2)
     .        - cos((-w_k(d)+w_k(c)+wj)*t3
     .        - (w_k(e)+wjp)*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        - (w_k(e)+wjp)*t1)
     .        + cos((-w_k(d)+w_k(c)+wj)*t3
     .        - (w_k(e)+wjp)*t1)
     .        + i*sin((w_k(e)+wjp)*t2
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        + cos((w_k(e)+wjp)*t2
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t2
     .        - (w_k(e)+wjp)*t1)
     .        - cos((-w_k(d)+w_k(c)+wj)*t2
     .        - (w_k(e)+wjp)*t1)
     .        - i*sin((w_k(e)+w_k(d)+wjp-w_k(c)-wj)*t2)
     .        - cos((w_k(e)+w_k(d)+wjp-w_k(c)-wj)*t2)
     .        - i*sin((-w_k(e)-w_k(d)-wjp+w_k(c)+wj)*t2)
     .        + cos((-w_k(e)-w_k(d)-wjp+w_k(c)+wj)*t2))
     .        / ((w_k(d)-w_k(c)-wj)*(w_k(e)+wjp))

      pmaca = pmaca
     .     - (i*sin((wjp-w_k(e))*t3
     .        + (w_k(d)-w_k(c)-wj)*t2)
     .        + cos((wjp-w_k(e))*t3
     .        + (w_k(d)-w_k(c)-wj)*t2)
     .        - i*sin((wjp-w_k(e))*t3
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        - cos((wjp-w_k(e))*t3
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)-wjp)*t2)
     .        - cos((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)-wjp)*t2)
     .        - i*sin((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)-wjp)*t1)
     .        + cos((-w_k(d)+w_k(c)+wj)*t3
     .        + (w_k(e)-wjp)*t1)
     .        + i*sin((wjp-w_k(e))*t2
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        + cos((wjp-w_k(e))*t2
     .        + (w_k(d)-w_k(c)-wj)*t1)
     .        + i*sin((-w_k(d)+w_k(c)+wj)*t2
     .        + (w_k(e)-wjp)*t1)
     .        - cos((-w_k(d)+w_k(c)+wj)*t2
     .        + (w_k(e)-wjp)*t1)
     .        - i*sin((w_k(e)-w_k(d)-wjp+w_k(c)+wj)*t2)
     .        + cos((w_k(e)-w_k(d)-wjp+w_k(c)+wj)*t2)
     .        - i*sin((-w_k(e)+w_k(d)+wjp-w_k(c)-wj)*t2)
     .        - cos((-w_k(e)+w_k(d)+wjp-w_k(c)-wj)*t2))
     .        / ((w_k(d)-w_k(c)-wj)*(w_k(e)-wjp))

      H = (0.d0,0.d0)

      do s1 = 0,1
      do s2 = 0,1
      do n1 = 0,3
      do n2 = 0,3
      do n3 = 0,3

         col = s1*128 + s2*64 + n1*16 + n2*4 + n3

C=======================================================================
C        ppaac = p x p x a x a x c
C=======================================================================

         if (s1 .eq. 1 .and. s2 .eq. 1 .and.
     .       n1 .ge. 1 .and. n2 .ge. 1 .and. n3 .le. 2) then

            row = (n1-1)*16 + (n2-1)*4 + (n3+1)

            amp = sqrt(dble(n1))
     .           *sqrt(dble(n2))
     .           *sqrt(dble(n3+1))

            H(row,col) = H(row,col) + amp*ppaac
            H(col,row) = H(col,row) + amp*conjg(ppaac)

         endif


C=======================================================================
C        pmaca = p x m x a x c x a
C=======================================================================

         if (s1 .eq. 1 .and. s2 .eq. 0 .and.
     .       n1 .ge. 1 .and. n2 .le. 2 .and. n3 .ge. 1) then

            row = 64 + (n1-1)*16 + (n2+1)*4 + (n3-1)

            amp = sqrt(dble(n1))
     .           *sqrt(dble(n2+1))
     .           *sqrt(dble(n3))

            H(row,col) = H(row,col) + amp*pmaca
            H(col,row) = H(col,row) + amp*conjg(pmaca)

         endif


C=======================================================================
C        ppcac = p x p x c x a x c
C=======================================================================

         if (s1 .eq. 1 .and. s2 .eq. 1 .and.
     .       n1 .le. 2 .and. n2 .ge. 1 .and. n3 .le. 2) then

            row = (n1+1)*16 + (n2-1)*4 + (n3+1)

            amp = sqrt(dble(n1+1))
     .           *sqrt(dble(n2))
     .           *sqrt(dble(n3+1))

            H(row,col) = H(row,col) + amp*ppcac
            H(col,row) = H(col,row) + amp*conjg(ppcac)

         endif


C=======================================================================
C        pmcca = p x m x c x c x a
C=======================================================================

         if (s1 .eq. 1 .and. s2 .eq. 0 .and.
     .       n1 .le. 2 .and. n2 .le. 2 .and. n3 .ge. 1) then

            row = 64 + (n1+1)*16 + (n2+1)*4 + (n3-1)

            amp = sqrt(dble(n1+1))
     .           *sqrt(dble(n2+1))
     .           *sqrt(dble(n3))

            H(row,col) = H(row,col) + amp*pmcca
            H(col,row) = H(col,row) + amp*conjg(pmcca)

         endif


C=======================================================================
C        ppaca = p x p x a x c x a
C=======================================================================

         if (s1 .eq. 1 .and. s2 .eq. 1 .and.
     .       n1 .ge. 1 .and. n2 .le. 2 .and. n3 .ge. 1) then

            row = (n1-1)*16 + (n2+1)*4 + (n3-1)

            amp = sqrt(dble(n1))
     .           *sqrt(dble(n2+1))
     .           *sqrt(dble(n3))

            H(row,col) = H(row,col) + amp*ppaca
            H(col,row) = H(col,row) + amp*conjg(ppaca)

         endif


C=======================================================================
C        ppacc = p x p x a x c x c
C=======================================================================

         if (s1 .eq. 1 .and. s2 .eq. 1 .and.
     .       n1 .ge. 1 .and. n2 .le. 2 .and. n3 .le. 2) then

            row = (n1-1)*16 + (n2+1)*4 + (n3+1)

            amp = sqrt(dble(n1))
     .           *sqrt(dble(n2+1))
     .           *sqrt(dble(n3+1))

            H(row,col) = H(row,col) + amp*ppacc
            H(col,row) = H(col,row) + amp*conjg(ppacc)

         endif


C=======================================================================
C        pmacc = p x m x a x c x c
C=======================================================================

         if (s1 .eq. 1 .and. s2 .eq. 0 .and.
     .       n1 .ge. 1 .and. n2 .le. 2 .and. n3 .le. 2) then

            row = 64 + (n1-1)*16 + (n2+1)*4 + (n3+1)

            amp = sqrt(dble(n1))
     .           *sqrt(dble(n2+1))
     .           *sqrt(dble(n3+1))

            H(row,col) = H(row,col) + amp*pmacc
            H(col,row) = H(col,row) + amp*conjg(pmacc)

         endif

      enddo
      enddo
      enddo
      enddo
      enddo

      CALL EXACT_TAYLOR(255,H,H_C)
      
      CALL ssmmm_operators(H_C,PSI_IN,j,k,a,b,c,d,e)

      return
      end

      subroutine ssmmm_operators(H,PSI,j,k,ion1,ion2,
     .                          mode1,mode2,mode3)
      IMPLICIT NONE
      INTEGER a,b,c,d,e,f,g,i,j,k,ion1,ion2,mode1,mode2,mode3
      INTEGER ii,base_bits,idx,bit_combo
      INTEGER mp1,mp2,mp3,mp4,mp5,mp6,ip1,ip2
      INTEGER IND(0:255)
      COMPLEX*16 OUTP(0:255),H(0:255,0:255),PSI(0:2**j*4**k-1)

      mp1 = k*2-mode1*2
      mp2 = k*2-mode1*2+1
      mp3 = k*2-mode2*2
      mp4 = k*2-mode2*2+1
      mp5 = k*2-mode3*2
      mp6 = k*2-mode3*2+1

      ip1 = k*2+j-ion1
      ip2 = k*2+j-ion2

      do ii = 0, 2**(j+2*k-3)-1
         call gen_other_bits_ssmmm(j+2*k-1,mp1,mp2,mp3,mp4,mp5,mp6,
     .                             ip1,ip2,ii,base_bits)
            
         do a = 0,1
           do b = 0,1
             do c = 0,1
               do d = 0,1
                 do e = 0,1
                   do f = 0,1
                     do g = 0,1
                       do i = 0,1
                       bit_combo = base_bits
                       if (a .EQ. 1) bit_combo = IBSET(bit_combo, ip1)
                       if (b .EQ. 1) bit_combo = IBSET(bit_combo, ip2)
                       if (c .EQ. 1) bit_combo = IBSET(bit_combo, mp6)     
                       if (d .EQ. 1) bit_combo = IBSET(bit_combo, mp5) 
                       if (e .EQ. 1) bit_combo = IBSET(bit_combo, mp4) 
                       if (f .EQ. 1) bit_combo = IBSET(bit_combo, mp3)
                       if (g .EQ. 1) bit_combo = IBSET(bit_combo, mp2)
                       if (i .EQ. 1) bit_combo = IBSET(bit_combo, mp1)

                       IND(idx) = bit_combo
                       idx = idx + 1 
                       end do
                     end do          
                   end do     
                 end do 
               end do
             end do
           end do
         end do
      
        do a = 0,255
          OUTP(a) = PSI(IND(a))
        end do   
      
        OUTP = MATMUL(H,OUTP)

        do a = 0,255
          PSI(IND(a)) = OUTP(a)
        end do

      end do
      
      return
      end


      subroutine gen_other_bits_ssmmm(nbits,mp1,mp2,mp3,mp4,mp5,mp6,
     .                                ip1,ip2,num,result)
      IMPLICIT NONE
      INTEGER nbits,mp1,mp2,mp3,mp4,mp5,mp6,ip1,ip2,num,result
      INTEGER bit_pos,source_bit,temp

      result = 0
      source_bit = 0
      temp = num

      do bit_pos = 0, nbits
         if (bit_pos .NE. mp1 .AND. bit_pos .NE. mp2
     .       .AND. bit_pos .NE. mp3 .AND. bit_pos .NE. mp4
     .       .AND. bit_pos .NE. mp5 .AND. bit_pos .NE. mp6 
     .       .AND. bit_pos .NE. ip1 .AND. bit_pos .NE. ip2)
     .      then
            if (BTEST(temp, source_bit)) then
               result = IBSET(result, bit_pos)
            end if
            source_bit = source_bit + 1
          end if
      end do


      return
      end

C=====================================================================
C                      INTEGER x MODE x MODE x MODE
C=====================================================================

c      subroutine integer_mode_mode_mode(j,k,a,c,d,e,t3,t2,t1,PSI_IN,
c     .                                  LD,RABI,w,wqbt,w_k)
c      INTEGER j,k,a,c,d,e
c      DOUBLE PRECISION t3,t2,t1,RABI(7),LD(j,k),w(7),wqbt(7),w_k(7)
c      COMPLEX*16 pm
      

C=====================================================================

C=====================================================================

 
c==== create psi input vector      
      subroutine gen_psi(j,k,SPIN,MODE,VEC)
      IMPLICIT NONE
      INTEGER a,b,c,d,e,f,g,h,i,j,k,x,y,z,t,length
      COMPLEX*16 VEC(0:2**j*4**k-1)
      DOUBLE PRECISION SPIN(j,2), MODE(k,4)
      length = 2**j*4**k

      VEC = 1.d0
      do a = 1, k
         e = 1
         do b = 1, length/(4**(a-1))
            if (b .EQ. 1) then
                c = 1
            else
               c = 4**(a-1)*(b-1)+1
            end if

            do d = c, c + 4**(a-1)-1
               VEC(d-1) = VEC(d-1)*MODE(a,e)
            end do
            e = e + 1
            if (e .GT. 4) then
               e = 1
            end if
         end do
      end do

      do f = 1, j
         z = 1
         do g = 1, length/(4**k*2**(f-1))
            if (g .EQ. 1) then
               x = 1
            else
               x = 4**k*2**(f-1)*(g-1)+1
            end if

            do y = x, x + (4**k)*2**(f-1)-1
               VEC(y-1) = VEC(y-1)*SPIN(j-(f-1),z)
            end do
            z = z + 1
            if (z .GT. 2) then
               z = 1
            end if
         end do
      end do

      return
      end
       
      

c==== find exact taylor series using scaling trick
      SUBROUTINE EXACT_TAYLOR(size,MAT,CP)
      IMPLICIT NONE
      INTEGER size
      COMPLEX*16 MAT(0:size,0:size), CP(0:size,0:size), 
     .           MAT_REP(0:size,0:size)

      MAT_REP = MAT

c---- compute scaled exponential approximation into CP
      CALL TAYLOR_APPROXIMATION(size, MAT_REP, 50, -(0.D0,1.D0), CP)

      RETURN
      END

c==== taylor series with convergence check
      SUBROUTINE TAYLOR_APPROXIMATION(size,MAT,pow,co,APROX)
      IMPLICIT NONE
      INTEGER size, pow, i, j, k
      COMPLEX*16 MAT(0:size,0:size), APROX(0:size,0:size)
      COMPLEX*16 RES(0:size,0:size)
      COMPLEX*16 co, last
      DOUBLE PRECISION tol, norm_term, fact

      APROX = 0.d0
      RES = 0.d0
      fact = 1.d0
      last = 1.d0
      tol = 1.d-10

      do i = 0, size
         RES(i,i) = 1.d0
      end do

      APROX = RES
      
      do k = 1, pow
         
         RES = matmul(MAT,RES)
         last = last * co/fact
         APROX = APROX + last*RES
         fact = fact + 1.d0 

c------- compute frobenius norm of APROX
         norm_term = 0.D0
         do i = 0, size
            do j = 0, size
               norm_term=norm_term+
     .         abs(RES(i,j))**2
            end do
         end do
         norm_term = SQRT(norm_term)

c------- convergence check
         if (abs(last)*norm_term .lt. tol) then
            exit
         end if
      end do

      RETURN
      END

      SUBROUTINE check_hermitian(A,N)
      IMPLICIT NONE
      INTEGER N, i, j, notherm
      COMPLEX*16 A(0:N-1,0:N-1)
      DOUBLE PRECISION tol, diff

      tol = 1.0D-12
      notherm = 0

      DO i = 0, N-1
         DO j = 0, N-1
            diff = ABS(A(i,j) - CONJG(A(j,i)))
            print *, diff
            IF (diff .GT. tol) THEN
               notherm = notherm + 1
               PRINT *, 'Non-Herm at (', i, ',', j, ')'
            END IF
         END DO
      END DO

      IF (notherm .EQ. 0) THEN
         PRINT *, 'Matrix is Hermitian'
      ELSE
         PRINT *, 'Matrix not Hermitian, count=', notherm
      END IF
      
      RETURN
      END

      SUBROUTINE CHECK_UNITARY(U,N)
      IMPLICIT NONE
      INTEGER N
      COMPLEX*16 U(0:N,0:N), UDAG(0:N,0:N), PROD(0:N,0:N)
      DOUBLE PRECISION TOL, ERR
      INTEGER I, J

      TOL = 1.0D-12

C     Conjugate transpose
      DO I = 0, N
         DO J = 0, N
            UDAG(I,J) = CONJG(U(J,I))
         END DO
      END DO

C     Matrix product U† U
      PROD = MATMUL(UDAG, U)

C     Check deviation from identity
      ERR = 0.0D0
      DO I = 0, N
         DO J = 0, N
            IF (I .EQ. J) THEN
               write(6,*) I,J,PROD(I,J)
               ERR = MAX(ERR, ABS(PROD(I,J) - (1.0D0,0.0D0)))
            ELSE
               write(6,*) I,J,PROD(I,J)
               ERR = MAX(ERR, ABS(PROD(I,J)))
            END IF
         END DO
      END DO

      IF (ERR .LT. TOL) THEN
         PRINT *, 'Matrix is unitary. Error =', ERR
      ELSE
         PRINT *, 'Matrix is NOT unitary. Error =', ERR
      END IF

      END

      subroutine find_prob(PSI,j,k,ion,spin,PROB)
      IMPLICIT NONE
      INTEGER j,k,n,n1,n2,i,bi,n0,ion,spin
      COMPLEX*16 PSI(0:2**j*4**k-1)
      DOUBLE PRECISION PROB
      
      bi = spin
      n0 = bi*2**(2*k+(j-ion))
      PROB = 0.d0

      do i = 0, 2**(j+2*k-1)-1
         n = n0
         n1 = mod(i,2**(2*k+(j-ion)))
         n2 = i/(2**(j+2*k-1-(ion-1)))
         n = n+n1+n2*2**(j+2*k-(ion-1))
         PROB = PROB + cdabs(PSI(n))**2
      end do
      
      return
      end          
      

      subroutine sum_diff_sq(j,k,n,sum,MAG,l,MAG_C)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER n,l,j,k
      DOUBLE PRECISION MAG(0:n),MAG_C(0:j*k-1,0:n),sum
  
      sum = 0.d0
      
      do i = 0,n
         sum = sum + abs(MAG(i)-MAG_C(l,i))**2
      end do   

      return
      end

      subroutine find_mode_prob(PSI,j,k,mode,phonon,PROB)
      IMPLICIT NONE
      INTEGER bm0,bm1,n0,n1,n2,n,i,j,k,mode,phonon
      DOUBLE PRECISION PROB
      COMPLEX*16 PSI(0:2**j*4**k-1)
      
      bm0 = ibits(phonon,0,1)
      bm1 = ibits(phonon,1,1)

      n0 = bm0*2**(2*(k-mode))+bm1*2**(2*(k-mode)+1)

      PROB = 0.d0

      do i = 0, 2**(j+2*k-2)-1
         n = n0
         n1 = mod(i,2**(2*(k-mode)))
         n2 = i/(2**(2*(k-mode))) 
         n = n+n1+n2*2**(2*(k-mode+1))
         PROB = PROB + cdabs(PSI(n))**2
      end do
      return
      end

      subroutine find_joint_prob(VEC,j,k,ion,spin,mode,phonon,PROB)
      IMPLICIT NONE
      INTEGER bm0,bm1,bi,n0,n1,n2,n3,n,i,j,k,ion,mode,phonon,spin
      DOUBLE PRECISION PROB
      COMPLEX*16 VEC(0:2**j*4**k-1)

      bm0 = ibits(phonon,0,1)
      bm1 = ibits(phonon,1,1)
      bi = spin

      n0 = bm0*2**(2*(k-mode))+bm1*2**(2*(k-mode)+1)+
     .     bi*2**(2*k+(j-ion))

      PROB = 0.d0

      do i = 0, 2**(j+2*k-3)-1
         n = n0
         n1 = mod(i,2**(2*(k-mode)))
         n3 = i/(2**(j+2*k-3-(ion-1)))
         n2 = ibits(i,2*(k-mode),j+2*k-3-(ion-1)-(2*(k-mode)))

         n = n+n1+n3*2**(j+2*k-(ion-1))+n2*2**(2*(k-mode)+2)
         PROB = PROB + cdabs(VEC(n))**2

      end do

      return
      end
