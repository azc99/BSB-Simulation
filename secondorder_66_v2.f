      program neighbor_modes
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER config_num,data_interval,num_points
      CHARACTER*30 ARG
      DOUBLE PRECISION LD(6,6)
      DOUBLE PRECISION tau,del_t,take_data
      INTEGER start,end,count_rate
      DOUBLE PRECISION elapsed_time
      
      
c---- Get configuration number from command line
      call getarg(1, ARG)
      read(ARG, *) config_num
      
      write(6,*) 'Config', config_num+1, 'starting...'

      LD(1,1) =  0.0453
      LD(1,2) = -0.0658
      LD(1,3) = -0.0592
      LD(1,4) =  0.0395
      LD(1,5) = -0.0192
      LD(1,6) = -0.00603

      LD(2,1) =  0.0439
      LD(2,2) = -0.0375
      LD(2,3) =  0.0146
      LD(2,4) = -0.0591
      LD(2,5) =  0.0624
      LD(2,6) =  0.0338

      LD(3,1) =  0.0432
      LD(3,2) = -0.0123
      LD(3,3) =  0.0473
      LD(3,4) = -0.0314
      LD(3,5) = -0.0433
      LD(3,6) = -0.0711

      LD(4,1) =  0.0432
      LD(4,2) =  0.0122
      LD(4,3) =  0.0473
      LD(4,4) =  0.0314
      LD(4,5) = -0.0432
      LD(4,6) =  0.0711

      LD(5,1) =  0.0439
      LD(5,2) =  0.0375
      LD(5,3) =  0.0146
      LD(5,4) =  0.0591
      LD(5,5) =  0.0624
      LD(5,6) = -0.0338

      LD(6,1) =  0.0453
      LD(6,2) =  0.0658
      LD(6,3) = -0.0592
      LD(6,4) = -0.0395
      LD(6,5) = -0.0192
      LD(6,6) = 0.00603

      tau = 800.d-6
      del_t = 1.d-6
      take_data = 25.d-6
      data_interval = NINT(take_data/del_t)
      num_points = NINT(tau/take_data) + 1

      call system_clock(start,count_rate)
      call specific_time(6,6,tau,config_num,LD,del_t,data_interval,
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

      detune = 1.d-5

      RABI(1) = 1/(11.32)
      RABI(2) = 1/(9.38)
      RABI(3) = 1/(9.59)
      RABI(4) = 1/(9.88)
      RABI(5) = 1/(8.77)
      RABI(6) = 1/(ABI(

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

         do x = 1,j
           CALL spin_only(j,k,PSI_IN,x,m,t2,t1,LD,RABI,w,wqbt,w_k)
         enddo

         do x = 1,j
            CALL sigma_z(j,k,x,t2,t1,PSI_IN,RABI,w,wqbt,w_k)
         end do

         do x = 1,j
            
c            driving_mode = MOD((k-m+x-1),k)+1
c            low = MAX(driving_mode-1,1)
c            high = MIN(driving_mode+1,k)
            do y = 1,k
              CALL spin_mode(j,k,x,y,t2,t1,PSI_IN,LD,RABI,w,wqbt,w_k)
            end do
         end do

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
