      program neighbor_modes
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER config_num,data_interval,num_points
      CHARACTER*30 ARG
      DOUBLE PRECISION LD(7,7)
      DOUBLE PRECISION tau,del_t,take_data
      INTEGER start,end,count_rate
      DOUBLE PRECISION elapsed_time
      
      
c---- Get configuration number from command line
      call getarg(1, ARG)
      read(ARG, *) config_num
      
      write(6,*) 'Config', config_num+1, 'starting...'

      LD(1,1) =  0.0440d0
      LD(1,2) =  0.0627d0
      LD(1,3) = -0.0577d0
      LD(1,4) = -0.0431d0
      LD(1,5) = -0.0253d0
      LD(1,6) = -0.0113d0
      LD(1,7) = -0.00341d0

      LD(2,1) =  0.0407d0
      LD(2,2) =  0.0396d0
      LD(2,3) =  0.00367d0
      LD(2,4) =  0.0466d0
      LD(2,5) =  0.0623d0
      LD(2,6) =  0.0478d0
      LD(2,7) =  0.0216d0

      LD(3,1) =  0.0389d0
      LD(3,2) =  0.0193d0
      LD(3,3) =  0.0375d0
      LD(3,4) =  0.0443d0
      LD(3,5) = -0.0100d0
      LD(3,6) = -0.0612d0
      LD(3,7) = -0.0552d0

      LD(4,1) =  0.0383d0
      LD(4,2) =  1.10d-5
      LD(4,3) =  0.0486d0
      LD(4,4) =  5.40d-5
      LD(4,5) = -0.0538d0
      LD(4,6) = -9.10d-5
      LD(4,7) =  0.0738d0

      LD(5,1) =  0.0389d0
      LD(5,2) = -0.0193d0
      LD(5,3) =  0.0376d0
      LD(5,4) = -0.0443d0
      LD(5,5) = -0.0102d0
      LD(5,6) =  0.0612d0
      LD(5,7) = -0.0551d0

      LD(6,1) =  0.0407d0
      LD(6,2) = -0.0396d0
      LD(6,3) =  0.00373d0
      LD(6,4) = -0.0467d0
      LD(6,5) =  0.0623d0
      LD(6,6) = -0.0477d0
      LD(6,7) =  0.0216d0

      LD(7,1) =  0.0440d0
      LD(7,2) = -0.0627d0
      LD(7,3) = -0.0577d0
      LD(7,4) =  0.0430d0
      LD(7,5) = -0.0252d0
      LD(7,6) =  0.0113d0
      LD(7,7) = -0.00339d0

      tau = 800.d-6
      del_t = 25.d-7
      take_data = 25.d-6
      data_interval = NINT(take_data/del_t)
      num_points = NINT(tau/take_data) + 1

      call system_clock(start,count_rate)
      call specific_time(7,7,tau,config_num,LD,del_t,data_interval,
     .                   num_points)
      call system_clock(end)
      write (6,*) "time: ", real(end - start) / real(count_rate), " s"
      end 
     
      subroutine specific_time(j,k,tau,m,LD,del_t,data_interval,
     .                         num_points)
      IMPLICIT NONE
      INTEGER j,k,x,y,z,ii,iii,m,ts,total_ts,ion
      COMPLEX*16 PSI_IN(0:2**j*4**k-1),A(0:2**j*4**k-1)
     .           ,B(0:2**j*4**k-1),i
      DOUBLE PRECISION LD(j,k),pi,MODE(k,4),
     .                 SPIN(j,2),tau,P(j),PROB
      INTEGER rate,low,high,driving_mode,cnt,num_points,
     .        data_interval
      DOUBLE PRECISION del_t,t1,t2,detune(j,j)
      DOUBLE PRECISION RABI(j),w_k(k),w(j),wqbt(j),w_access(j)
      CHARACTER*50 F4
      INTEGER s1,s2,s3,s4,s5,s6,s7,m1,m2,m3,m4,m5,m6,m7,xx,cur,curp
     

 
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
c     MODE(x,y) = mode k+1-x, phonon y-1
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
      RABI(6) = 1/(9.34)
      RABI(7) = 1/(9.41)

      RABI = RABI*(pi/2)
 

      w_access(1) = (208.499094 - 3.1471)*2*pi
      w_access(2) = (208.499094 - 3.1312)*2*pi
      w_access(3) = (208.499094 - 3.1020)*2*pi
      w_access(4) = (208.499094 - 3.0679)*2*pi
      w_access(5) = (208.499094 - 3.0295)*2*pi            
      w_access(6) = (208.499094 - 2.9880)*2*pi
      w_access(7) = (208.499094 - 2.9444)*2*pi

      w(1) = (w_access(MOD(k-m,k)+1)) - detune(m+1,1)
      w(2) = (w_access(MOD(k-m+1,k)+1)) - detune(m+1,2)
      w(3) = (w_access(MOD(k-m+2,k)+1)) - detune(m+1,3)
      w(4) = (w_access(MOD(k-m+3,k)+1)) - detune(m+1,4)
      w(5) = (w_access(MOD(k-m+4,k)+1)) - detune(m+1,5)
      w(6) = (w_access(MOD(k-m+5,k)+1)) - detune(m+1,6)
      w(7) = (w_access(MOD(k-m+6,k)+1)) - detune(m+1,7)

      do ii = 0, j-1
         write(6,*) MOD(k-m+ii,k)+1
      end do

      wqbt(1) = 208.499094*2*pi
      wqbt(2) = 208.499094*2*pi
      wqbt(3) = 208.499094*2*pi
      wqbt(4) = 208.499094*2*pi
      wqbt(5) = 208.499094*2*pi
      wqbt(6) = 208.499094*2*pi
      wqbt(7) = 208.499094*2*pi

      w_k(1) = abs(w_access(1)-wqbt(1))
      w_k(2) = abs(w_access(2)-wqbt(2))
      w_k(3) = abs(w_access(3)-wqbt(3))
      w_k(4) = abs(w_access(4)-wqbt(4))
      w_k(5) = abs(w_access(5)-wqbt(5))
      w_k(6) = abs(w_access(6)-wqbt(6))
      w_k(7) = abs(w_access(7)-wqbt(7))

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
     
      ! open config files 
      write(F4, '(A,I0,A)') 'config',m,'_total.txt'
      open(unit=11, file=F4, status='unknown')

      ! find spin prob for t = 0
      do x=1,j
         CALL find_prob(PSI_IN,j,k,x,1,PROB)
         P(x) = PROB
      end do
      
      if (m .EQ. 0) then
         write(11,*) t1,(P(x), x = 1, j)
      else
         write(11,*) (P(x), x = 1, j)
      end if
      write(30,*) t1,P(1)

      ! rearrange state vector to contiguous format
      do x = 0,2**(3*j)-1
          m7 = IAND(x,3)
          m6 = IAND(ISHFT(x,-2),3)
          m5 = IAND(ISHFT(x,-4),3)
          m4 = IAND(ISHFT(x,-6),3)
          m3 = IAND(ISHFT(x,-8),3)
          m2 = IAND(ISHFT(x,-10),3)
          m1 = IAND(ISHFT(x,-12),3)

          s7 = IAND(ISHFT(x,-14),1)  
          s6 = IAND(ISHFT(x,-15),1)
          s5 = IAND(ISHFT(x,-16),1)
          s4 = IAND(ISHFT(x,-17),1)
          s3 = IAND(ISHFT(x,-18),1)
          s2 = IAND(ISHFT(x,-19),1)
          s1 = IAND(ISHFT(x,-20),1)

          xx = 0
          xx = IOR(xx,ISHFT(m1,19))
          xx = IOR(xx,ISHFT(m2,17))
          xx = IOR(xx,ISHFT(m3,15))
          xx = IOR(xx,ISHFT(m4,13))
          xx = IOR(xx,ISHFT(m5,11))
          xx = IOR(xx,ISHFT(m6,9))
          xx = IOR(xx,ISHFT(m7,7))
          
          xx = IOR(xx,ISHFT(s7,6))
          xx = IOR(xx,ISHFT(s6,5))
          xx = IOR(xx,ISHFT(s5,4))
          xx = IOR(xx,ISHFT(s4,3))
          xx = IOR(xx,ISHFT(s3,2))
          xx = IOR(xx,ISHFT(s2,1))
          xx = IOR(xx,s1)
          
          A(xx) = PSI_IN(x)
            
      end do

      cur = 1
      
      ! evolve state
C$OMP TARGET DATA MAP(TO:A(0:2097151))
C$OMP& MAP(ALLOC:B(0:2097151),PSI_IN(0:2097151))
      curp = cur
      do ts = 1,total_ts
         t2 = t1 + del_t
         CALL spin_only(j,k,m,curp,A,B,t2,t1,LD,RABI,
     .                  w,wqbt,w_k)

         CALL sigma_z(j,k,curp,A,B,t2,t1,RABI,w,wqbt,w_k)
         CALL spin_mode(j,k,m,curp,A,B,t2,t1,LD,RABI,
     .                  w,wqbt,w_k)          


         if (MOD(ts,data_interval) .EQ. 0) then

            ! back into original basis
C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TOFROM:A,B,PSI_IN)
C$OMP& PRIVATE(s1,s2,s3,s4,s5,s6,s7)
C$OMP& PRIVATE(m1,m2,m3,m4,m5,m6,m7,xx)
            do x = 0,2**(3*j)-1
               s1 = IAND(x,1)
               s2 = IAND(ISHFT(x,-1),1)
               s3 = IAND(ISHFT(x,-2),1)
               s4 = IAND(ISHFT(x,-3),1)
               s5 = IAND(ISHFT(x,-4),1)
               s6 = IAND(ISHFT(x,-5),1)
               s7 = IAND(ISHFT(x,-6),1)
                  

               m7 = IAND(ISHFT(x,-7),3)   
               m6 = IAND(ISHFT(x,-9),3)   
               m5 = IAND(ISHFT(x,-11),3)
               m4 = IAND(ISHFT(x,-13),3)
               m3 = IAND(ISHFT(x,-15),3)
               m2 = IAND(ISHFT(x,-17),3)
               m1 = IAND(ISHFT(x,-19),3)               

               xx = 0
               xx = IOR(xx,ISHFT(s1,20))
               xx = IOR(xx,ISHFT(s2,19))
               xx = IOR(xx,ISHFT(s3,18))
               xx = IOR(xx,ISHFT(s4,17))
               xx = IOR(xx,ISHFT(s5,16))
               xx = IOR(xx,ISHFT(s6,15))
               xx = IOR(xx,ISHFT(s7,14))

               xx = IOR(xx,ISHFT(m1,12))
               xx = IOR(xx,ISHFT(m2,10))
               xx = IOR(xx,ISHFT(m3,8))
               xx = IOR(xx,ISHFT(m4,6))
               xx = IOR(xx,ISHFT(m5,4))
               xx = IOR(xx,ISHFT(m6,2))
               xx = IOR(xx,m7)
               
               if (curp .EQ. 1) then  
                  PSI_IN(xx) = A(x)
               else
                  PSI_IN(xx) = B(x)
               end if

            end do
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP TARGET UPDATE FROM(PSI_IN(0:2097151))

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

            write(30,*) t2*1e6, P(1)
            cnt = cnt + 1
         end if
         t1 = t1 + del_t  
            
      end do
C$OMP END TARGET DATA

      close(11)

      return
      end

c===================================================================
c                 APPLY SPIN ONLY TERMS TO PSI
C===================================================================

      subroutine spin_only(j,k,m,cur,A,B,t2,t1,LD,RABI
     .                     ,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y,z,mode,xx,ion,m,spin_order(j,j),cur
      INTEGER s_swap,diff,s_c,s_b,s_a,m1,m2,m3,state
      DOUBLE PRECISION LD(j,k),RABI(j),wj,w_k(k),w(j),wqbt(j),t2,t1
      DOUBLE PRECISION plusx,minusy
      COMPLEX*16 i,SIG(0:1,0:1)
      COMPLEX*16 A(0:2**j*4**k-1), B(0:2**j*4**k-1),o0,o1,a0,a1
      
      i = (0.d0,1.d0)

      do z = 1,j

         ion = z 
            
         wj = wqbt(ion)-w(ion)
         plusx = RABI(ion)*(1/wj)
     .           *(sin(wj*t2)-sin(wj*t1))
         minusy = RABI(ion)*(1/wj)
     .           *(cos(wj*t1)-cos(wj*t2))
      
        SIG(0,0) = cos(sqrt(plusx**2+minusy**2))
        SIG(1,1) = cos(sqrt(plusx**2+minusy**2))
        SIG(1,0) =
     .  (sin(sqrt(plusx**2+minusy**2))/sqrt(plusx**2+minusy**2))*
     .  (-i*plusx+minusy)
        SIG(0,1) =
     .  (sin(sqrt(plusx**2+minusy**2))/sqrt(plusx**2+minusy**2))*
     .  (-i*plusx-minusy)

        if (z .EQ. j) GOTO 100
      
        if (cur .EQ. 1) then
          CALL x_carrier(z,A,B,SIG)
        else
          CALL x_carrier(z,B,A,SIG)
        end if
        cur = IEOR(cur,1)

      end do
  100 CONTINUE

      if (cur .EQ. 1) then


C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TO:SIG)
C$OMP& MAP(TOFROM:A,B)
C$OMP& PRIVATE(x,a0,a1,o0,o1)
         do x = 0,1048575
            a0 = A(2*x)
            a1 = A(2*x+1)

            o0 = SIG(0,0)*a0 + SIG(0,1)*a1
            o1 = SIG(1,0)*a0 + SIG(1,1)*a1

            A(2*x)   = o0
            A(2*x+1) = o1
         end do
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      else

C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TO:SIG)
C$OMP& MAP(TOFROM:A,B)
C$OMP& PRIVATE(x,a0,a1,o0,o1)
         do x = 0,1048575
            a0 = B(2*x)
            a1 = B(2*x+1)

            o0 = SIG(0,0)*a0 + SIG(0,1)*a1
            o1 = SIG(1,0)*a0 + SIG(1,1)*a1

            B(2*x)   = o0
            B(2*x+1) = o1
         end do
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      end if

      return      
      end

      subroutine x_carrier(z,A,B,SIG)
      INTEGER z,x,s_swap,diff,xx
      COMPLEX*16 A(0:2097151),B(0:2097151),OUTP(0:1),SIG(0:1,0:1)
      COMPLEX*16 a0,a1

C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TO:SIG)
C$OMP& MAP(TOFROM:A,B)
C$OMP& PRIVATE(x,a0,a1,OUTP,s_swap,diff,xx)
      do x=0,1048575
            a0 = A(2*x); a1 = A(2*x+1)
            OUTP(0) = SIG(0,0)*a0 + SIG(0,1)*a1
            OUTP(1) = SIG(1,0)*a0 + SIG(1,1)*a1

            ! s_last = 0,1
            s_swap = IAND(ISHFT(2*x,-z),1)
            diff = IEOR(s_swap, 0)
            xx = IEOR(2*x, IOR(ISHFT(diff,z), diff))

            B(xx) = OUTP(0)
            B(xx+2**(z)) = OUTP(1);
      end do
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      return
      end


C=================================================================
C                         APPLY SIGMA Z
C=================================================================

      subroutine sigma_z(j,k,cur,A,B,t2,t1,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,z,ion,cur
      DOUBLE PRECISION t2,t1,RABI(j),w(j),wqbt(j),w_k(k),wj
      COMPLEX*16 i,pm,hc(0:1,j)
      COMPLEX*16 A(0:2097151),B(0:2097151)

      i = (0.d0,1.d0)

      do ion=1,j
         wj = wqbt(ion) - w(ion)
         pm = -5.d-1*i*RABI(ion)*RABI(ion)*
     .      (2.d0*i*(sin(wj*t2 - wj*t1) - wj*t2 + wj*t1))
     .      / (wj**2)
         hc(0,ion) = exp(-i*pm)
         hc(1,ion) = exp(i*pm)
      end do


      if (cur .EQ. 1) then
         CALL to_spin_mode(j,k,A,B,hc)
      else
         CALL to_spin_mode(j,k,B,A,hc)
      end if
      cur = IEOR(cur,1)

      return
      end

      subroutine to_spin_mode(j,k,A,B,hc)
      INTEGER z,x,s7,s6,s5,s4,s3,s2,s1,m7,m6,m5,m4,m3,m2,m1,state,xx
      INTEGER ion,sbit,j,k,m
      COMPLEX*16 A(0:2097151),B(0:2097151),hc(0:1,j)
      COMPLEX*16 output,sz_phases
      
C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TO:hc)
C$OMP& MAP(TOFROM:A,B)
C$OMP& PRIVATE(output,sz_phases,ion,sbit,xx)
C$OMP& PRIVATE(s1,s2,s3,s4,s5,s6,s7)
C$OMP& PRIVATE(m1,m2,m3,m4,m5,m6,m7)
      do x=0,2097151
          sz_phases = (1.d0,0.d0)
  
          do ion = 1,j
            sbit = IAND(ISHFT(x,-(ion-1)),1)
            if (sbit .EQ. 0) then
               sz_phases = sz_phases * hc(0,MOD(k+ion-2,k)+1)
            else
               sz_phases = sz_phases * hc(1,MOD(k+ion-2,k)+1)
            end if
         end do
     
        
         output = A(x)*sz_phases 

         s7 = IAND(x,1)   
         s1 = IAND(ISHFT(x, -1) ,1)
         s2 = IAND(ISHFT(x, -2), 1)
         s3 = IAND(ISHFT(x, -3), 1)
         s4 = IAND(ISHFT(x, -4), 1)
         s5 = IAND(ISHFT(x, -5), 1)
         s6 = IAND(ISHFT(x, -6), 1)         

         m7 = IAND(ISHFT(x, -7), 3)   
         m6 = IAND(ISHFT(x, -9), 3)              
         m5 = IAND(ISHFT(x, -11), 3)
         m4 = IAND(ISHFT(x, -13), 3)
         m3 = IAND(ISHFT(x, -15), 3)
         m2 = IAND(ISHFT(x, -17), 3)
         m1 = IAND(ISHFT(x, -19), 3)

         xx = 0
         xx = IOR(xx, ISHFT(s7,20)) 
         xx = IOR(xx, ISHFT(m7,18))
         xx = IOR(xx, ISHFT(s6,17))
         xx = IOR(xx ,ISHFT(m6,15))
         xx = IOR(xx, ISHFT(s5,14))
         xx = IOR(xx, ISHFT(m5,12))
         xx = IOR(xx, ISHFT(s4,11))
         xx = IOR(xx, ISHFT(m4,9))
         xx = IOR(xx, ISHFT(s3,8))
         xx = IOR(xx, ISHFT(m3,6))
         xx = IOR(xx, ISHFT(s2,5))
         xx = IOR(xx, ISHFT(m2,3))
         xx = IOR(xx, ISHFT(s1,2))
         xx = IOR(xx, m1)


         B(xx) = output
      enddo
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      return
      end


c===============================================================

c===============================================================
c                 APPLY SPIN x MODE TERMS
c===============================================================

      subroutine spin_mode(j,k,m,cur,A,B,t2,t1,LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y,z,ion,mode,xx,ii,cur,m,s_c,s_b,s_a,tate
      INTEGER spin_order(j,k),swap(12,3),s(j),m1,m2,m3,m_swap1,m_swap2
      INTEGER counter,m_swap,diff,spin_order2(j,k)
      DOUBLE PRECISION wj,w(j),w_k(k),RABI(j),LD(j,k),wqbt(j),t2,t1
      DOUBLE PRECISION plusannx,plusanny,pluscrex,pluscrey,minannx,
     .                 minanny,mincrex,mincrey
      COMPLEX*16 plusann,pluscre,minann,mincre
      COMPLEX*16 i,A(0:2**j*4**k-1),B(0:2**j*4**k-1)
      COMPLEX*16 ODD(0:7,0:7),EVEN(0:7,0:7),OUTP(0:7),NOUTP(0:7)
      
      i = (0.d0,1.d0)

      do y = 1,j
         do z = 1,k
            
            ion = y
            mode = MOD(k-(y-1)+(z-1),k)+1
  
            wj = wqbt(ion)-w(ion)

            plusann = i*RABI(ion)*LD(ion,mode)*(i/(w_k(mode)+wj))
     .            *(exp((-i)*(w_k(mode)+wj)*t2)-
     .              exp((-i)*(w_k(mode)+wj)*t1))
            pluscre = i*RABI(ion)*LD(ion,mode)*((-i)/(w_k(mode)-wj))
     .           *(exp(i*((w_k(mode)-wj)*t2))-
     .             exp(i*((w_k(mode)-wj)*t1)))
            minann = (-i)*RABI(ion)*LD(ion,mode)*((-i)/(wj-w_k(mode)))*
     .            (exp(i*(wj-w_k(mode))*t2)-
     .             exp(i*(wj-w_k(mode))*t1))
            mincre = (-i)*RABI(ion)*LD(ion,mode)*((-i)/(w_k(mode)+wj))*
     .            (exp(i*(w_k(mode)+wj)*t2)-
     .             exp(i*(w_k(mode)+wj)*t1))

            plusannx = -REAL(plusann)
            plusanny = -AIMAG(plusann)
            pluscrex = -REAL(pluscre)
            pluscrey = -AIMAG(pluscre)
            minannx = -REAL(minann)
            minanny = -AIMAG(minann)
            mincrex = -REAL(mincre)
            mincrey = -AIMAG(mincre)

            EVEN = 0.d0
            ODD = 0.d0

            do x = 0,7
              EVEN(x,x) = 1.d0
              ODD(x,x) = 1.d0
            end do

            EVEN(0,0) = cos(sqrt((pluscrex/2)**2+(minanny/2)**2))
            EVEN(5,5) = cos(sqrt((pluscrex/2)**2+(minanny/2)**2))
            EVEN(5,0) = 
     .      (sin(sqrt((pluscrex/2)**2+(minanny/2)**2))/
     .      sqrt((pluscrex/2)**2+(minanny/2)**2))*
     .      (i*(pluscrex/2)-(minanny/2))
            EVEN(0,5) = 
     .      (sin(sqrt((pluscrex/2)**2+(minanny/2)**2))/
     .      sqrt((pluscrex/2)**2+(minanny/2)**2))*
     .      (i*(pluscrex/2)+(minanny/2))

            EVEN(1,1) = cos(sqrt((plusannx/2)**2+(mincrey/2)**2))
            EVEN(4,4) = cos(sqrt((plusannx/2)**2+(mincrey/2)**2))
            EVEN(4,1) =
     .      (sin(sqrt((plusannx/2)**2+(mincrey/2)**2))/
     .      sqrt((plusannx/2)**2+(mincrey/2)**2))*
     .      (i*(plusannx/2)-(mincrey/2))
            EVEN(1,4) =
     .      (sin(sqrt((plusannx/2)**2+(mincrey/2)**2))/
     .      sqrt((plusannx/2)**2+(mincrey/2)**2))*
     .      (i*(plusannx/2)+(mincrey/2))

            EVEN(2,2) = cos(sqrt((sqrt(3.d0)*pluscrex/2)**2+
     .            (sqrt(3.d0)*minanny/2)**2))
            EVEN(7,7) = cos(sqrt((sqrt(3.d0)*pluscrex/2)**2+
     .            (sqrt(3.d0)*minanny/2)**2))
            EVEN(7,2) =
     .(sin(sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))/
     .sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))*
     .      (i*(sqrt(3.d0)*pluscrex/2)-(sqrt(3.d0)*minanny/2))
            EVEN(2,7) =
     .(sin(sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))/
     .sqrt((sqrt(3.d0)*pluscrex/2)**2+(sqrt(3.d0)*minanny/2)**2))*
     .      (i*(sqrt(3.d0)*pluscrex/2)+(sqrt(3.d0)*minanny/2))

            EVEN(3,3) = cos(sqrt((sqrt(3.d0)*plusannx/2)**2+
     .            (sqrt(3.d0)*mincrey/2)**2))
            EVEN(6,6) = cos(sqrt((sqrt(3.d0)*plusannx/2)**2+
     .            (sqrt(3.d0)*mincrey/2)**2))
            EVEN(6,3) =
     .(sin(sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))/
     .sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))*
     .      (i*(sqrt(3.d0)*plusannx/2)-(sqrt(3.d0)*mincrey/2))
            EVEN(3,6) =
     .(sin(sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))/
     .sqrt((sqrt(3.d0)*plusannx/2)**2+(sqrt(3.d0)*mincrey/2)**2))*
     .      (i*(sqrt(3.d0)*plusannx/2)+(sqrt(3.d0)*mincrey/2))
    

            ODD(1,1) = cos(sqrt((sqrt(2.d0)*pluscrex)**2+
     .            (sqrt(2.d0)*minanny)**2))
            ODD(6,6) = cos(sqrt((sqrt(2.d0)*pluscrex)**2+
     .            (sqrt(2.d0)*minanny)**2))
            ODD(6,1) =
     .(sin(sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))/
     .sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))*
     .      (i*(sqrt(2.d0)*pluscrex)-(sqrt(2.d0)*minanny))
            ODD(1,6) =
     .(sin(sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))/
     .sqrt((sqrt(2.d0)*pluscrex)**2+(sqrt(2.d0)*minanny)**2))*
     .      (i*(sqrt(2.d0)*pluscrex)+(sqrt(2.d0)*minanny))


            ODD(2,2) = cos(sqrt((sqrt(2.d0)*plusannx)**2+
     .            (sqrt(2.d0)*mincrey)**2))
            ODD(5,5) = cos(sqrt((sqrt(2.d0)*plusannx)**2+
     .            (sqrt(2.d0)*mincrey)**2))
            ODD(5,2) =
     .(sin(sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))/
     .sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))*
     .      (i*(sqrt(2.d0)*plusannx)-(sqrt(2.d0)*mincrey))
            ODD(2,5) =
     .(sin(sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))/
     .sqrt((sqrt(2.d0)*plusannx)**2+(sqrt(2.d0)*mincrey)**2))*
     .      (i*(sqrt(2.d0)*plusannx)+(sqrt(2.d0)*mincrey)) 

            if (y .EQ. j .AND. z. EQ. k) GOTO 100

            if (z .EQ. k) then

               if (cur .EQ. 1) then
                  CALL x_spin_swap(y,A,B,EVEN,ODD)
               else
                  CALL x_spin_swap(y,B,A,EVEN,ODD)
               end if
               cur = IEOR(cur,1)   

            else

               if (cur .EQ. 1) then
                  CALL x_mode_swap(z,A,B,EVEN,ODD)
               else 
                  CALL x_mode_swap(z,B,A,EVEN,ODD)
               end if
               cur = IEOR(cur,1)

            end if

        end do
      end do
  100 CONTINUE

       
C     PARALLEL

      if (cur .EQ. 1) then
         CALL to_spin_only(A,B,EVEN,ODD)
      else
         CALL to_spin_only(B,A,EVEN,ODD)
      end if
      cur = IEOR(cur,1)

      return
      end

      subroutine x_mode_swap(z,A,B,EVEN,ODD)
      INTEGER x,xx,state,m_swap,diff,z
      COMPLEX*16 A(0:2097151),B(0:2097151),EVEN(0:7,0:7),
     .           ODD(0:7,0:7)
      COMPLEX*16 o0,o1,o2,o3,o4,o5,o6,o7,n0,n1,n2,n3,n4,n5,n6,n7


C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TO:EVEN,ODD)
C$OMP& MAP(TOFROM:A,B)
C$OMP& PRIVATE(x,xx,state,m_swap,diff)
C$OMP& PRIVATE(o0,o1,o2,o3,o4,o5,o6,o7)
C$OMP& PRIVATE(n0,n1,n2,n3,n4,n5,n6,n7)
      do x=0,262143
        o0 = A(8*x); o1 = A(8*x+1);
        o2 = A(8*x+2); o3 = A(8*x+3)
        o4 = A(8*x+4); o5 = A(8*x+5)
        o6 = A(8*x+6); o7 = A(8*x+7)

        n0 = EVEN(0,0)*o0 + EVEN(0,5)*o5
        n1 = EVEN(1,1)*o1 + EVEN(1,4)*o4
        n2 = EVEN(2,2)*o2 + EVEN(2,7)*o7
        n3 = EVEN(3,3)*o3 + EVEN(3,6)*o6
        n4 = EVEN(4,1)*o1 + EVEN(4,4)*o4
        n5 = EVEN(5,0)*o0 + EVEN(5,5)*o5
        n6 = EVEN(6,3)*o3 + EVEN(6,6)*o6
        n7 = EVEN(7,2)*o2 + EVEN(7,7)*o7

        o0 = n0; o1 = n1; o2 = n2; o3 = n3
        o4 = n4; o5 = n5; o6 = n6; o7 = n7

        n1 = ODD(1,1)*o1 + ODD(1,6)*o6
        n2 = ODD(2,2)*o2 + ODD(2,5)*o5
        n5 = ODD(5,2)*o2 + ODD(5,5)*o5
        n6 = ODD(6,1)*o1 + ODD(6,6)*o6

        o0 = n0; o1 = n1; o2 = n2; o3 = n3
        o4 = n4; o5 = n5; o6 = n6; o7 = n7

        n0 = EVEN(0,0)*o0 + EVEN(0,5)*o5
        n1 = EVEN(1,1)*o1 + EVEN(1,4)*o4
        n2 = EVEN(2,2)*o2 + EVEN(2,7)*o7
        n3 = EVEN(3,3)*o3 + EVEN(3,6)*o6
        n4 = EVEN(4,1)*o1 + EVEN(4,4)*o4
        n5 = EVEN(5,0)*o0 + EVEN(5,5)*o5
        n6 = EVEN(6,3)*o3 + EVEN(6,6)*o6
        n7 = EVEN(7,2)*o2 + EVEN(7,7)*o7

        CALL swap_helper1(0,x,z,n0,B)
        CALL swap_helper1(1,x,z,n1,B)
        CALL swap_helper1(2,x,z,n2,B)
        CALL swap_helper1(3,x,z,n3,B)
        CALL swap_helper1(4,x,z,n4,B)
        CALL swap_helper1(5,x,z,n5,B)
        CALL swap_helper1(6,x,z,n6,B)
        CALL swap_helper1(7,x,z,n7,B)    

      end do
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      
      return
      end

      subroutine swap_helper1(n,x,z,res,B)
      INTEGER m_swap,diff,n,z,x,xx
      COMPLEX*16 B(0:2097151),res
C$OMP DECLARE TARGET
      m_swap = IAND(ISHFT(8*x+n,-(3*z)),3)
      diff = IEOR(m_swap, IAND(n,3))
      xx = IEOR(8*x+n, IOR(ISHFT(diff, 3*z), ISHFT(diff, 0)))

      B(xx) = res
      return
      end

      subroutine x_spin_swap(y,A,B,EVEN,ODD)
      INTEGER x,xx,state,m_swap,diff,y
      COMPLEX*16 A(0:2097151),B(0:2097151),EVEN(0:7,0:7),ODD(0:7,0:7),
     .           o0,o1,o2,o3,o4,o5,o6,o7,n0,n1,n2,n3,n4,n5,n6,n7

C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TO:EVEN,ODD)
C$OMP& MAP(TOFROM:A,B)
C$OMP& PRIVATE(x,xx,state,m_swap,diff)
C$OMP& PRIVATE(o0,o1,o2,o3,o4,o5,o6,o7)
C$OMP& PRIVATE(n0,n1,n2,n3,n4,n5,n6,n7)
      do x=0,262143
        o0 = A(8*x); o1 = A(8*x+1);
        o2 = A(8*x+2); o3 = A(8*x+3)
        o4 = A(8*x+4); o5 = A(8*x+5)
        o6 = A(8*x+6); o7 = A(8*x+7)

        n0 = EVEN(0,0)*o0 + EVEN(0,5)*o5
        n1 = EVEN(1,1)*o1 + EVEN(1,4)*o4
        n2 = EVEN(2,2)*o2 + EVEN(2,7)*o7
        n3 = EVEN(3,3)*o3 + EVEN(3,6)*o6
        n4 = EVEN(4,1)*o1 + EVEN(4,4)*o4
        n5 = EVEN(5,0)*o0 + EVEN(5,5)*o5
        n6 = EVEN(6,3)*o3 + EVEN(6,6)*o6
        n7 = EVEN(7,2)*o2 + EVEN(7,7)*o7

        o0 = n0; o1 = n1; o2 = n2; o3 = n3
        o4 = n4; o5 = n5; o6 = n6; o7 = n7

        n1 = ODD(1,1)*o1 + ODD(1,6)*o6
        n2 = ODD(2,2)*o2 + ODD(2,5)*o5
        n5 = ODD(5,2)*o2 + ODD(5,5)*o5
        n6 = ODD(6,1)*o1 + ODD(6,6)*o6

        o0 = n0; o1 = n1; o2 = n2; o3 = n3
        o4 = n4; o5 = n5; o6 = n6; o7 = n7

        n0 = EVEN(0,0)*o0 + EVEN(0,5)*o5
        n1 = EVEN(1,1)*o1 + EVEN(1,4)*o4
        n2 = EVEN(2,2)*o2 + EVEN(2,7)*o7
        n3 = EVEN(3,3)*o3 + EVEN(3,6)*o6
        n4 = EVEN(4,1)*o1 + EVEN(4,4)*o4
        n5 = EVEN(5,0)*o0 + EVEN(5,5)*o5
        n6 = EVEN(6,3)*o3 + EVEN(6,6)*o6
        n7 = EVEN(7,2)*o2 + EVEN(7,7)*o7

        CALL swap_helper2(0,x,y,B,n0)
        CALL swap_helper2(1,x,y,B,n1)
        CALL swap_helper2(2,x,y,B,n2)
        CALL swap_helper2(3,x,y,B,n3)
        CALL swap_helper2(4,x,y,B,n4)
        CALL swap_helper2(5,x,y,B,n5)
        CALL swap_helper2(6,x,y,B,n6)
        CALL swap_helper2(7,x,y,B,n7)

      end do 
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      return
      end

      subroutine swap_helper2(n,x,y,B,res)
      INTEGER state,y,n,m_swap1,m_swap2,diff,xx,x
      COMPLEX*16 B(0:2097151),res
C$OMP DECLARE TARGET
      state = 8*x+n
      m_swap1 = IAND(ISHFT(state,-(3*y+2)),1)
      m_swap2 = IAND(ISHFT(state,-2),1)

      diff = IEOR(m_swap1,m_swap2)

      xx = IEOR(state,
     .     IOR(ISHFT(diff,3*y+2),ISHFT(diff,2)))
      
      B(xx) = res

      return
      end

      subroutine to_spin_only(A,B,EVEN,ODD)
      INTEGER x,xx,y,state,m1,m2,m3,m4,m5,m6,m7,s1,s2,s3,s4,s5,s6,s7
      COMPLEX*16 A(0:2097151),B(0:2097151),OUTP(0:7),NOUTP(0:7),
     .           EVEN(0:7,0:7),ODD(0:7,0:7)
C$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
C$OMP& MAP(TO:EVEN,ODD)
C$OMP& MAP(TOFROM:A,B)
C$OMP& PRIVATE(x,y,state,xx,OUTP,NOUTP)
C$OMP& PRIVATE(m1,m2,m3,m4,m5,m6,m7,s1,s2,s3,s4,s5,s6,s7)
      do x=0,262143
        OUTP(0) = A(8*x); OUTP(1) = A(8*x+1);
        OUTP(2) = A(8*x+2); OUTP(3) = A(8*x+3)
        OUTP(4) = A(8*x+4); OUTP(5) = A(8*x+5)
        OUTP(6) = A(8*x+6); OUTP(7) = A(8*x+7)
      
        NOUTP(0) = EVEN(0,0)*OUTP(0) + EVEN(0,5)*OUTP(5)
        NOUTP(1) = EVEN(1,1)*OUTP(1) + EVEN(1,4)*OUTP(4)
        NOUTP(2) = EVEN(2,2)*OUTP(2) + EVEN(2,7)*OUTP(7)
        NOUTP(3) = EVEN(3,3)*OUTP(3) + EVEN(3,6)*OUTP(6)
        NOUTP(4) = EVEN(4,1)*OUTP(1) + EVEN(4,4)*OUTP(4)
        NOUTP(5) = EVEN(5,0)*OUTP(0) + EVEN(5,5)*OUTP(5)
        NOUTP(6) = EVEN(6,3)*OUTP(3) + EVEN(6,6)*OUTP(6)
        NOUTP(7) = EVEN(7,2)*OUTP(2) + EVEN(7,7)*OUTP(7)

        OUTP = NOUTP

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

        do y = 0,7
            state = 8*x + y
            m1 = IAND(y,3)
            s7 = IAND(ISHFT(state,-2), 1)
            m2 = IAND(ISHFT(state,-3), 3)
            s1 = IAND(ISHFT(state,-5), 1)
            m3 = IAND(ISHFT(state,-6), 3)
            s2 = IAND(ISHFT(state,-8), 1)
            m4 = IAND(ISHFT(state,-9),3)
            s3 = IAND(ISHFT(state,-11),1)
            m5 = IAND(ISHFT(state,-12),3)
            s4 = IAND(ISHFT(state,-14),1)
            m6 = IAND(ISHFT(state,-15),3)
            s5 = IAND(ISHFT(state,-17),1)
            m7 = IAND(ISHFT(state,-18),3)
            s6 = IAND(ISHFT(state,-20),1)

            xx = 0
            xx = IOR(xx,ISHFT(m1,19))
            xx = IOR(xx, ISHFT(m2, 17))
            xx = IOR(xx, ISHFT(m3, 15))
            xx = IOR(xx, ISHFT(m4, 13))
            xx = IOR(xx, ISHFT(m5, 11))
            xx = IOR(xx, ISHFT(m6, 9))
            xx = IOR(xx, ISHFT(m7,7))
            xx = IOR(xx, ISHFT(s7,6))
            xx = IOR(xx, ISHFT(s6,5))
            xx = IOR(xx, ISHFT(s5, 4))
            xx = IOR(xx, ISHFT(s4, 3))
            xx = IOR(xx, ISHFT(s3, 2))
            xx = IOR(xx, ISHFT(s2, 1))
            xx = IOR(xx, s1)
            
            B(xx) = OUTP(y)
        end do 

           
      end do
C$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

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
