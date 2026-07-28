      subroutine specific_time(j,k,tau,m,LD,del_t,points,pi_times,
     .                         w_access,detune)
      IMPLICIT NONE
      INTEGER j,k,x,y,z,ii,iii,m,ts,total_ts
      COMPLEX*16 i,PSI_IN(0:2**j*4**k-1),A(0:2**j*4**k-1),
     .          B(0:2**j*4**k-1)
      DOUBLE PRECISION LD(j,k),pi,MODE(k,4),
     .                 SPIN(j,2),tau,P(7),PROB
      INTEGER rate,low,high,driving_mode,cnt,points,
     .        data_interval, mod_result
      DOUBLE PRECISION del_t,t1,t2,detune(7,7),pi_times(7)
      DOUBLE PRECISION RABI(7),w_k(7),w(7),wqbt(7),w_access(7)
      CHARACTER*50 F4
      INTEGER s1,s2,s3,s4,s5,m1,m2,m3,m4,m5,m6,m7,xx,cur,curp

 
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
      
      RABI = 0.d0
      w = 0.d0
      wqbt = 0.d0
      w_k = 0.d0        
   
      do ii = 1,j
         RABI(ii) = 1.d0/(pi_times(ii))
      end do

      RABI = RABI*(pi/2)

      w(1) = (w_access(MOD(k-m,k)+1)) - detune(m+1,1)
      w(2) = (w_access(MOD(k-m+1,k)+1)) - detune(m+1,2)
      w(3) = (w_access(MOD(k-m+2,k)+1)) - detune(m+1,3)
      w(4) = (w_access(MOD(k-m+3,k)+1)) - detune(m+1,4)
      w(5) = (w_access(MOD(k-m+4,k)+1)) - detune(m+1,5)
      w(6) = (w_access(MOD(k-m+5,k)+1)) - detune(m+1,6)
      w(7) = (w_access(MOD(k-m+6,k)+1)) - detune(m+1,7)


      do ii = 1,j
         wqbt(ii) = 208.499094d0*2.d0*pi
      end do

      w_k(1) = abs(w_access(1)-wqbt(1))
      w_k(2) = abs(w_access(2)-wqbt(2))
      w_k(3) = abs(w_access(3)-wqbt(3))
      w_k(4) = abs(w_access(4)-wqbt(4))
      w_k(5) = abs(w_access(5)-wqbt(5))
      w_k(6) = abs(w_access(6)-wqbt(6))
      w_k(7) = abs(w_access(7)-wqbt(7))

c--- set to proper units

      RABI(1) = RABI(1) * (1.d6)
      RABI(2) = RABI(2) * (1.d6)
      RABI(3) = RABI(3) * (1.d6)
      RABI(4) = RABI(4) * (1.d6)
      RABI(5) = RABI(5) * (1.d6)
      RABI(6) = RABI(6) * (1.d6)
      RABI(7) = RABI(7) * (1.d6)

      w_k(1) = w_k(1) * (1.d6)
      w_k(2) = w_k(2) * (1.d6)
      w_k(3) = w_k(3) * (1.d6)
      w_k(4) = w_k(4) * (1.d6)
      w_k(5) = w_k(5) * (1.d6)
      w_k(6) = w_k(6) * (1.d6)
      w_k(7) = w_k(7) * (1.d6)

      w(1) = w(1) * (1.d6)
      w(2) = w(2) * (1.d6)
      w(3) = w(3) * (1.d6)
      w(4) = w(4) * (1.d6)
      w(5) = w(5) * (1.d6)
      w(6) = w(6) * (1.d6)
      w(7) = w(7) * (1.d6)

      wqbt(1) = wqbt(1) * (1.d6)
      wqbt(2) = wqbt(2) * (1.d6)
      wqbt(3) = wqbt(3) * (1.d6)
      wqbt(4) = wqbt(4) * (1.d6)
      wqbt(5) = wqbt(5) * (1.d6)
      wqbt(6) = wqbt(6) * (1.d6)
      wqbt(7) = wqbt(7) * (1.d6)

c--------------------------------------
       
      total_ts = NINT(tau/del_t)
      data_interval = NINT((tau/points)/del_t)

      cnt = 0
      t1 = 0.d0
     
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
      cnt = cnt + 1
      
      ! rearrange state vector to contiguous format
      do x = 0,2**(3*j)-1
          m5 = IAND(x,3)
          m4 = IAND(ISHFT(x,-2),3)
          m3 = IAND(ISHFT(x,-4),3)
          m2 = IAND(ISHFT(x,-6),3)
          m1 = IAND(ISHFT(x,-8),3)
          
          s5 = IAND(ISHFT(x,-10),1)
          s4 = IAND(ISHFT(x,-11),1)
          s3 = IAND(ISHFT(x,-12),1)
          s2 = IAND(ISHFT(x,-13),1)
          s1 = IAND(ISHFT(x,-14),1)

          xx = 0
          xx = IOR(xx,ISHFT(m1,13))
          xx = IOR(xx,ISHFT(m2,11))
          xx = IOR(xx,ISHFT(m3,9))
          xx = IOR(xx,ISHFT(m4,7))
          xx = IOR(xx,ISHFT(m5,5))
            
          xx = IOR(xx,ISHFT(s5,4))
          xx = IOR(xx,ISHFT(s4,3))
          xx = IOR(xx,ISHFT(s3,2))
          xx = IOR(xx,ISHFT(s2,1))
          xx = IOR(xx,s1)
          
          A(xx) = PSI_IN(x)
            
      end do

      cur = 1
      
      ! evolve state
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP& PRIVATE(ts,t2,x,ii,PROB,s1,s2,s3,s4,s5,m1,m2,m3,m4,m5,xx,curp)

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
C$OMP DO SCHEDULE(STATIC)
            do x = 0,2**(3*j)-1
               s1 = IAND(x,1)
               s2 = IAND(ISHFT(x,-1),1)
               s3 = IAND(ISHFT(x,-2),1)
               s4 = IAND(ISHFT(x,-3),1)
               s5 = IAND(ISHFT(x,-4),1)

                  
               m5 = IAND(ISHFT(x,-5),3)   
               m4 = IAND(ISHFT(x,-7),3)
               m3 = IAND(ISHFT(x,-9),3)
               m2 = IAND(ISHFT(x,-11),3)
               m1 = IAND(ISHFT(x,-13),3)

               xx = 0
               xx = IOR(xx,ISHFT(s1,14))
               xx = IOR(xx,ISHFT(s2,13))
               xx = IOR(xx,ISHFT(s3,12))
               xx = IOR(xx,ISHFT(s4,11))
               xx = IOR(xx,ISHFT(s5,10))

               xx = IOR(xx,ISHFT(m1,8))
               xx = IOR(xx,ISHFT(m2,6))
               xx = IOR(xx,ISHFT(m3,4))
               xx = IOR(xx,ISHFT(m4,2))
               xx = IOR(xx,m5)

               PSI_IN(xx) = A(x)
            end do
C$OMP END DO

C$OMP SINGLE
            write (6,*) "config: ",m,"t: ",t2

            do ii = 1,j
               CALL find_prob(PSI_IN,j,k,ii,1,PROB)
               P(ii) = PROB
            end do

            if (m .EQ. 0) then
               write(11,*) t2*1.d6,(P(ii),ii=1,j)
            else
               write(11,*) (P(ii),ii=1,j)
            end if

            cnt = cnt + 1
C$OMP END SINGLE
            ! back into shifted basis
C$OMP DO SCHEDULE(STATIC)
            do x = 0,2**(3*j)-1
               m5 = IAND(x,3)
               m4 = IAND(ISHFT(x,-2),3)
               m3 = IAND(ISHFT(x,-4),3)
               m2 = IAND(ISHFT(x,-6),3)
               m1 = IAND(ISHFT(x,-8),3)

               s5 = IAND(ISHFT(x,-10),1)
               s4 = IAND(ISHFT(x,-11),1)
               s3 = IAND(ISHFT(x,-12),1)
               s2 = IAND(ISHFT(x,-13),1)
               s1 = IAND(ISHFT(x,-14),1)

               xx = 0
               xx = IOR(xx,ISHFT(m1,13))
               xx = IOR(xx,ISHFT(m2,11))
               xx = IOR(xx,ISHFT(m3,9))
               xx = IOR(xx,ISHFT(m4,7))
               xx = IOR(xx,ISHFT(m5,5))

               xx = IOR(xx,ISHFT(s5,4))
               xx = IOR(xx,ISHFT(s4,3))
               xx = IOR(xx,ISHFT(s3,2))
               xx = IOR(xx,ISHFT(s2,1))
               xx = IOR(xx,s1)

              A(xx) = PSI_IN(x)

            end do
C$OMP END DO
            
         end if
C$OMP SINGLE
         t1 = t1 + del_t  
C$OMP END SINGLE
            
      end do
C$OMP END PARALLEL

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
  
C$OMP DO PRIVATE(x,a0,a1,o0,o1) SCHEDULE(STATIC)
      do x = 0,16383
        a0 = A(2*x); a1 = A(2*x+1)

        o0 = SIG(0,0)*a0 + SIG(0,1)*a1
        o1 = SIG(1,0)*a0 + SIG(1,1)*a1
      
        A(2*x) = o0; A(2*x+1) = o1 
      end do
C$OMP END DO

      return      
      end

      subroutine x_carrier(z,A,B,SIG)
      INTEGER z,x,s_swap,diff,xx
      COMPLEX*16 A(0:32767),B(0:32767),OUTP(0:1),SIG(0:1,0:1)
      COMPLEX*16 a0,a1

C$OMP DO PRIVATE(x,a0,a1,OUTP,s_swap,diff,xx)
C$OMP& SCHEDULE(STATIC)
      do x=0,16383
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
C$OMP END DO
      return
      end


C=================================================================
C                         APPLY SIGMA Z
C=================================================================

      subroutine sigma_z(j,k,cur,A,B,t2,t1,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,z,ion,cur
      DOUBLE PRECISION t2,t1,RABI(5),w(5),wqbt(5),w_k(5),wj
      COMPLEX*16 i,pm,H(0:1,0:1),hc1,hc2
      COMPLEX*16 A(0:32767),B(0:32767)

      i = (0.d0,1.d0)

      do z = 1,j
         
         ion = j-z+1
            
      
         wj = wqbt(ion) - w(ion)

         pm = -5.d-1*i*RABI(ion)*RABI(ion)*
     .      (2.d0*i*(sin(wj*t2 - wj*t1) - wj*t2 + wj*t1))
     .      / (wj**2)


         hc1 = exp(-i*pm)
         hc2 = exp(i*pm)
            
         if (z .EQ. j) GOTO 100            
   
         if (cur .EQ. 1) then
            CALL x_sigma_z(j-z,A,B,hc1,hc2)
         else
            CALL x_sigma_z(j-z,B,A,hc1,hc2)
         end if
         cur = IEOR(cur,1)

      end do
  100 CONTINUE

         if (cur .EQ. 1) then
            CALL to_spin_mode(A,B,hc1,hc2)
         else
            CALL to_spin_mode(B,A,hc1,hc2)
         end if
         cur = IEOR(cur,1)

      return
      end


      subroutine x_sigma_z(z,A,B,hc1,hc2)
      INTEGER z,x,s_swap,diff,xx
      COMPLEX*16 A(0:32767),B(0:32767),OUTP(0:1),hc1,hc2
      COMPLEX*16 a0,a1
C$OMP DO PRIVATE(x,a0,a1,OUTP,s_swap,diff,xx)
C$OMP& SCHEDULE(STATIC)

      do x = 0,16383
         a0 = A(2*x); a1 = A(2*x+1)
         OUTP(0) = hc1*a0
         OUTP(1) = hc2*a1
      
         s_swap = IAND(ISHFT(2*x,-z),1)
         diff = IEOR(s_swap,0)
         xx = IEOR(2*x, IOR(ISHFT(diff,z), diff))

         B(xx) = OUTP(0)
         B(xx+2**(z)) = OUTP(1)

      end do
C$OMP END DO
      return
      end
      

      subroutine to_spin_mode(A,B,hc1,hc2)
      INTEGER z,x,s5,s4,s3,s2,m5,m4,m3,m2,m1,state,xx
      COMPLEX*16 A(0:32767),B(0:32767),OUTP(0:1),hc1,hc2
      COMPLEX*16 a0,a1
      
C$OMP DO SCHEDULE(STATIC)
C$OMP& PRIVATE(x,a0,a1,OUTP,state,xx)
C$OMP& PRIVATE(s2,s3,s4,s5,m1,m2,m3,m4,m5)
      do x=0,16383
         
         a0 = A(2*x)
         a1 = A(2*x+1)

         OUTP(0) = hc1*a0
         OUTP(1) = hc2*a1

         state = 2*x

         s2 = IAND(ISHFT(state, -1) ,1)
         s3 = IAND(ISHFT(state, -2), 1)
         s4 = IAND(ISHFT(state, -3), 1)
         s5 = IAND(ISHFT(state, -4), 1)

         m5 = IAND(ISHFT(state, -5), 3)
         m4 = IAND(ISHFT(state, -7), 3)
         m3 = IAND(ISHFT(state, -9), 3)
         m2 = IAND(ISHFT(state, -11), 3)
         m1 = IAND(ISHFT(state, -13), 3)

         xx = 0
         xx = IOR(xx, ISHFT(s5,14))
         xx = IOR(xx, ISHFT(m5,12))
         xx = IOR(xx, ISHFT(s4,11))
         xx = IOR(xx, ISHFT(m4,9))
         xx = IOR(xx, ISHFT(s3, 8))
         xx = IOR(xx, ISHFT(m3, 6))
         xx = IOR(xx, ISHFT(s2,5))
         xx = IOR(xx, ISHFT(m2,3))
c         xx = IOR(xx, ISHFT(s1,2))
         xx = IOR(xx, m1)


      B(xx) = OUTP(0)
      B(xx + 4) = OUTP(1)
      enddo
C$OMP END DO
      return
      end


c===============================================================

c===============================================================
c                 APPLY SPIN x MODE TERMS
c===============================================================

      subroutine spin_mode(j,k,m,cur,A,B,t2,t1,LD,RABI,w,wqbt,w_k)
      IMPLICIT NONE
      INTEGER j,k,x,y,z,ion,mode,xx,ii,cur,m,s_c,s_b,s_a,state
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
      COMPLEX*16 A(0:32767),B(0:32767),EVEN(0:7,0:7),
     .           ODD(0:7,0:7)
      COMPLEX*16 o0,o1,o2,o3,o4,o5,o6,o7,n0,n1,n2,n3,n4,n5,n6,n7
C$OMP DO PRIVATE(x)
C$OMP& SCHEDULE(STATIC)
      do x=0,4095
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
C$OMP END DO
      
      return
      end

      subroutine swap_helper1(n,x,z,res,B)
      INTEGER m_swap,diff,n,z,x,xx
      COMPLEX*16 B(0:32767),res

      m_swap = IAND(ISHFT(8*x+n,-(3*z)),3)
      diff = IEOR(m_swap, IAND(n,3))
      xx = IEOR(8*x+n, IOR(ISHFT(diff, 3*z), ISHFT(diff, 0)))

      B(xx) = res
      return
      end

      subroutine x_spin_swap(y,A,B,EVEN,ODD)
      INTEGER x,xx,state,m_swap,diff,y
      COMPLEX*16 A(0:32767),B(0:32767),EVEN(0:7,0:7),ODD(0:7,0:7),
     .           o0,o1,o2,o3,o4,o5,o6,o7,n0,n1,n2,n3,n4,n5,n6,n7

C$OMP DO PRIVATE(x)
C$OMP& SCHEDULE(STATIC)
      do x=0,4095
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
C$OMP END DO

      return
      end

      subroutine swap_helper2(n,x,y,B,res)
      INTEGER state,y,n,m_swap1,m_swap2,diff,xx,x
      COMPLEX*16 B(0:32767),res

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
      INTEGER x,xx,y,state,m1,m2,m3,m4,m5,s1,s2,s3,s4,s5
      COMPLEX*16 A(0:32767),B(0:32767),OUTP(0:7),NOUTP(0:7),
     .           EVEN(0:7,0:7),ODD(0:7,0:7)
C$OMP DO SCHEDULE(STATIC)
C$OMP& PRIVATE(x,y,state,xx,OUTP,NOUTP)
C$OMP& PRIVATE(m1,m2,m3,m4,m5,s1,s2,s3,s4,s5)
      do x=0,4095
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
            s5 = IAND(ISHFT(state,-2), 1)
            m2 = IAND(ISHFT(state,-3), 3)
            s1 = IAND(ISHFT(state,-5), 1)
            m3 = IAND(ISHFT(state,-6), 3)
            s2 = IAND(ISHFT(state,-8), 1)
            m4 = IAND(ISHFT(state,-9),3)
            s3 = IAND(ISHFT(state,-11),1)
            m5 = IAND(ISHFT(state,-12),3)
            s4 = IAND(ISHFT(state,-14),1)

            xx = 0
            xx = IOR(xx,ISHFT(m1,13))
            xx = IOR(xx, ISHFT(m2, 11))
            xx = IOR(xx, ISHFT(m3, 9))
            xx = IOR(xx, ISHFT(m4, 7))
            xx = IOR(xx, ISHFT(m5, 5))
            xx = IOR(xx, ISHFT(s5, 4))
            xx = IOR(xx, ISHFT(s4, 3))
            xx = IOR(xx, ISHFT(s3, 2))
            xx = IOR(xx, ISHFT(s2, 1))
            xx = IOR(xx, s1)
            
            B(xx) = OUTP(y)
        end do 

           
      end do
C$OMP END DO

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

      subroutine run_model(j,k,LD,config,points,
     .                     tau,w_access,pi_times,detune)

      implicit none
      INTEGER j,k,config,points
      double precision LD(j,k)
      double precision tau,del_t,pi_times(7)
      double precision detune(7,7),w_access(7)

      del_t = 125.d-7

      call specific_time(j,k,tau,config-1,LD,del_t,
     .                   points,pi_times,w_access,detune)

      return
      end
