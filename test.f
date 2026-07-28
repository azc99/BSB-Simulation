      program test
      IMPLICIT NONE
      INTEGER k,m,z,j
      
      j = 6; k = 6; m = 0;

      do z = 1,j
         write(6,*) MOD(k-(m+1)+(z-1),k)+1
      end do
      
      do z = 1,j
         write(6,*) MOD(k-(m+2)+(z-1),k)+1
      end do

      return
      end
