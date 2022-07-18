c----------------------------------------------------------------------
c
        module rotatingBandForce
          real*8, allocatable  :: rotBandForce(:,:)
        end module
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine malloc_rotBandForce
c
        use rotatingBandForce
c
        allocate( rotBandForce (numRotBands, 3) )
c
        rotBandForce  = zero
        write(*,*) 'inside malloc_rotBandForce'
        return

        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine init_rotBandForce
c
        use rotatingBandForce
        use number_def_m
c
        rotBandForce  = zero
c
        return
        end
c
c----------------------------------------------------------------------
c




        subroutine commu_rotBandForce
c
        use rotatingBandForce
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        real*8, dimension(3) :: Forin, Forout
cRane   
        write(*,*) 'inside commu_rotBandForce'
c
        do i = 1, numRotBands
          if (iter .eq. nitr) then
            Forin  = (/ rotBandForce(i,1), rotBandForce(i,2), rotBandForce(i,3) /)
            if (numpe > 1) then
              call MPI_ALLREDUCE ( Forin,  Forout, 3,
     &                             MPI_DOUBLE_PRECISION,
     &                             MPI_SUM, MPI_COMM_WORLD, ierr)
              rotBandForce(i,1:3) = Forout(1:3)
            endif
          endif
        enddo
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        subroutine release_rotBandForce
c
        use rotatingBandForce
c
        include "common.h"
cRane   
        write(*,*) 'inside release_rotBandForce'
c
        if (allocated(rotBandForce))
     &    deallocate( rotBandForce )
c
        return
        end
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
