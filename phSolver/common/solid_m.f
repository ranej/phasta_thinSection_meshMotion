      module solid_data_m
        use pointer_data
        use outpar_m
        implicit none
c
        type (r3d), dimension(MAXBLK2) ::  b !for solid,left Cauchy_green tensor,added
        type (r3d), dimension(MAXBLK2) ::  b_dot!time derivative of b,added
        type (r3d), dimension(MAXBLK2) ::  b_af! b at time step n+af,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b !left Cauchy_green tensor on the boundary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_dot!time derivative of b on the boudary,added
        type (r3d), dimension(MAXBLK2) ::  bdy_b_af! b at time step n+af on the boudary,added
c
        integer, pointer :: is_solid(:)
        integer, parameter :: b_size = 6
c
        type solid_t
          logical :: is_active, restart
          integer :: nel,nelb
        end type solid_t
c
        integer :: b_array_size
        real*8, dimension(:), pointer :: temp_b, temp_b_dot, temp_b_af  ! temp arrays to read b, b_dot and b_af
        real*8, dimension(:), pointer :: bdy_temp_b, bdy_temp_b_dot, bdy_temp_b_af  ! on the boundary...
c
        type(solid_t) :: solid_p
c
      end module solid_data_m
c
      module solid_m
c
        use solid_data_m
        implicit none
c
        contains
c
        subroutine init_block(b,bdot,baf,tempb,tempbdot,tempbaf)
          use intpt_m
          use propar_m
          implicit none
          real*8, dimension(:,:,:), pointer, intent(out) :: b,bdot,baf
          real*8, dimension(:), pointer, intent(in) :: tempb,tempbdot,tempbaf
          integer :: ib,ifirst,ilast,incr
          do intp = 1,ngauss
            do ib = 1,b_size
              ifirst = ib + (intp-1)*ngauss
              incr   = ngauss*b_size
              ilast  = ifirst + (npro-1)*incr
              b(1:npro,intp,ib)    = tempb(ifirst:ilast:incr)
              bdot(1:npro,intp,ib) = tempbdot(ifirst:ilast:incr)
              baf(1:npro,intp,ib)  = tempbaf(ifirst:ilast:incr)
            enddo
          enddo
        end subroutine init_block
c
        subroutine dump_block(tempb,tempbdot,tempbaf,b,bdot,baf)
          use intpt_m
          use propar_m
          implicit none
          real*8, dimension(:), pointer, intent(out) :: tempb,tempbdot,tempbaf
          real*8, dimension(:,:,:), pointer, intent(in) :: b,bdot,baf
          integer :: ib,ifirst,ilast,incr
          do intp = 1,ngauss
            do ib = 1,b_size
              ifirst = ib + (intp-1)*ngauss
              incr   = ngauss*b_size
              ilast  = ifirst + (npro-1)*incr
              tempb(ifirst:ilast:incr)    = b(1:npro,intp,ib)
              tempbdot(ifirst:ilast:incr) = bdot(1:npro,intp,ib)
              tempbaf(ifirst:ilast:incr)  = baf(1:npro,intp,ib)
            enddo
          enddo
        end subroutine dump_block
c
        subroutine malloc_solid
          use matdat_def_m
          use number_def_m
          use elmpar_m
          use blkdat_m
          use intpt_m
          implicit none
          integer :: mattype, iblk, npro
c
c.... allocate space for solid arrays 
c
          solid_p%nel = 0
          b_array_size = 0  ! keep this size for later write_field...
c
          blocks_loop: do iblk = 1, nelblk
c
            mattype = lcblk(iblk_mattype,iblk)
            lcsyst = lcblk(3,iblk)
            ngauss = nint(lcsyst)
            npro = lcblk(1,iblk+1) - lcblk(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
c
              allocate (b(iblk)%p(npro,ngauss,b_size))
              allocate (b_dot(iblk)%p(npro,ngauss,b_size))
              allocate (b_af(iblk)%p(npro,ngauss,b_size))
c
              solid_p%nel = solid_p%nel + npro
              b_array_size = b_array_size + npro*ngauss*b_size
c
              if (solid_p%restart) then
c
                call init_block(b(iblk)%p,b_dot(iblk)%p,b_af(iblk)%p,temp_b,temp_b_dot,temp_b_af)
c
              else
c
                b(iblk)%p(:,:,:)= one
                b_af(iblk)%p(:,:,:) = one
c
              endif
            else
              cycle   
            endif
c
          enddo blocks_loop
c
          if (solid_p%restart) then
            deallocate(temp_b)
            deallocate(temp_b_dot)
            deallocate(temp_b_af)
          endif
c
          solid_p%nelb = 0
c
          boundary_blocks_loop: do iblk = 1, nelblb
c
            mattype = lcblkb(iblk_mattype,iblk)
            lcsyst = lcblkb(3,iblk)
            ngauss = nintb(lcsyst)
            npro = lcblkb(1,iblk+1) - lcblkb(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
c
              allocate (bdy_b(iblk)%p(npro,ngaussb,b_size))
              allocate (bdy_b_dot(iblk)%p(npro,ngaussb,b_size))
              allocate (bdy_b_af(iblk)%p(npro,ngaussb,b_size))
c
              solid_p%nelb = solid_p%nelb + npro
c
              if (solid_p%restart) then
c
                call init_block(bdy_b(iblk)%p,bdy_b_dot(iblk)%p,bdy_b_af(iblk)%p,bdy_temp_b,bdy_temp_b_dot,bdy_temp_b_af)
c
              else
c
                bdy_b(iblk)%p(:,:,:)= one
                bdy_b_af(iblk)%p(:,:,:) = one
c
              endif
            else
              cycle
            endif
c..
            if (solid_p%restart) then
              deallocate(bdy_temp_b)
              deallocate(bdy_temp_b_dot)
              deallocate(bdy_temp_b_af)
            endif
c
          enddo boundary_blocks_loop
c
c
        end subroutine malloc_solid
c
        subroutine free_solid
c
        end subroutine free_solid
c
        subroutine read_field(fieldtag,b)
          use iso_c_binding 
          use phio
          implicit none
          character(len=*), intent(in) :: fieldtag
          real*8, pointer, intent(out) :: b(:)
          character(len=1024) :: dataInt, dataDbl
          integer, target :: intfromfile(50) ! integers read from headers
          intfromfile=0
c          call phio_readheader(fhandle,
c     &     c_char_'trim(fieldtag)' // char(0),
c     &     c_loc(intfromfile), 1, dataInt, iotype)
c          if (intfromfile(1) > 0) then
c            allocate(b(intfromfile(1)))
c            call phio_readdatablock(fhandle,
c     &       c_char_'trim(fieldtag)' // char(0),
c     &       c_loc(b), intfromfile(1), dataDbl, iotype)
c          endif
        end subroutine read_field
c
        subroutine read_restart_solid
          use iso_c_binding 
          use phio
          implicit none
          call read_field('solid b',temp_b)
          call read_field('solid b_dot',temp_b_dot)
          call read_field('solid b_af',temp_b_af)
          if (associated(temp_b)) then
            solid_p%restart   = .true.
          endif
        end subroutine read_restart_solid
c
        subroutine write_restart_solid
          use matdat_def_m
          use number_def_m
          use elmpar_m
          use blkdat_m
          use intpt_m
          use workfc_m
          use timdat_m
          implicit none
          integer :: mattype, iblk, npro, b_len
c
          allocate(temp_b(b_array_size))
          allocate(temp_b_dot(b_array_size))
          allocate(temp_b_af(b_array_size))
c
          b_len = 0
          blocks_loop: do iblk = 1, nelblk
c
            mattype = lcblk(iblk_mattype,iblk)
            lcsyst = lcblk(3,iblk)
            ngauss = nint(lcsyst)
            npro = lcblk(1,iblk+1) - lcblk(1,iblk)
c
            if (mat_eos(mattype,1).eq.ieos_solid_1)then
              call dump_block(temp_b,temp_b_dot,temp_b_af,b(iblk)%p,b_dot(iblk)%p,b_af(iblk)%p)
            endif
c
          enddo blocks_loop
c
          call write_field(myrank,'a'//char(0),'solid b'//char(0),7,
     &     temp_b,'d'//char(0),b_array_size,1,lstep)
c
          deallocate(temp_b)
          deallocate(temp_b_dot)
          deallocate(temp_b_af)
c
        end subroutine write_restart_solid
c
      end module solid_m
