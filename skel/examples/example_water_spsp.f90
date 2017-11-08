program SimintTest
use SimintFortran
use iso_c_binding
implicit none
type(c_simint_shell), target :: s_shell(4),p_shell(1)
type(c_simint_multi_shellpair), target :: s2p_msh, s1p_msh

double precision :: alpha(1), coef(1)
double precision :: integrals(1000)
double precision, allocatable :: work(:)

double precision, pointer :: p1(:), p2(:)
integer :: i, ncomputed
integer :: worksize

call simint_init()


worksize = simint_eri_worksize(0, 3)
allocate(work(worksize))



call simint_initialize_shell(s_shell(1))
call simint_initialize_shell(s_shell(2))
call simint_initialize_shell(s_shell(3))
call simint_initialize_shell(s_shell(4))
call simint_initialize_shell(p_shell(1))
call simint_initialize_multi_shellpair(s2p_msh)
call simint_initialize_multi_shellpair(s1p_msh)
! oxygen s
alpha(1)=1d0
coef(1) =1d0
call simint_create_shell(1, 0, 0.0d0, 0.0d0, -0.02d0, &
                          alpha(1), coef(1), s_shell(1)) 
! oxygen sp
alpha(1)=0.4d0
coef(1) =1d0
call simint_create_shell(1, 0, 0.0d0, 0.0d0, -0.02d0, &
                         alpha(1), coef(1), s_shell(2)) 
call simint_create_shell(1, 1, 0.0d0, 0.0d0, -0.02d0, &
                        alpha(1), coef(1), p_shell(1)) 
!hydrogen s
alpha(1)=0.5d0
coef(1) =1d0
call simint_create_shell(1, 0, -0.74d0, 0.0d0, -0.76d0, &
                         alpha(1), coef(1), s_shell(3)) 
call simint_create_shell(1, 0,  0.74d0, 0.0d0, -0.76d0, &
                         alpha(1), coef(1), s_shell(4)) 

call C_F_POINTER(s_shell(1)%alpha, p1, shape=[s_shell(1)%nprim])
call C_F_POINTER(s_shell(1)%coef, p2, shape=[s_shell(1)%nprim])
write(*,*) "s_shell1  info"
do i = 1, s_shell(1)%nprim
  write(*,*) p1(i), p2(i)
end do
call C_F_POINTER(s_shell(2)%alpha, p1, shape=[s_shell(2)%nprim])
call C_F_POINTER(s_shell(2)%coef, p2, shape=[s_shell(2)%nprim])
write(*,*) "s_shell1  info"
do i = 1, s_shell(2)%nprim
  write(*,*) p1(i), p2(i)
end do
call C_F_POINTER(p_shell(1)%alpha, p1, shape=[p_shell(1)%nprim])
call C_F_POINTER(p_shell(1)%coef, p2, shape=[p_shell(1)%nprim])
write(*,*) "p_shell  info"
do i = 1, s_shell(1)%nprim
  write(*,*) p1(i), p2(i)
end do

call simint_normalize_shells(4, s_shell)
call simint_normalize_shells(1, p_shell(1))

call simint_create_multi_shellpair(1, s_shell(4), 1, p_shell(1), s2p_msh, 0)
call simint_create_multi_shellpair(1, s_shell(3), 1, p_shell(1), s1p_msh, 0)

call C_F_POINTER(s1p_msh%alpha, p1, shape=[s1p_msh%nprim])
call C_F_POINTER(s1p_msh%prefac, p2, shape=[s1p_msh%nprim])
write(*,*) "s2p Shell Pair info"
do i = 1, s1p_msh%nprim
  write(*,*) p1(i), p2(i)
end do

call C_F_POINTER(s2p_msh%alpha, p1, shape=[s2p_msh%nprim])
call C_F_POINTER(s2p_msh%prefac, p2, shape=[s2p_msh%nprim])
write(*,*) "s2p Shell Pair info"
do i = 1, s2p_msh%nprim
  write(*,*) p1(i), p2(i)
end do


ncomputed = simint_compute_eri(s2p_msh, s1p_msh, 0.0d0, work, integrals)
write(6,*) ' ncomputed ',ncomputed
do i=1,9
   write(*,'(i4,f20.6)') i,integrals(i)
enddo

deallocate(work)
call simint_free_shell(s_shell(1))
call simint_free_shell(s_shell(2))
call simint_free_shell(s_shell(3))
call simint_free_shell(s_shell(4))
call simint_free_shell(p_shell(1))
call simint_free_multi_shellpair(s2p_msh)
call simint_free_multi_shellpair(s1p_msh)

call simint_finalize()

end program
