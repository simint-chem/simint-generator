module SimintFortran
    use iso_c_binding
    implicit none


    type, bind(C) :: c_simint_shell
      integer(C_INT) :: am
      integer(C_INT) :: nprim
      real(C_DOUBLE) :: x
      real(C_DOUBLE) :: y
      real(C_DOUBLE) :: z
      type(C_PTR) :: alpha
      type(C_PTR) :: coef
      integer(C_SIZE_T) :: memsize
      type(C_PTR) :: ptr
    end type c_simint_shell

    type, bind(C) :: c_simint_multi_shellpair
      integer(C_INT) :: am1
      integer(C_INT) :: am2
      integer(C_INT) :: nprim
      integer(C_INT) :: nshell12
      integer(C_INT) :: nshell12_clip
      type(C_PTR) :: nprim12
      type(C_PTR) :: AB_x
      type(C_PTR) :: AB_y
      type(C_PTR) :: AB_z
      type(C_PTR) :: x
      type(C_PTR) :: y
      type(C_PTR) :: z
      type(C_PTR) :: PA_x
      type(C_PTR) :: PA_y
      type(C_PTR) :: PA_z
      type(C_PTR) :: PB_x
      type(C_PTR) :: PB_y
      type(C_PTR) :: PB_z
      type(C_PTR) :: alpha

#if SIMINT_OSTEI_MAXDER > 0
        type(C_PTR) :: alpha2
        type(C_PTR) :: beta2
#endif

      type(C_PTR) :: prefac
      type(C_PTR) :: screen
      real(C_DOUBLE) :: screen_max
      integer(C_SIZE_T) :: memsize
      type(C_PTR) :: ptr
    end type

  interface

    subroutine simint_init() bind(C, name="simint_init")
    end subroutine

    subroutine simint_finalize() bind(C, name="simint_finalize")
    end subroutine

    subroutine c_simint_initialize_shell(G) bind(C, name="simint_initialize_shell")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: G
    end subroutine

    subroutine c_simint_initialize_shells(n, G) bind(C, name="simint_initialize_shells")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: G
      integer(C_INT), intent(in), value :: n
    end subroutine

    subroutine c_simint_allocate_shell(nprim, G) bind(C, name="simint_allocate_shell")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: nprim
      type(C_PTR), intent(in), value :: G
    end subroutine

    subroutine c_simint_free_shell(G) bind(C, name="simint_free_shell")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: G
    end subroutine

    subroutine c_simint_free_shells(n, G) bind(C, name="simint_free_shells")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: G
      integer(C_INT), intent(in), value :: n
    end subroutine

    subroutine c_simint_normalize_shells(n, G) bind(C, name="simint_normalize_shells")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: G
      integer(C_INT), intent(in), value :: n
    end subroutine

    subroutine c_simint_copy_shell(src, dest) bind(C, name="simint_copy_shell")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: src, dest
    end subroutine

    subroutine c_simint_create_shell(nprim, am, &
               x, y, z, alpha, coef, G) bind(C, name="simint_create_shell")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: nprim, am
      real(C_DOUBLE), intent(in), value :: x, y, z
      type(C_PTR), intent(in), value :: alpha, coef, G
    end subroutine

    subroutine c_simint_create_zero_shell(G) bind(C, name="simint_create_zero_shell")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: G
    end subroutine

    subroutine c_simint_initialize_multi_shellpair(P) &
               bind(C, name="simint_initialize_multi_shellpair")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: P
    end subroutine

    subroutine c_simint_initialize_multi_shellpairs(n, P) bind(C, name="simint_initialize_multi_shellpairs")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: P
      integer(C_INT), intent(in), value :: n
    end subroutine

    subroutine c_simint_free_multi_shellpair(P) &
               bind(C, name="simint_free_multi_shellpair")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: P
    end subroutine

    subroutine c_simint_free_multi_shellpairs(n, P) bind(C, name="simint_free_multi_shellpairs")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: P
      integer(C_INT), intent(in), value :: n
    end subroutine

    subroutine c_simint_create_multi_shellpair( &
               na, A, nb, B, P, screen_method) &
               bind(C, name="simint_create_multi_shellpair")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: A, B
      type(C_PTR), intent(in), value :: P
      integer(C_INT), intent(in), value :: na, nb, screen_method
    end subroutine

    subroutine c_simint_create_multi_shellpair2( &
               npair, AB, P, screen_method) &
               bind(C, name="simint_create_multi_shellpair2")

      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: AB, P
      integer(C_INT), intent(in), value :: npair, screen_method
    end subroutine

    function c_simint_compute_eri(P, Q, screen_tol, work, integrals) &
             result(res) bind(C, name="simint_compute_eri")
      use iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: P, Q, work, integrals
      real(C_DOUBLE), intent(in), value :: screen_tol
      integer(C_INT) :: res
    end function

    function c_simint_eri_worksize(derorder, maxam) &
             result(res) bind(C, name="simint_eri_worksize")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: derorder, maxam
      integer(C_SIZE_T) :: res
    end function

    function c_simint_eri_workmem(derorder, maxam) &
             result(res) bind(C, name="simint_eri_workmem")
      use iso_c_binding
      implicit none
      integer(C_INT), intent(in), value :: derorder, maxam
      integer(C_SIZE_T) :: res
    end function

  end interface

  contains

    subroutine simint_initialize_shell(G)
      implicit none
      type(c_simint_shell), intent(inout), target :: G

      call c_simint_initialize_shell(C_LOC(G))
    end subroutine

    subroutine simint_initialize_shells(n, G)
      implicit none
      integer, intent(in) :: n
      type(c_simint_shell), intent(inout), target :: G(n)

      call c_simint_initialize_shells(INT(n, C_INT), C_LOC(G))
    end subroutine

    subroutine simint_allocate_shell(nprim, G)
      implicit none
      integer :: nprim
      type(c_simint_shell), intent(inout), target :: G

      call c_simint_allocate_shell(INT(nprim, C_INT), C_LOC(G))
    end subroutine

    subroutine simint_free_shell(G)
      implicit none
      type(c_simint_shell), intent(inout), target :: G

      call c_simint_free_shell(C_LOC(G))
    end subroutine

    subroutine simint_free_shells(n, G)
      implicit none
      integer, intent(in) :: n
      type(c_simint_shell), intent(inout), target :: G(n)

      call c_simint_free_shells(INT(n, C_INT), C_LOC(G))
    end subroutine

    subroutine simint_normalize_shells(n, G)
      implicit none
      integer, intent(in) :: n
      type(c_simint_shell), intent(inout), target :: G(n)

      call c_simint_normalize_shells(INT(n, C_INT), C_LOC(G))
    end subroutine

    subroutine simint_copy_shell(src, dest)
      implicit none
      type(c_simint_shell), intent(in), target :: src
      type(c_simint_shell), intent(out), target :: dest

      call c_simint_copy_shell(C_LOC(src), C_LOC(dest))
    end subroutine

    subroutine simint_create_shell(nprim, am, &
               x, y, z, alpha, coef, G)
      implicit none
      integer, intent(in) :: nprim, am
      double precision, intent(in) :: x, y, z
      double precision, intent(in), target :: alpha(nprim), coef(nprim)
      type(c_simint_shell), intent(inout), target :: G

      call c_simint_create_shell(INT(nprim, C_INT), &
                                 INT(am, C_INT), &
                                 REAL(x, C_DOUBLE), &
                                 REAL(y, C_DOUBLE), &
                                 REAL(z, C_DOUBLE), &
                                 C_LOC(alpha), &
                                 C_LOC(coef), &
                                 C_LOC(G))

    end subroutine

    subroutine simine_create_zero_shell(G)
      implicit none
      type(c_simint_shell), intent(inout), target :: G

      call c_simint_create_zero_shell(C_LOC(G))
    end subroutine

    subroutine simint_initialize_multi_shellpair(P)
      implicit none
      type(c_simint_multi_shellpair), intent(inout), target :: P

      call c_simint_initialize_multi_shellpair(C_LOC(P))
    end subroutine

    subroutine simint_initialize_multi_shellpairs(n, P)
      implicit none
      integer, intent(in) :: n
      type(c_simint_multi_shellpair), intent(inout), target :: P(n)

      call c_simint_initialize_multi_shellpairs(INT(n, C_INT), C_LOC(P))
    end subroutine

    subroutine simint_free_multi_shellpair(P)
      implicit none
      type(c_simint_multi_shellpair), intent(inout), target :: P

      call c_simint_free_multi_shellpair(C_LOC(P))
    end subroutine

    subroutine simint_free_multi_shellpairs(n, P)
      implicit none
      integer, intent(in) :: n
      type(c_simint_multi_shellpair), intent(inout), target :: P(n)

      call c_simint_free_multi_shellpairs(INT(n, C_INT), C_LOC(P))
    end subroutine

    subroutine simint_create_multi_shellpair( &
               na, A, nb, B, P, screen_method)
      implicit none
      integer, intent(in) :: na, nb, screen_method
      type(c_simint_multi_shellpair), intent(inout), target :: P
      type(c_simint_shell), intent(in), target :: A(na), B(nb)

      call c_simint_create_multi_shellpair(INT(na, C_INT), &
                                           C_LOC(A), &
                                           INT(nb, C_INT), &
                                           C_LOC(B), &
                                           C_LOC(P), &
                                           INT(screen_method, C_INT))


    end subroutine

    subroutine simint_create_multi_shellpair2( &
               npair, AB, P, screen_method)
      implicit none
      integer, intent(in) :: npair, screen_method
      type(c_simint_shell), intent(in), target :: AB(npair)
      type(c_simint_multi_shellpair), intent(inout), target :: P

      call c_simint_create_multi_shellpair2(INT(npair, C_INT), &
                                            C_LOC(AB), C_LOC(P), &
                                            INT(screen_method, C_INT))
    end subroutine

    function simint_compute_eri(P, Q, screen_tol, work, integrals) &
             result(res)
      implicit none

      type(c_simint_multi_shellpair), intent(in), target :: P, Q
      double precision, intent(inout), target :: work(*), integrals(*)
      real(C_DOUBLE), intent(in) :: screen_tol
      integer :: res
      integer(C_INT) :: res2

      res2 = c_simint_compute_eri(C_LOC(P), C_LOC(Q), &
                                  REAL(screen_tol, C_DOUBLE), &
                                  C_LOC(work), C_LOC(integrals))
      res = INT(res2)
    end function

    function simint_eri_worksize(derorder, maxam) result(res)
      implicit none
      integer, intent(in) :: derorder, maxam
      integer :: res
      integer(C_INT) :: res2

      res2 = c_simint_eri_worksize(INT(derorder, C_INT), INT(maxam, C_INT))
      res = INT(res2)
    end function

    function simint_eri_workmem(derorder, maxam) result(res)
      implicit none
      integer, intent(in) :: derorder, maxam
      integer :: res
      integer(C_INT) :: res2

      res2 = c_simint_eri_workmem(INT(derorder, C_INT), INT(maxam, C_INT))
      res = INT(res2)
    end function

end module
