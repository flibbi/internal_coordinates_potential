  PROGRAM forcefield
!
!
!
!
  IMPLICIT none

  INTEGER :: n_atx, n_at, n_atsc, n_atscx, n_bon, n_bonx,    &
             n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo, nrwsx,    &
             n_l1x, n_l2x, n_l3x, nqptx, nbndx
  PARAMETER (n_atx = 2, n_atscx = 600, n_bonx = 800, n_fc_nnx = 4,   &
             n_fc_box = 9, nrwsx=200, n_l1x=16, n_l2x=16, n_l3x=1,   &
             nqptx=1000, nbndx=6)

  INTEGER :: n_l1, n_l2, n_l3

  INTEGER :: i0_d, i1_d, i2_d, i3_d, i1_at, i2_at, i1_l1, i1_l2,  i2_l1, i2_l2,   &
             il, il1, shift1, shift2, i_nn, ia,ib,ic, i_bon, i1_uc, i1_sc, i2_sc,       &
             m1_uc, m2_uc, i2_uc, m1_sc, m2_sc, iq, nqpt
  INTEGER :: ind_il_nn(n_atscx,3), ind_sh_nn(n_atscx,3,2), ind_il_bon(n_bonx,6),   &
             ind_sh_bon(n_bonx,6,2)

  REAL(kind=8) :: a_lat(3,3), a_latd(3,3),r_uc(3,n_atx), r_sc(3,n_atscx), vect_1(3), vect_2(3),   &
                  vect_3(3), vect_4(3), vect_bond(3,6),  vect_nn(3,4), pos(3,n_atscx), &
                  pos_pp(3,n_atscx),pos_mp(3,n_atscx),pos_pm(3,n_atscx),pos_mm(3,n_atscx), &
                  xxq(3,nqptx), eV_to_J, A_to_m, uma_to_kg, THz_to_cmm1
  REAL(kind=8) :: k_nn(n_fc_nnx), k_bo(n_fc_box), e_nn, e_bo, e_bond, e_tot, ene, ene_pp, ene_pm, &
                  ene_mp, ene_mm, summ, ened(3), strain
  REAL(kind=8) :: p_sc(3), dist, dist_nn, prod12(3), delta, fc_mat(n_atx,3,n_atscx,3),    &
                  frc(n_l1x,n_l2x,n_l3x,3,3,n_atx,n_atx), a_modulus, mass

  REAL(kind=8) :: atws(3,3),rws2(0:3,nrwsx),xq(3), bg(3,3), epsil(3,3), zeu(3,3,n_atx),    &
                  freq2(3*n_atx), freq_ifc(nqptx,3*n_atx), const, omega
  INTEGER :: nrws2, nsc, itau(n_atx)
  COMPLEX(kind=8) :: dyn(3,3,n_atx,n_atx), dyn_(3,3,n_atx,n_atx), dync(3*n_atx,3*n_atx),   &
                     z(3*n_atx,3*n_atx)

  CHARACTER*80 file_cell, file_fc, file_qp 

  namelist /input/ file_cell, file_fc, file_qp, n_l1, n_l2, n_l3
  file_cell = ' '
  file_fc = ' '
  file_qp = ' '
  n_l1 = 4
  n_l2 = 4
  n_l3 = 1

  read(5,input)

!
! Few constants
!

  eV_to_J = 1.60218d-19
  A_to_m = 1.d-10
  uma_to_kg = 1.66054d-27
  THz_to_cmm1 = 33.35641
  a_modulus = 2.467 
  mass = 12.d0

  const = dsqrt(eV_to_J / (uma_to_kg * A_to_m**2) ) / (2*3.14159)* 1d-12 * THz_to_cmm1 


  WRITE(*,*) ' '
  WRITE(*,*) '------------------------------------------------------------------------------'
  WRITE(*,*) '|     RUNNING THE ICP CODE                                                   |'
  WRITE(*,*) '------------------------------------------------------------------------------'  
!
! Read the information about the unit cell
!
  OPEN(unit=50,file=file_cell,status='old',form='formatted')
  DO i1_d = 1,3
    READ(50,*)(a_lat(i2_d, i1_d), i2_d = 1,3)
  ENDDO

  READ(50,*) n_at
  DO i1_at = 1, n_at
    READ(50,*) (r_uc(i2_d, i1_at), i2_d = 1,3)
  ENDDO

!
! Read the force constant 
!
  n_fc_nn = 4
  n_fc_bo = 9
  OPEN(unit=55,file=file_fc,status='old',form='formatted')
  DO i1_d = 1,n_fc_nn
    READ(55,*)k_nn(i1_d)
  ENDDO
  DO i1_d = 1,n_fc_bo
    READ(55,*)k_bo(i1_d)
  ENDDO

!
! Read the qpoints 
!
  OPEN(unit=60,file=file_qp,status='old',form='formatted')
  READ(60,*)nqpt
  DO iq = 1, nqpt
    READ(60,*)(xxq(ia,iq), ia=1,3)
  ENDDO



!
! Look for the 3 nearest neighbor (nn) in the supercell
!

  n_atsc = n_at*n_l1*n_l2  ! atoms in the supercell
  WRITE(*,*) ' '
  WRITE(*,*) '1) BUILDING GEOMETRY '
  WRITE(*,*) 'Number of atoms:', n_atsc

! Atoms in the supercell
  OPEN(unit=300,file='supercell.dat',form = 'formatted')
  WRITE(300,*) 'Identity      x              y              z'
  DO i1_at = 1, n_at
    DO i1_l1 = 0, n_l1 - 1
      DO i1_l2 = 0, n_l2 - 1
         il = i1_at + n_at*i1_l1 + n_at*n_l1*i1_l2
         DO ia = 1,3
           r_sc(ia,il) = r_uc(ia,i1_at) + a_lat(ia,1)*DFLOAT(i1_l1)  &
                                        + a_lat(ia,2)*DFLOAT(i1_l2)
         ENDDO !ia
         WRITE(300,'(i4,2f16.8)') il, r_sc(1,il), r_sc(2,il)
      ENDDO ! i1_l1
    ENDDO ! i1_l2
  ENDDO ! i1_at
  CLOSE(unit=300)

! Find the nearest neighbors
  dist_nn = (1.d0 + 0.05d0)/DSQRT(3.d0)
  ind_il_nn(:,:) = 0
  ind_sh_nn(:,:,:) = 0

  DO i1_at = 1, n_atsc
    i_nn = 0
    DO i2_at = 1, n_at
      DO i1_l1 = -1, n_l1
        DO i1_l2 = -1, n_l2
          il = i2_at + n_at*i1_l1 + n_at*n_l1*i1_l2
          IF(il.EQ.i1_at) GOTO 99
          dist = 0.d0
          DO ia = 1,3
            p_sc(ia) = r_uc(ia,i2_at) + a_lat(ia,1)*DFLOAT(i1_l1)   &
                                      + a_lat(ia,2)*DFLOAT(i1_l2)
            dist = dist + ( r_sc(ia,i1_at) - p_sc(ia) )**2
          ENDDO !ia
          dist = DSQRT(dist)
          IF(dist.LT.dist_nn) THEN
            i_nn = i_nn + 1
            i2_l1 = i1_l1
            i2_l2 = i1_l2
            shift1 = 0
            shift2 = 0
            IF(i1_l1.eq.-1) THEN
              i2_l1 = i1_l1 + n_l1
              shift1 = n_l1
              GOTO 77
            ENDIF
            IF(i1_l1.eq.n_l1) THEN
              i2_l1 = i1_l1 - n_l1
              shift1 = -n_l1
            ENDIF
77          CONTINUE
            IF(i1_l2.eq.-1) THEN
              i2_l2 = i1_l2 + n_l2
              shift2 = n_l2
              goto 88
            ENDIF
            IF(i1_l2.eq.n_l2) THEN
              i2_l2 = i1_l2 - n_l2
              shift2 = -n_l2
            ENDIF
88          CONTINUE
            il1 = i2_at + n_at*i2_l1 + n_at*n_l1*i2_l2
            ind_il_nn(i1_at,i_nn) = il1
            ind_sh_nn(i1_at,i_nn,1) = shift1
            ind_sh_nn(i1_at,i_nn,2) = shift2 
          ENDIF          
99        CONTINUE
        ENDDO ! i1_l1
      ENDDO ! i1_l2
    ENDDO !i2_at  
  ENDDO !i1_at

  DO i_bon = 1, n_atsc
    DO i1_at = 1, 3
      DO ia = 1, 3
        vect_bond(ia,i1_at) = r_sc(ia, ind_il_nn(i_bon,i1_at))                  &
                              - a_lat(ia,1)*DFLOAT( ind_sh_nn(i_bon,i1_at,1) )  &
                              - a_lat(ia,2)*DFLOAT( ind_sh_nn(i_bon,i1_at,2) )  
      ENDDO ! ia
    ENDDO !i1_at
  ENDDO !i_bon


  n_bon = 0
  ind_il_bon(:,:) = 0
  ind_sh_bon(:,:,:) = 0
  DO i1_at = 1, n_atsc
    DO i_nn = 1, 3

      DO i_bon = 1, n_bon
        IF(i1_at.EQ.ind_il_bon(i_bon,1) .OR. i1_at.EQ.ind_il_bon(i_bon,2)) THEN
          IF(ind_il_nn(i1_at,i_nn).EQ.ind_il_bon(i_bon,1) .OR. ind_il_nn(i1_at,i_nn).EQ.ind_il_bon(i_bon,2)) THEN
            GOTO 109
          ENDIF
        ENDIF
      ENDDO !i_bon

        n_bon = n_bon + 1
        ind_il_bon(n_bon,1) = i1_at
        ind_il_bon(n_bon,2) = ind_il_nn(i1_at,i_nn)

        ind_sh_bon(n_bon,2,1) = ind_sh_nn(i1_at,i_nn,1)
        ind_sh_bon(n_bon,2,2) = ind_sh_nn(i1_at,i_nn,2)        

        i1_d = mod(i_nn,3) + 1
        i2_d = mod(i1_d,3) + 1


        DO ia = 1, 3
          vect_1(ia) = r_sc(ia, ind_il_nn(i1_at,i_nn))                      &
                       - a_lat(ia,1)*DFLOAT( ind_sh_nn(i1_at,i_nn,1) )     &
                       - a_lat(ia,2)*DFLOAT( ind_sh_nn(i1_at,i_nn,2) ) - r_sc(ia,i1_at)
          vect_2(ia) = r_sc(ia, ind_il_nn(i1_at,i1_d))                      &
                       - a_lat(ia,1)*DFLOAT( ind_sh_nn(i1_at,i1_d,1) )     &
                       - a_lat(ia,2)*DFLOAT( ind_sh_nn(i1_at,i1_d,2) ) - r_sc(ia,i1_at)

        ENDDO !ia
        CALL prod_vect(vect_1, vect_2, prod12)
        IF(prod12(3).GT.0.d0) THEN
          ind_il_bon(n_bon,3) = ind_il_nn(i1_at,i1_d)
          ind_il_bon(n_bon,4) = ind_il_nn(i1_at,i2_d)

          ind_sh_bon(n_bon,3,1) = ind_sh_nn(i1_at,i1_d,1)
          ind_sh_bon(n_bon,3,2) = ind_sh_nn(i1_at,i1_d,2)
          ind_sh_bon(n_bon,4,1) = ind_sh_nn(i1_at,i2_d,1)
          ind_sh_bon(n_bon,4,2) = ind_sh_nn(i1_at,i2_d,2)

        ELSE
          ind_il_bon(n_bon,4) = ind_il_nn(i1_at,i1_d)
          ind_il_bon(n_bon,3) = ind_il_nn(i1_at,i2_d)

          ind_sh_bon(n_bon,4,1) = ind_sh_nn(i1_at,i1_d,1)
          ind_sh_bon(n_bon,4,2) = ind_sh_nn(i1_at,i1_d,2)
          ind_sh_bon(n_bon,3,1) = ind_sh_nn(i1_at,i2_d,1)
          ind_sh_bon(n_bon,3,2) = ind_sh_nn(i1_at,i2_d,2)

        ENDIF

        i0_d = ind_il_nn(i1_at,i_nn)
        DO ia = 1, 3
          IF( ind_il_nn(i0_d,ia) .EQ. i1_at) i1_d = ia           
        ENDDO !ia
        i2_d = mod(i1_d,3) + 1
        i3_d = mod(i2_d,3) + 1

        DO ia = 1, 3
          vect_1(ia) = r_sc(ia,ind_il_nn(i0_d,i1_d))                      &
                       - a_lat(ia,1)*DFLOAT( ind_sh_nn(i0_d,i1_d,1) )     &
                       - a_lat(ia,2)*DFLOAT( ind_sh_nn(i0_d,i1_d,2) ) - r_sc(ia,i0_d)
          vect_2(ia) = r_sc(ia,ind_il_nn(i0_d,i2_d))                      &
                       - a_lat(ia,1)*DFLOAT( ind_sh_nn(i0_d,i2_d,1) )     &
                       - a_lat(ia,2)*DFLOAT( ind_sh_nn(i0_d,i2_d,2) ) - r_sc(ia,i0_d)


        ENDDO !ia


        CALL prod_vect(vect_1, vect_2, prod12)       
        IF(prod12(3).GT.0.d0) THEN
          ind_il_bon(n_bon,5) = ind_il_nn(i0_d,i2_d)
          ind_il_bon(n_bon,6) = ind_il_nn(i0_d,i3_d)

          ind_sh_bon(n_bon,5,1) = ind_sh_nn(i0_d,i2_d,1) + ind_sh_bon(n_bon,2,1)
          ind_sh_bon(n_bon,5,2) = ind_sh_nn(i0_d,i2_d,2) + ind_sh_bon(n_bon,2,2)
          ind_sh_bon(n_bon,6,1) = ind_sh_nn(i0_d,i3_d,1) + ind_sh_bon(n_bon,2,1) 
          ind_sh_bon(n_bon,6,2) = ind_sh_nn(i0_d,i3_d,2) + ind_sh_bon(n_bon,2,2)
        ELSE
          ind_il_bon(n_bon,6) = ind_il_nn(i0_d,i2_d)
          ind_il_bon(n_bon,5) = ind_il_nn(i0_d,i3_d)

          ind_sh_bon(n_bon,6,1) = ind_sh_nn(i0_d,i2_d,1) + ind_sh_bon(n_bon,2,1)
          ind_sh_bon(n_bon,6,2) = ind_sh_nn(i0_d,i2_d,2) + ind_sh_bon(n_bon,2,2)
          ind_sh_bon(n_bon,5,1) = ind_sh_nn(i0_d,i3_d,1) + ind_sh_bon(n_bon,2,1) 
          ind_sh_bon(n_bon,5,2) = ind_sh_nn(i0_d,i3_d,2) + ind_sh_bon(n_bon,2,2)
        ENDIF

109   CONTINUE

    ENDDO ! i_nn
  ENDDO !i1_at

  WRITE(*,*)'Number of bonds:', n_bon

  DO i_bon = 1, n_bon
    DO i1_at = 1, 6
      DO ia = 1,3
        vect_bond(ia,i1_at) = r_sc(ia,ind_il_bon(i_bon,i1_at))                  &
                              - a_lat(ia,1)*DFLOAT( ind_sh_bon(i_bon,i1_at,1) )  &
                              - a_lat(ia,2)*DFLOAT( ind_sh_bon(i_bon,i1_at,2) )  
      ENDDO ! ia
    ENDDO !i1_at

  ENDDO !i_bon



!
! Second derivative calculation (IFC)
!
  WRITE(*,*)' '
  WRITE(*,*)'2) CALCULATING IFCs'
  OPEN(unit=1000,file='force-constants.dat',form = 'formatted')
  WRITE(1000,*) 'l1   l2   coord_supercell   coord_unitcell   ind_supercell   ind_unitcell   IFC'

  delta = 0.0001
  DO i1_uc = 1, n_at
    DO m1_uc = 1, 3
      DO i1_at = 1, n_at
        DO i1_l1 = 0, n_l1 - 1
          DO i1_l2 = 0, n_l2 - 1
            i1_sc = i1_at + n_at*i1_l1 + n_at*n_l1*i1_l2
            DO m1_sc = 1, 3
              IF( i1_uc.NE.i1_sc .OR. m1_uc.NE.m1_sc ) THEN
                DO i2_sc = 1, n_atsc
                  DO m2_sc = 1, 3            
                    pos_pp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
                    pos_pm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
                    pos_mp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
                    pos_mm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
                    IF( i2_sc.EQ.i1_uc .AND. m2_sc.EQ.m1_uc ) THEN
                      pos_pp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) + delta
                      pos_pm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) + delta
                      pos_mp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) - delta
                      pos_mm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) - delta          
                    ENDIF
                    IF( i2_sc.EQ.i1_sc .AND. m2_sc.EQ.m1_sc ) THEN
                      pos_pp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) + delta
                      pos_pm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) - delta
                      pos_mp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) + delta
                      pos_mm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) - delta
                    ENDIF
                  ENDDO !m2_sc
                ENDDO !i2_sc
          
                CALL energy(ene_pp,pos_pp,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                            n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,             &
                            ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)          
                CALL energy(ene_pm,pos_pm,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                            n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,             &
                            ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)
                CALL energy(ene_mp,pos_mp,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                            n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,             &
                            ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)
                CALL energy(ene_mm,pos_mm,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                            n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,             &
                            ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)
                frc(i1_l1+1,i1_l2+1,1,m1_sc,m1_uc,i1_at,i1_uc) = ( (ene_pp - ene_pm) - (ene_mp - ene_mm) )/   &
                                                                (4.d0*delta*delta) / a_modulus**2 / mass
                WRITE(1000,'(6i5,f16.8)') i1_l1+1,i1_l2+1, m1_sc, m1_uc, i1_at, i1_uc, &
                                           frc(i1_l1+1,i1_l2+1,1,m1_sc,m1_uc,i1_at,i1_uc)
              ENDIF
            ENDDO !m1_sc
          ENDDO !i1_l2
        ENDDO !i1_l1
      ENDDO !i1_at
    ENDDO !m1_uc
  ENDDO !i1_uc

  DO i1_uc = 1, n_at
    DO m1_uc = 1, 3
      DO i2_sc = 1, n_atsc
        DO m2_sc = 1, 3            
          pos_pp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
          pos_pm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
          pos_mp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
          pos_mm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
          IF( i2_sc.EQ.i1_uc .AND. m2_sc.EQ.m1_uc ) THEN
            pos_pp(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) + 2.d0*delta
            pos_pm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc)
            pos_mm(m2_sc,i2_sc) = r_sc(m2_sc,i2_sc) - 2.d0*delta         
          ENDIF
        ENDDO !m2_sc
      ENDDO !i2_sc
      CALL energy(ene_pp,pos_pp,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                  n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,             &
                  ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)          
      CALL energy(ene_pm,pos_pm,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                  n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,             &
                  ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)
      CALL energy(ene_mm,pos_mm,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                  n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,             &
                  ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)

      frc(1,1,1,m1_uc,m1_uc,i1_uc,i1_uc) = ( ene_pp - 2.d0*ene_pm + ene_mm ) / &
                                           (4.d0*delta*delta)/a_modulus**2 /mass
      WRITE(1000,'(6i5,f16.8)') i1_l1+1,i1_l2+1, m1_sc, m1_uc, i1_at, i1_uc, &
                               frc(i1_l1+1,i1_l2+1,1,m1_sc,m1_uc,i1_at,i1_uc)
    ENDDO !m1_uc
  ENDDO !i1_uc

  WRITE(*,*)' '
  WRITE(*,*)'3) TESTING ACOUSTIC SUM RULE SIMPLE'
  DO i1_uc = 1, n_at
    DO m1_uc = 1, 3
      DO m1_sc = 1, 3 
        summ = 0.d0
        DO i1_at = 1, n_at
          DO i1_l1 = 1, n_l1
            DO i1_l2 = 1, n_l2 
              summ = summ + frc(i1_l1,i1_l2,1,m1_sc,m1_uc,i1_at,i1_uc)*mass
            ENDDO !i1_l2
          ENDDO !i1_l1
        ENDDO !i1_at
        WRITE(*,'(3i5,f16.8)') i1_uc, m1_uc, m1_sc, summ
      ENDDO !m1_sc
    ENDDO !m1_uc
  ENDDO !i1_uc
  CLOSE(unit=1000)


!
! Calculation of the dynamical matrix
!
  WRITE(*,*) ' '    
  WRITE(*,*) '4) CALCULATING PHONONS'
  CALL recips(a_lat(1,1),a_lat(1,2),a_lat(1,3),   &
              bg(1,1),bg(1,2),bg(1,3))

  itau(1) = 1
  itau(2) = 2

  omega = a_lat(1,1)*a_lat(2,2)*a_lat(3,3)
  epsil(:,:)=0.d0
  zeu(:,:,:)=0.d0

      call DCOPY(9,a_lat,1,atws,1)
      call DSCAL(3,dfloat(n_l1),atws(1,1),1)
      call DSCAL(3,dfloat(n_l2),atws(1,2),1)
      call DSCAL(3,dfloat(n_l3),atws(1,3),1)
      call wsinit(rws2,nrwsx,nrws2,atws) ! 

  OPEN(unit=7777,file='all-modes.dat',form = 'formatted')
  OPEN(unit=800,file='frequencies.dat',form = 'formatted')
  OPEN(unit=801,file='frequencies.gnu',form = 'formatted')

  DO iq = 1, nqpt
    dyn(:,:,:,:) = (0.0d0,0.0d0)
    nsc = 1
    call setupmat (xxq(1,iq),dyn,n_at,n_atx,a_lat,bg,r_uc,itau,nsc, &
                     1.d0,dyn_, n_at,n_atx,a_lat,bg,r_uc,omega, &
                     epsil,zeu,frc,n_l1,n_l2,n_l3,n_l1x,n_l2x,n_l3x, &
                     .false.,rws2,nrws2)
    !
    ! Dyagonalize
    !
    DO i1_uc = 1, n_at
      DO m1_uc = 1,3
        i1_d = m1_uc + 3*(i1_uc - 1)
        DO i2_uc = 1, n_at
          DO m2_uc = 1,3
            i2_d = m2_uc + 3*(i2_uc - 1)
            dync(i1_d,i2_d) = dyn(m1_uc,m2_uc,i1_uc,i2_uc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO 

    CALL cdiagh2(6,dync,6,freq2,z)
    WRITE(7777,'(A5,2f16.8)') 'q = ', xxq(1,iq), xxq(2,iq)      
    WRITE(7777,'(A78)') ' *****************************************************************************'
    DO ia=1,6
      IF(freq2(ia).lt.0.d0)THEN
        freq_ifc(iq,ia) = -1.d0*DSQRT(-1.d0*freq2(ia))*const
      ELSE
        freq_ifc(iq,ia) = DSQRT(freq2(ia))*const
      ENDIF
      WRITE(7777,'(A9,i2,A4,f16.6,A8)') ' freq(', ia, ')=', freq_ifc(iq,ia), '[cm-1]'
      DO i1_uc = 1, n_at
        WRITE(7777,'(A3,f14.6,f10.6,f14.6,f10.6,f14.6,f10.6,A3)') &
        '(', (DREAL(z(ib+(i1_uc-1)*3,ia)), DIMAG(z(ib+(i1_uc-1)*3,ia)),ib=1,3), ')'
      ENDDO
    ENDDO  

    WRITE(800,'(2f8.4,6f10.4)')xxq(1,iq),xxq(2,iq),(freq_ifc(iq,ia), ia=1,6)
  ENDDO !iq    

  DO ia=1,6
    DO iq = 1, nqpt
      WRITE(801,*)iq, freq_ifc(iq,ia)
    ENDDO !iq
  ENDDO   
  
  CLOSE(unit=7777)
  CLOSE(unit=800)
  CLOSE(unit=801)

  WRITE(*,*) '  '
  WRITE(*,*) 'FINISHED'
  WRITE(*,*) '  '
  STOP
  END

!######################################################################
! Subroutine energy
!######################################################################
  SUBROUTINE energy(ene,pos,k_nn,k_bo,a_lat,n_atsc,n_atscx,n_bon,n_bonx,   &
                    n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo,            &
                    ind_il_nn,ind_sh_nn,ind_il_bon,ind_sh_bon)
  IMPLICIT NONE
  INTEGER :: n_atsc, n_atscx, n_bon, n_bonx, n_fc_nnx, n_fc_box, n_fc_nn, n_fc_bo
  INTEGER :: ind_il_nn(n_atscx,3), ind_sh_nn(n_atscx,3,2), ind_il_bon(n_bonx,6),   &
             ind_sh_bon(n_bonx,6,2)
  REAL(kind=8) :: pos(3,n_atscx), ene, k_nn(n_fc_nnx), k_bo(n_fc_box), a_lat(3,3)

  INTEGER :: i1_at, ia, i_nn, i_bon
  REAL(kind=8) :: e_nn, e_bo, vect_nn(3,4), vect_bond(3,6) 

!
! ENERGY CALCULATION
!
  ene = 0.d0



! nn-body term
  DO i1_at = 1, n_atsc
    DO ia = 1, 3
      vect_nn(ia,1) = pos(ia,i1_at) 
    ENDDO !ia
    DO i_nn = 1, 3
      DO ia = 1, 3
        vect_nn(ia,i_nn+1) = pos(ia,ind_il_nn(i1_at,i_nn))     &
                             - a_lat(ia,1)*DFLOAT( ind_sh_nn(i1_at,i_nn,1) )  &
                             - a_lat(ia,2)*DFLOAT( ind_sh_nn(i1_at,i_nn,2) ) 
      ENDDO !ia
    ENDDO !i_nn

    CALL pot_nn(k_nn, n_fc_nnx, n_fc_nn, vect_nn, e_nn)
    ene = ene + e_nn

  ENDDO !i1_at

! bo term


  DO i_bon = 1, n_bon
    DO i1_at = 1, 6
      DO ia = 1,3
        vect_bond(ia,i1_at) = pos(ia,ind_il_bon(i_bon,i1_at))                  &
                              - a_lat(ia,1)*DFLOAT( ind_sh_bon(i_bon,i1_at,1) )  &
                              - a_lat(ia,2)*DFLOAT( ind_sh_bon(i_bon,i1_at,2) )  
      ENDDO ! ia
    ENDDO !i1_at

   CALL pot_bo(k_bo, n_fc_box, n_fc_bo, vect_bond, e_bo)
   ene = ene + e_bo

  ENDDO !i_bon


  RETURN
  END SUBROUTINE energy


!######################################################################
! Subroutine bo potential
!######################################################################
  SUBROUTINE pot_bo(k_bo, n_fc_box, n_fc_bo, vect_bond, e_bo)
  IMPLICIT NONE
  INTEGER :: n_fc_box, n_fc_bo
  REAL(kind=8) :: vect_bond(3,6), k_bo(n_fc_box), e_bo, d(4), ang(4),    &
                  prod, v(3,4), chi(6),      &
                  v_chi(3,4), mv_chi(4), v12(3), v21(3)
  REAL(kind=8) :: pi, d0, ang0, chi0, ec_1, ec_2, ec_3, ec_4,   &
                  ec_5, ec_6, ec_7, ec_8, ec_9, att

  INTEGER :: ia, ib, ic

  pi = DATAN(1.d0)*4.d0
  d0 = 1.d0/DSQRT(3.d0)
  ang0 = 2.d0/3.d0*pi
  chi0 = 0.d0

  DO ia = 1, 2
    DO ib = 1, 3
      v(ib,ia) = vect_bond(ib,ia+2) - vect_bond(ib,1)
    ENDDO !ib
  ENDDO !ia
  DO ia = 3, 4
    DO ib = 1, 3
      v(ib,ia) = vect_bond(ib,ia+2) - vect_bond(ib,2)
    ENDDO !ib
  ENDDO !ia


  DO ib = 1, 3
    v12(ib) = vect_bond(ib,2) - vect_bond(ib,1)
    v21(ib) = -v12(ib)
  ENDDO !ib

! lenght
  DO ia = 1, 4
    d(ia) = 0.d0
    DO ib = 1, 3
      d(ia) = d(ia) + v(ib,ia)**2
    ENDDO !ib
    d(ia) = DSQRT( d(ia) )
  ENDDO !ia

! angle

  DO ia = 1, 2
    prod = 0.d0
    DO ib = 1, 3
      prod = prod + v(ib,ia)*v12(ib)
    ENDDO !ic
    ang(ia) = DACOS( prod/(d(ia)*DSQRT( v12(1)**2+v12(2)**2+v12(3)**2 ) ) )
  ENDDO !ia



  DO ia = 3, 4
    prod = 0.d0
    DO ib = 1, 3
      prod = prod + v(ib,ia)*v21(ib)
    ENDDO !ic
    ang(ia) = DACOS( prod/(d(ia)*DSQRT( v21(1)**2+v21(2)**2+v21(3)**2 ) ) ) 
  ENDDO !ia


! dihedral angle


  CALL prod_vect(v12, v(1,1), v_chi(1,1))
  CALL prod_vect(v(1,2), v12, v_chi(1,2))
  CALL prod_vect(v21, v(1,3), v_chi(1,3))
  CALL prod_vect(v(1,4), v21, v_chi(1,4))



  DO ia = 1, 4
    prod = 0.d0
    DO ib = 1, 3 
      prod = prod + v_chi(ib,ia)*v_chi(ib,ia)
    ENDDO !ib
    mv_chi(ia) = DSQRT(prod) 
  ENDDO !ia

  chi(:) = 0.d0

! m=2
  prod = 0.d0
  DO ic = 1, 3
    prod = prod + v_chi(ic,1)*v_chi(ic,2)
  ENDDO !ic
  att = prod/(mv_chi(1)*mv_chi(2))
  IF(att.gt.1.d0) THEN
    !WRITE(9999,*)'Warning'
    !WRITE(9999,'(f20.16)')att
    IF(att-1.d0 .LT. 1.d-14) THEN
      att = 1.d0
    ELSE
      WRITE(9999,*)'STOPPING'
      STOP
    ENDIF
  ENDIF
  chi(1) = DACOS( att )   

  prod = 0.d0
  DO ic = 1, 3
    prod = prod + v_chi(ic,3)*v_chi(ic,4)
  ENDDO !ic
  att = prod/(mv_chi(3)*mv_chi(4))
  IF(att.gt.1.d0) THEN
    !WRITE(9999,*)'Warning'
    !WRITE(9999,'(f20.16)')att
    IF(att-1.d0 .LT. 1.d-14) THEN
      att = 1.d0
    ELSE
      WRITE(9999,*)'STOPPING'
      STOP
    ENDIF
  ENDIF
  chi(2) = DACOS( att )



! m=3
  prod = 0.d0
  DO ic = 1, 3
    prod = prod + v_chi(ic,2)*v_chi(ic,3)
  ENDDO !ic
  att = prod/(mv_chi(2)*mv_chi(3))
  IF(att.gt.1.d0) THEN
    !WRITE(9999,*)'Warning'
    !WRITE(9999,'(f20.16)')att
    IF(att-1.d0 .LT. 1.d-14) THEN
      att = 1.d0
    ELSE
      WRITE(9999,*)'STOPPING'
      STOP
    ENDIF
  ENDIF
  chi(3) = DACOS( att )

   

  prod = 0.d0
  DO ic = 1, 3
    prod = prod + v_chi(ic,1)*v_chi(ic,4)
  ENDDO !ic
  att = prod/(mv_chi(1)*mv_chi(4))
  IF(att.gt.1.d0) THEN
    !WRITE(9999,*)'Warning'
    !WRITE(9999,'(f20.16)')att
    IF(att-1.d0 .LT. 1.d-14) THEN
      att = 1.d0
    ELSE
      WRITE(9999,*)'STOPPING'
      STOP
    ENDIF
  ENDIF
  chi(4) = DACOS( att )



! m=4
  prod = 0.d0
  DO ic = 1, 3
    prod = prod + v_chi(ic,1)*v_chi(ic,3)
  ENDDO !ic
  att = prod/(mv_chi(1)*mv_chi(3))
  IF(att.gt.1.d0) THEN
    !WRITE(9999,*)'Warning'
    !WRITE(9999,'(f20.16)')att
    IF(att-1.d0 .LT. 1.d-14) THEN
      att = 1.d0
    ELSE
      WRITE(9999,*)'STOPPING'
      STOP
    ENDIF
  ENDIF
  chi(5) = DACOS( att )



  prod = 0.d0
  DO ic = 1, 3
    prod = prod + v_chi(ic,2)*v_chi(ic,4)
  ENDDO !ic 
  att = prod/(mv_chi(2)*mv_chi(4))
  IF(att.gt.1.d0) THEN
    !WRITE(9999,*)'Warning'
    !WRITE(9999,'(f20.16)')att
    IF(att-1.d0 .LT. 1.d-14) THEN
      att = 1.d0
    ELSE
      WRITE(9999,*)'STOPPING'
      STOP
    ENDIF
  ENDIF
  chi(6) = DACOS( att )



! Delta
  DO ia = 1, 4
    d(ia) = d(ia) - d0
    ang(ia) = ang(ia) - ang0
  ENDDO !ia
  DO ia = 1, 6
    chi(ia) = chi(ia) - chi0
  ENDDO !ia

  ec_1 = k_bo(1)*( d(2)*d(3) + d(1)*d(4) )
  ec_2 = k_bo(2)*( ang(2)*ang(3) + ang(1)*ang(4) )
  ec_3 = k_bo(3)*( d(2)*ang(3) + d(1)*ang(4)  +  &
                   d(3)*ang(2) + d(4)*ang(1)  )

  ec_4 = k_bo(4)*( d(2)*d(4) + d(1)*d(3) )
  ec_5 = k_bo(5)*( ang(2)*ang(4) + ang(1)*ang(3) )
  ec_6 = k_bo(6)*( d(2)*ang(4) + d(1)*ang(3)  +  &
                   d(4)*ang(2) + d(3)*ang(1)  )

  ec_7 = k_bo(7)*(chi(1)**2 + chi(2)**2)
  ec_8 = k_bo(8)*(chi(3)**2 + chi(4)**2)
  ec_9 = k_bo(9)*(chi(5)**2 + chi(6)**2)


  e_bo = ec_1 + ec_2 + ec_3 +   &
         ec_4 + ec_5 + ec_6 +   &
         ec_7 + ec_8 + ec_9


  RETURN
  END SUBROUTINE pot_bo




!######################################################################
! Subroutine nn potential
!######################################################################
  SUBROUTINE pot_nn(k_nn, n_fc_nnx, n_fc_nn, vect_nn, e_nn)
  IMPLICIT NONE
  INTEGER :: n_fc_nnx, n_fc_nn
  REAL(kind=8) :: vect_nn(3,4), k_nn(n_fc_nnx), e_nn, d(3), ang(3), prod, v(3,3)
  REAL(kind=8) :: pi, d0, ang0, ec_1, ec_2, ec_3, ec_4, ec_corr

  INTEGER :: ia, ib, ic

  pi = DATAN(1.d0)*4.d0
  d0 = 1.d0/DSQRT(3.d0)
  ang0 = 2.d0/3.d0*pi

  DO ia = 1, 3
    DO ib = 1, 3
      v(ib,ia) = vect_nn(ib,ia+1) - vect_nn(ib,1)
    ENDDO !ib
  ENDDO !ia

! lenght
  DO ia = 1, 3
    d(ia) = 0.d0
    DO ib = 1, 3
      d(ia) = d(ia) + v(ib,ia)**2
    ENDDO !ib
    d(ia) = DSQRT( d(ia) )
  ENDDO !ia

! angle
  DO ia = 1, 3
    ib = mod(ia,3) + 1
    prod = 0.d0
    DO ic = 1, 3
      prod = prod + v(ic,ib)*v(ic,ia)
    ENDDO !ic
    ang(ia) = DACOS( prod/(d(ia)*d(ib)) )
  ENDDO !ia

  DO ia = 1, 3
    d(ia) = d(ia) -d0
    ang(ia) = ang(ia) - ang0
  ENDDO !ia

  ec_1 = 2.d0*k_nn(1)*( d(1)**2 + d(2)**2 + d(3)**2 )
  ec_2 = k_nn(2)*( ang(1)**2 + ang(2)**2 + ang(3)**2 )
  ec_3 = k_nn(3)*( d(1)*d(2) + d(1)*d(3) + d(3)*d(2) )
  ec_4 = k_nn(4)*( ang(1)*(d(1)+d(2)) + ang(2)*(d(2)+d(3)) + ang(3)*(d(3)+d(1)) )
  ec_corr = 2 * (-117.60) *  ( d(1)**3 + d(2)**3 + d(3)**3 )

  e_nn = ec_1 + ec_2 + ec_3 + ec_4 + ec_corr


  RETURN
  END SUBROUTINE pot_nn




!######################################################################
! Subroutine nn potential
!######################################################################
   SUBROUTINE calc_d_ang_nn(vect_nn)
   IMPLICIT NONE
   REAL(kind=8) :: vect_nn(3,4), d(3), ang(3), prod, v(3,3)
   INTEGER :: ia, ib, ic

   DO ia = 1, 3
     DO ib = 1, 3
       v(ib,ia) = vect_nn(ib,ia+1) - vect_nn(ib,1)
     ENDDO !ib
   ENDDO !ia


! lenght
   DO ia = 1, 3
     d(ia) = 0.d0
     DO ib = 1, 3
       d(ia) = d(ia) + v(ib,ia)**2
     ENDDO !ib
     d(ia) = DSQRT( d(ia) )
   ENDDO !ia

! angle
   DO ia = 1, 3
     ib = mod(ia,3) + 1
     prod = 0.d0
     DO ic = 1, 3
       prod = prod + v(ic,ib)*v(ic,ia)
     ENDDO !ic
     ang(ia) = DACOS( prod/(d(ia)*d(ib)) )
   ENDDO !ia

   DO ia = 1, 3
     WRITE(*,'(6f12.6)')d(1), d(2), d(3), ang(1), ang(2), ang(3)
   ENDDO !ia   
   


   RETURN
   END SUBROUTINE calc_d_ang_nn


!######################################################################
! Subroutine vectorial product
!######################################################################
  SUBROUTINE prod_vect(a, b, axb)
  IMPLICIT NONE
  REAL(kind=8) a(3), b(3), axb(3)

!  WRITE(*,*)'prod vect'
!  WRITE(*,*)a(1), a(2), a(3)
!  WRITE(*,*)b(1), b(2), b(3)
!  WRITE(*,*)'fine'

  axb(1) = a(2)*b(3) - a(3)*b(2)
  axb(2) = a(3)*b(1) - a(1)*b(3)
  axb(3) = a(1)*b(2) - a(2)*b(1)

   RETURN
  END SUBROUTINE prod_vect


!
!-----------------------------------------------------------------------
subroutine setupmat (q,dyn,nat,nax,at,bg,tau,itau_blk,nsc,alat, &
                     dyn_blk,nat_blk,nax_blk,at_blk,bg_blk,tau_blk,omega_blk, &
                     epsil,zeu,frc,nr1,nr2,nr3,nrx1,nrx2,nrx3,has_zstar,rws,nrws)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  implicit none
  real(kind=8), parameter :: tpi=2.d0*3.14159265358979d0
  !
  ! I/O variables
  !
  integer:: nr1, nr2, nr3, nrx1,nrx2,nrx3, nax, nat, nat_blk, nax_blk, &
            nsc, nrws, itau_blk(nat)
  real(kind=8) :: q(3), tau(3,nax), at(3,3), bg(3,3), alat,      &
                  epsil(3,3), zeu(3,3,nax_blk), rws(0:3,nrws),   &
                  frc(nrx1,nrx2,nrx3,3,3,nax_blk,nax_blk)
  real(kind=8) :: tau_blk(3,nax_blk), at_blk(3,3), bg_blk(3,3), omega_blk
  complex(kind=8) dyn_blk(3,3,nax_blk,nax_blk)
  complex(kind=8) ::  dyn(3,3,nax,nax)
  logical has_zstar
  !
  ! local variables
  !
  real(kind=8) :: arg
  complex(kind=8) :: cfac(nat)
  integer :: i,j,k, na,nb, na_blk, nb_blk, iq
  real(kind=8) qp(3), qbid(3,nsc) ! automatic array
  !
  !
  call q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !
  do iq=1,nsc
     !
     do k=1,3
        qp(k)= q(k) + qbid(k,iq)
     end do
     !
     dyn_blk(:,:,:,:) = (0.d0,0.d0)
     call frc_blk (nax_blk,dyn_blk,qp,tau_blk,nat_blk,              &
                         nr1,nr2,nr3,nrx1,nrx2,nrx3,frc,at_blk,bg_blk,rws,nrws)
!     if (has_zstar) &
!          call rgd_blk(nax_blk,nat_blk,dyn_blk,qp,tau_blk,   &
!                       epsil,zeu,bg_blk,omega_blk,+1.d0)
     !
     do na=1,nat
        na_blk = itau_blk(na)
        do nb=1,nat
           nb_blk = itau_blk(nb)
           !
           arg=tpi* ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
           !
           cfac(nb) = cmplx(cos(arg),sin(arg))/nsc
           !
        end do ! nb
        !
        do i=1,3
           do j=1,3
              !
              do nb=1,nat
                 nb_blk = itau_blk(nb)
                 dyn(i,j,na,nb) = dyn(i,j,na,nb) + cfac(nb) * &
                      dyn_blk(i,j,na_blk,nb_blk)
              end do ! nb
              !
           end do ! j
        end do ! i
     end do ! na
     !
  end do ! iq
  !
  return
end subroutine setupmat


!-----------------------------------------------------------------------
subroutine frc_blk(nax,dyn,q,tau,nat,                             &
                   nr1,nr2,nr3,nrx1,nrx2,nrx3,frc,at,bg,rws,nrws)
  !-----------------------------------------------------------------------
  ! calculates the dynamical matrix at q from the (short-range part of the)
  ! force constants 
  !
  implicit none
  integer nr1, nr2, nr3, nrx1,nrx2,nrx3, nax, nat, n1, n2, n3, &
          ipol, jpol, na, nb, m1, m2, m3, nint, i,j, nrws
  complex(kind=8) dyn(3,3,nax,nax), cmplx
  real(kind=8) frc(nrx1,nrx2,nrx3,3,3,nax,nax), tau(3,nax), q(3), arg, &
               at(3,3), bg(3,3), r(3), weight, r_ws(3),  &
               total_weight, rws(0:3,nrws)
  real(kind=8), parameter:: tpi = 2.0*3.14159265358979d0
  real(kind=8), external :: wsweight
  !



  do na=1, nat
     do nb=1, nat
        total_weight=0.0d0
        do n1=-2*nrx1,2*nrx1
           do n2=-2*nrx2,2*nrx2
              do n3=-2*nrx3,2*nrx3
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 do i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                    r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                 end do

                 weight = wsweight(r_ws,rws,nrws)
                 if (weight .gt. 0.0) then
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = mod(n1+1,nr1)
                    if(m1.le.0) m1=m1+nr1
                    m2 = mod(n2+1,nr2)
                    if(m2.le.0) m2=m2+nr2
                    m3 = mod(n3+1,nr3)
                    if(m3.le.0) m3=m3+nr3
                    !
                    ! FOURIER TRANSFORM
                    !
                    arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                    do ipol=1, 3
                       do jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               frc(m1,m2,m3,ipol,jpol,na,nb)     &
                               *cmplx(cos(arg),-sin(arg))*weight
                       end do
                    end do
                 end if
                 total_weight=total_weight + weight
              end do
           end do
        end do
        if (abs(total_weight-nr1*nr2*nr3).gt.1.0d-8) then
           write(*,*) total_weight
           write(*,*)'ERRORE  frc_blk wrong total_weight'
        end if
     end do
  end do

  !
  return
end subroutine frc_blk
!

!
!-----------------------------------------------------------------------
subroutine q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !-----------------------------------------------------------------------
  ! generate list of q (qbid) that are G-vectors of the supercell
  ! but not of the bulk
  !
  implicit none
  integer :: nsc
  real(kind=8) qbid(3,nsc), at_blk(3,3), bg_blk(3,3), at(3,3), bg(3,3)
  !
  integer, parameter:: nr1=4, nr2=4, nr3=4, &
                       nrm=(2*nr1+1)*(2*nr2+1)*(2*nr3+1)
  real(kind=8), parameter:: eps=1.0e-7
  integer :: i, j, k,i1, i2, i3, idum(nrm), iq
  real(kind=8) :: qnorm(nrm), qbd(3,nrm) ,qwork(3), delta
  logical lbho
  !
  i = 0
  do i1=-nr1,nr1
     do i2=-nr2,nr2
        do i3=-nr3,nr3
           i = i + 1
           do j=1,3
              qwork(j) = i1*bg(j,1) + i2*bg(j,2) + i3*bg(j,3)
           end do ! j
           !
           qnorm(i)  = qwork(1)**2 + qwork(2)**2 + qwork(3)**2
           !
           do j=1,3
              !
              qbd(j,i) = at_blk(1,j)*qwork(1) + &
                         at_blk(2,j)*qwork(2) + &
                         at_blk(3,j)*qwork(3)
           end do ! j
           !
           idum(i) = 1
           !
        end do ! i3
     end do ! i2
  end do ! i1
  !
  do i=1,nrm-1
     if (idum(i).eq.1) then
        do j=i+1,nrm
           if (idum(j).eq.1) then
              lbho=.true.
              do k=1,3
                 delta = qbd(k,i)-qbd(k,j)
                 lbho = lbho.and. (abs(nint(delta)-delta).lt.eps)
              end do ! k
              if (lbho) then
                 if(qnorm(i).gt.qnorm(j)) then
                    qbd(1,i) = qbd(1,j)
                    qbd(2,i) = qbd(2,j)
                    qbd(3,i) = qbd(3,j)
                    qnorm(i) = qnorm(j)
                 end if
                 idum(j) = 0
              end if
           end if
        end do ! j
     end if
  end do ! i
  !
  iq = 0
  do i=1,nrm
     if (idum(i).eq.1) then
        iq=iq+1
        qbid(1,iq)= bg_blk(1,1)*qbd(1,i) +  &
                    bg_blk(1,2)*qbd(2,i) +  &
                    bg_blk(1,3)*qbd(3,i)
        qbid(2,iq)= bg_blk(2,1)*qbd(1,i) +  &
                    bg_blk(2,2)*qbd(2,i) +  &
                    bg_blk(2,3)*qbd(3,i)
        qbid(3,iq)= bg_blk(3,1)*qbd(1,i) +  &
                    bg_blk(3,2)*qbd(2,i) +  &
                    bg_blk(3,3)*qbd(3,i)
     end if
  end do ! i
  !
  if (iq.ne.nsc) WRITE(*,*)'ERRORE q_gen probably nr1,nr2,nr3 too small '
  return
end subroutine q_gen
!
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine wsinit(rws,nrwsx,nrws,atw)
!-----------------------------------------------------------------------
!
!  USE kinds, only : DP
  implicit none
  integer i, ii, ir, jr, kr, nrws, nrwsx, nx
  real(kind=8) rt, eps, rws(0:3,nrwsx), atw(3,3)
  parameter (eps=1.0d-6,nx=2)
  ii = 1
  do ir=-nx,nx
     do jr=-nx,nx
        do kr=-nx,nx
           do i=1,3
              rws(i,ii) = atw(i,1)*ir + atw(i,2)*jr + atw(i,3)*kr
           end do
           rws(0,ii)=rws(1,ii)*rws(1,ii)+rws(2,ii)*rws(2,ii)+            &
                               rws(3,ii)*rws(3,ii)
           rws(0,ii)=0.5d0*rws(0,ii)
           if (rws(0,ii).gt.eps) ii = ii + 1
           if (ii.gt.nrwsx) WRITE(*,*)'ERRORE wsinit ii.gt.nrwsx'
        end do
     end do
  end do
  nrws = ii - 1
  return
end subroutine wsinit
!
!-----------------------------------------------------------------------
function wsweight(r,rws,nrws)
!-----------------------------------------------------------------------
!
!  USE kinds, only : dp
  implicit none
  integer ir, nreq, nrws
  real(kind=8) r(3), rrt, ck, eps, rws(0:3,nrws), wsweight
  parameter (eps=1.0d-6)
!
  wsweight = 0.d0
  nreq = 1
  do ir =1,nrws
     rrt = r(1)*rws(1,ir) + r(2)*rws(2,ir) + r(3)*rws(3,ir)
     ck = rrt-rws(0,ir)
     if ( ck .gt. eps ) return
     if ( abs(ck) .lt. eps ) nreq = nreq + 1
  end do
  wsweight = 1.d0/dble(nreq)
  return
end function wsweight

!---------------------------------------------------------------------

subroutine recips (a1, a2, a3, b1, b2, b3)
  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1,b2,b3
  !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
  !
  !     first the input variables
  !
!  use kinds, ONLY: DP
  implicit none
  real(kind=8) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
  ! input: first direct lattice vector
  ! input: second direct lattice vector
  ! input: third direct lattice vector
  ! output: first reciprocal lattice vector
  ! output: second reciprocal lattice vector
  ! output: third reciprocal lattice vector
  !
  !   then the local variables
  !
  real(kind=8) :: den, s
  ! the denominator
  ! the sign of the permutations
  integer :: iperm, i, j, k, l, ipol
  ! counter on the permutations
  !\
  !  Auxiliary variables
  !/
  !
  ! Counter on the polarizations
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.d0
100 do iperm = 1, 3
     den = den + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s
  if (s.lt.0.d0) goto 100
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
     b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
     b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
     l = i
     i = j
     j = k
     k = l
  enddo
  return
end subroutine recips





!-----------------------------------------------------------------------------
  subroutine cdiagh2 (n,h,ldh,e,v)
!-----------------------------------------------------------------------------
  !
  !   calculates all the eigenvalues and eigenvectors of a complex
  !   hermitean matrix H . On output, the matrix is unchanged
  !
 implicit none
 !
 ! on INPUT
 integer          n,       &! dimension of the matrix to be diagonalized
      &           ldh       ! leading dimension of h, as declared
 ! in the calling pgm unit
 complex(kind=8)  h(ldh,n)  ! matrix to be diagonalized
 !
 ! on OUTPUT
 real(kind=8)     e(n)      ! eigenvalues
 complex(kind=8)  v(ldh,n)  ! eigenvectors (column-wise)
 !
 ! LOCAL variables (LAPACK version)
 !
 integer          lwork,   &! aux. var.
      &           ILAENV,  &! function which gives block size
      &           nb,      &! block size
      &           info      ! flag saying if the exec. of libr. routines was ok
 !
 real(kind=8), allocatable::   rwork(:)
 complex(kind=8), allocatable:: work(:)
 !
 !     check for the block size
 !
 nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
 if (nb.lt.1) nb=max(1,n)
 if (nb.eq.1.or.nb.ge.n) then
    lwork=2*n-1
 else
    lwork = (nb+1)*n
 endif
 !
 ! allocate workspace
 !
 call ZCOPY(n*ldh,h,1,v,1)
 allocate(work (lwork))
 allocate(rwork (3*n-2))
 call ZHEEV('V','U',n,v,ldh,e,work,lwork,rwork,info)
 IF(info.ne.0)WRITE(*,*)'ATT. info=',info
! call errore ('cdiagh2','info =/= 0',abs(info))
 ! deallocate workspace
 deallocate(rwork)
 deallocate(work)
 !
 return
end subroutine cdiagh2


