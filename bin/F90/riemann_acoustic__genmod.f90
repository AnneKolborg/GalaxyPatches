        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:02 2020
        MODULE RIEMANN_ACOUSTIC__genmod
          INTERFACE 
            SUBROUTINE RIEMANN_ACOUSTIC(QLEFT,QRIGHT,FGDNV,NGRID)
              REAL(KIND=8) :: QLEFT(1:32,1:8)
              REAL(KIND=8) :: QRIGHT(1:32,1:8)
              REAL(KIND=8) :: FGDNV(1:32,1:9)
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE RIEMANN_ACOUSTIC
          END INTERFACE 
        END MODULE RIEMANN_ACOUSTIC__genmod
