        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:01 2020
        MODULE TRACE3D__genmod
          INTERFACE 
            SUBROUTINE TRACE3D(Q,DQ,QM,QP,DX,DY,DZ,DT,NGRID)
              REAL(KIND=8) :: Q(1:32,-1:4,-1:4,-1:4,1:8)
              REAL(KIND=8) :: DQ(1:32,-1:4,-1:4,-1:4,1:8,1:3)
              REAL(KIND=8) :: QM(1:32,-1:4,-1:4,-1:4,1:8,1:3)
              REAL(KIND=8) :: QP(1:32,-1:4,-1:4,-1:4,1:8,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: DZ
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE TRACE3D
          END INTERFACE 
        END MODULE TRACE3D__genmod
