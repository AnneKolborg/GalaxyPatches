        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 17 10:25:11 2020
        MODULE GET_CELL_INDEX_FAST__genmod
          INTERFACE 
            SUBROUTINE GET_CELL_INDEX_FAST(INDP,CELL_LEV,XPART,IND_GRID,&
     &NBORS_FATHER_CELLS,NP,ILEVEL)
              INTEGER(KIND=4) :: INDP(1:99)
              INTEGER(KIND=4) :: CELL_LEV(1:99)
              REAL(KIND=8) :: XPART(1:99,1:3)
              INTEGER(KIND=4) :: IND_GRID
              INTEGER(KIND=4) :: NBORS_FATHER_CELLS(1:27)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: ILEVEL
            END SUBROUTINE GET_CELL_INDEX_FAST
          END INTERFACE 
        END MODULE GET_CELL_INDEX_FAST__genmod
