C     METROPOLIS SAMPLE GENERATOR
C     0!
C     Lluís Torres 
C     TFM
C     FORTRAN 95

      PROGRAM SAMPLE_GENERATOR

      USE MODEL

C-----(SYSTEM)------------------------------------------------
C     NODES, EDGES, CONNECTIVITY
      INTEGER N, M, z
C     R (:= m TROTTER INDEX)
      INTEGER R
C     +1 -1 EDGES RATIO (1 => ALL +1), (0 => ALL -1)
      REAL*8 p
C     TEMPERATURE (TEMP := k_B·T)
      REAL*8 TEMP
C     H (:= Γ = TRANSVERSE FIELD)
      REAL*8 H
C     p-LIST, TEMP-LIST, Γ-LIST
      INTEGER p_SIZE, TEMP_SIZE, H_SIZE
      REAL*8,ALLOCATABLE:: p_LIST(:), TEMP_LIST(:), H_LIST(:)
C     PRINCIPAL ARRAYS
      TYPE(MULTI_ARRAY),ALLOCATABLE:: NBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: INBR(:)
      TYPE(MULTI_ARRAY),ALLOCATABLE:: JJ(:)
C-----(SIMULATION)---------------------------------------------
C     SAMPLE SIZE (# OF SPIN CONFIGURATIONS)
      INTEGER C
C     TOTAL MONTE-CARLO STEPS (MCS)
      INTEGER MCTOT
C     MCS TILL WE CONSIDER EQUILIBRIUM
      INTEGER MCINI
C     SAVE SPIN CONFIGURATIONS EVERY SC (MCS)
      INTEGER SC
C     NUMBER OF GRAPHS TO SIMULATE FOR EVERY P VALUE
      INTEGER NSEEDS
C     SEED NUMBER, INITIAL SEED NUMBER
      INTEGER SEED, SEEDini
      PARAMETER(SEEDini = 100)
C     RANDOM NUMBER GENERATOR
      EXTERNAL r1279
C     ESTIMATE TIME VARIABLES
      REAL*4 TIME1, TIME2, time
C     SIMULATION VARIABLES
      INTEGER, ALLOCATABLE:: S1(:,:), S2(:,:)
      LOGICAL valid1, valid2
      REAL*8 DE1, DE2
      REAL*8 ENE1, ENE2
C-----(SPIN CONFIGURATION SAVING VARIABLES)-------------------
C     STORE SPIN CONFIGURATION AS N/zip_size INTEGERS
      INTEGER zip_size
      INTEGER, ALLOCATABLE:: bin1(:), bin2(:)
      INTEGER, ALLOCATABLE:: decimal1(:), decimal2(:)
      INTEGER, ALLOCATABLE:: array1(:), array2(:)
C-----(ORDER PARAMETER)---------------------------------------
C     EDWARDS-ANDERSON ORDER PARAMETER
      REAL*8 Q, QSELF
      INTEGER, ALLOCATABLE:: S1K(:,:), S2K(:,:)
      INTEGER loc
C-----(DUMMY)-------------------------------------------------
      INTEGER ITEMP, IH, Ip 
      INTEGER IMC
      CHARACTER(4) str
      CHARACTER(3) str1, str2, str3, str4
      INTEGER, ALLOCATABLE:: CHECK(:)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER TAU !NOT USED IN THIS PROGRAM

C-----------------------------------------------------------------------
C     START
C-----------------------------------------------------------------------

      PRINT*, 'METROPOLIS SAMPLE GENERATION'

C***********************************************************************
C     READ SIMULATION VARIABLES FROM INPUT FILE
      CALL READ_INPUT(N,z,R,TEMP_SIZE,TEMP_LIST,H_SIZE,H_LIST,
     . p_SIZE,p_LIST,C,MCINI,NSEEDS,SC,zip_size,TAU)
      M = z*N/2
      MCTOT = MCINI + C*SC/2
C     ALLOCATION
      ALLOCATE(decimal1(1:N/zip_size))
      ALLOCATE(decimal2(1:N/zip_size))
      ALLOCATE(array1(1:N))
      ALLOCATE(array2(1:N))
      ALLOCATE(bin1(1:N))
      ALLOCATE(bin2(1:N))
      ALLOCATE(S1(1:R,1:N))
      ALLOCATE(S2(1:R,1:N))
      
      ALLOCATE(S1K(1:R,1:N))
      ALLOCATE(S2K(1:R,1:N))

      ALLOCATE(CHECK(1:10))
      CHECK = (/0,10,30,100,300,1000,3000,10000,30000,100000/)

C***********************************************************************
      CALL CPU_TIME(TIME1)
C***********************************************************************

C     FOR ALL TEMP VALUES
      DO ITEMP = 1,TEMP_SIZE
      TEMP = TEMP_LIST(ITEMP)
      WRITE(str,'(f4.2)') TEMP
      str1 = str(1:1)//str(3:4)

C     FOR ALL Γ VALUES      
      DO IH = 1,H_SIZE
      H = H_LIST(IH)
      WRITE(str,'(f4.2)') H
      str2 = str(1:1)//str(3:4)

C***********************************************************************
C     CREATE DIRECTORY FOR EACH TEMP AND Γ VALUE IN THE 'sample' FOLDER
      CALL SYSTEM('mkdir -p results/sample/T'//str1//'_Γ'//str2)
C     COPY INPUT FILE TO THE RESULTS FOLDER
      CALL SYSTEM('cp input.txt results/input.txt')
C     CREATE DIRECTORY FOR EACH TEMP AND Γ VALUE IN THE 'order' FOLDER
      CALL SYSTEM('mkdir -p results/order/T'//str1//'_Γ'//str2)
C***********************************************************************

C     FOR ALL p VALUES
      DO Ip = 1,p_SIZE
      p = p_LIST(IP)
      WRITE(str,'(f4.2)') p
      str3 = str(1:1)//str(3:4)

C     FOR ALL SEEDS
      DO SEED = SEEDini,SEEDini+NSEEDS-1
      WRITE(str4,'(i3)') SEED
      CALL setr1279(SEED)

C***********************************************************************
C     INITIAL RANDOM SYSTEM (GRAPH+COUPLINGS)
      CALL IRS(N,M,p,NBR,INBR,JJ)
C***********************************************************************
C     SPIN CONFIGURATION FILE FOR EACH p VALUE AND SEED
      OPEN(UNIT=1,FILE='results/sample/T'//str1//'_Γ'//str2//
     .'/S_'//str3//'_'//str4//'.bin',FORM='UNFORMATTED')
C     ORDER PARAMETER FILE FOR EACH p VALUE AND SEED
      OPEN(UNIT=2,FILE='results/order/T'//str1//'_Γ'//str2//
     .'/q_'//str3//'_'//str4//'.dat')
C***********************************************************************
C     GENERATION OF TWO RANDOM INITIAL SPIN CONFIGURATIONS
      DO i = 1,R
            DO j = 1,N
                  S1(i,j) = INT(2*MOD(INT(2*r1279()),2) - 1)
                  S2(i,j) = INT(2*MOD(INT(2*r1279()),2) - 1)
            END DO
      END DO
C***********************************************************************
C     INITIAL ENERGY
      ENE1 = ENERG(N,R,S1,TEMP,H,NBR,JJ)
      ENE2 = ENERG(N,R,S2,TEMP,H,NBR,JJ)
C***********************************************************************
C     INITIAL ORDER PARAMETER
      S1K = S1
      S2K = S2

      Q = ABS(OP(R,N,S1,S2)) !OVERLAP BETWEEN TWO DIFFERENT REALIZATIONS
      QSELF = 0.5D0*(ABS(OP(R,N,S1,S1K))+ABS(OP(R,N,S2,S2K))) !OVERLAP WITH ITSELF

      WRITE(2,*) 0, Q, QSELF

      loc = 2
C***********************************************************************
C     MONTE-CARLO SIMULATION
      DO IMC = 1,MCTOT
            DO IPAS = 1,N*R
                  CALL METROPOLIS(N,R,S1,valid1,TEMP,H,DE1,NBR,JJ)
                  CALL METROPOLIS(N,R,S2,valid2,TEMP,H,DE2,NBR,JJ)
                  IF (valid1) THEN
                        ENE1 = ENE1 + DE1
                  END IF
                  IF (valid2) THEN
                        ENE2 = ENE2 + DE2
                  END IF
            END DO
C           EXTRACT THE SPIN CONFIGURATION EVERY SC MONTE-CARLO STEPS
            IF ((IMC.GT.MCINI).AND.(SC*(IMC/SC).EQ.IMC)) THEN
                  DO i = 1,R
                        array1 = S1(i,:)
                        CALL ARRAY2BIN(N,bin1,array1)
                        CALL BIN2DEC(N,zip_size,bin1,decimal1)
                        WRITE(1) decimal1
                  END DO
                  DO i = 1,R
                        array2 = S2(i,:)
                        CALL ARRAY2BIN(N,bin2,array2)
                        CALL BIN2DEC(N,zip_size,bin2,decimal2)
                        WRITE(1) decimal2
                  END DO
            END IF

C     EXTRACT THE ORDER PARAMETER EVERY CHECK MONTE-CARLO STEPS
      IF (IMC.EQ.CHECK(loc)) THEN
            GAP = CHECK(loc)-CHECK(loc-1)
            WRITE(2,*) IMC, Q/GAP, QSELF/GAP
            loc = loc + 1
            S1K = S1
            S2K = S2
            Q = ABS(OP(R,N,S1,S2)) 
            QSELF = 0.5D0*(ABS(OP(R,N,S1,S1K))+ABS(OP(R,N,S2,S2K)))
      END IF

      Q = Q + ABS(OP(R,N,S1,S2)) !OVERLAP BETWEEN TWO DIFFERENT REALIZATIONS
      QSELF = QSELF + 0.5D0*(ABS(OP(R,N,S1,S1K))+ABS(OP(R,N,S2,S2K))) !OVERLAP WITH ITSELF

      END DO !IMC
C***********************************************************************
      CLOSE(1)
      CLOSE(2)
C     DELLOCATE ARRAYS
      DO i = 1,N
            DEALLOCATE(NBR(i)%v)
            DEALLOCATE(INBR(i)%v)
            DEALLOCATE(JJ(i)%v)
      END DO
      DEALLOCATE(NBR)
      DEALLOCATE(INBR)
      DEALLOCATE(JJ)
C***********************************************************************
C     ESTIMATE EXECUTION TIME
200   FORMAT (A,I4,A,I3,A,I3,A,I3,A,I3,A,I3)
      IF ((SEED.EQ.SEEDini).AND.(p.EQ.p_LIST(1)).AND. 
     .(TEMP.EQ.TEMP_LIST(1)).AND.(H.EQ.H_LIST(1))) THEN
      CALL CPU_TIME(TIME2)
      time = (TIME2-TIME1)*NSEEDS*H_SIZE*p_SIZE*TEMP_SIZE
      WRITE(*,200) "ESTIMATED TIME: ", INT(time/3600), ' h',
     . INT((time/3600-INT(time/3600))*60), ' min', 
     . INT((time/60-INT(time/60))*60), ' s'
      END IF
C***********************************************************************

      END DO !SEED
      END DO !Ip
      END DO !IH
      END DO !ITEMP

C***********************************************************************
      CALL CPU_TIME(TIME2)
      time = (TIME2-TIME1)
      WRITE(*,200) "CPU TIME: ", INT(time/3600), ' h',
     . INT((time/3600-INT(time/3600))*60), ' min', 
     . INT((time/60-INT(time/60))*60), ' s'
C***********************************************************************

      END PROGRAM SAMPLE_GENERATOR