C LOCAL ARRAYS
C     EELAS - ELASTIC STRAINS
C     EPLAS - PLASTIC STRAINS
C     ALPHA - BACK STRESS
C     OLDS - STRESS OF START OF INCREMENT
C     OLDPL - PLASTIC STRAINS AT STRAT OF INCREMENT
C     OLDALPHA - BACK STRESS AT STRAT OF INCREMENT

C     P - EQUATION(2)
C     ANISO - IN EQUATION(1) TO GET F
C     A(3,1) - D F / D STRESS
C     DA(3,3) - D A / D STRESS

C     DELT_EPLAS(3) - INCREMENT OF PLASTIC STRAIN IN (5)
C     R_SIGMA(3,1) - EQUATION(17)
C     R_ALPHA(3,1) - EQUATION(18)     
C     D1(3,3) - FOR (24) (25) (26) ITERATION 


C     D2(3,3) - FOR (24) (25) (26) ITERATION
C     D3(3,3) - FOR (24) (25) (26) ITERATION 
C     BIG_FAI(3,3) - FOR (24) (25) (26) ITERATION  
C     HIN_BIG_FAI(3,3) - 
C     BIG_GAMMA(3) - FOR (24) (25) (26) ITERATION   
C     BIG_PI(3,3) - FOR (24) (25) (26) ITERATION  
C     R_HAT(3,1) - FOR (24) (25) (26) ITERATION
C     Q_WAVE - FOR (24) (25) (26) //JUST A SCALAR   

C     BIG_SIGMA(3,3) - FOR (36) DDSDDE
C     HIN_BIG_SIGMA(3,3)
C     BIG_R(3,3) - FOR (36) DDSDDE      
C     BIG_GAMMA(3,1) - FOR (36) DDSDDE      
C     BIG_KESI(1,1) - FOR (36) DDSDDE 

C     TEMP1(1,1) - FOR TEMP USE
C     TEMP31(3,1) - FOR TEMP USE
C     TEMP33(3,3) - FOR TEMP USE
C     HIN_TEMP33(3,3)
C     ONES(3,3) - IDENTITY MATRIX

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,    
     2CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,    
     3CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C      INTEGER L    !if this step is plastic

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4JSTEP(4)

      DIMENSION EELAS(3), EPLAS(3), ALPHA(3), OLDS(3), OLDPL(3),
     1OLDALPHA(3), P(3,3), ANISO(3,3), A(3,1), DA(3,3), DELT_EPLAS(3,1),
     2R_SIGMA(3,1), R_ALPHA(3,1), D1(3,3), D2(3,3), D3(3,3), 
     3BIG_FAI(3,3), BIG_GAMMA(3,1), BIG_PI(3,3),R_HAT(3,1), 
     4BIG_SIGMA(3,3),BIG_R(3,3), BIG_LAMDA(3,1), BIG_KESI(1,1), 
     5HIN_BIG_FAI(3,3),HIN_BIG_SIGMA(3,3),DOT_STRESS(3,1),ONES(3,3),
     6DOT_ALPHA(3,1),TEMP1(1,1),TEMP31(3,1),TEMP33(3,3),HIN_TEMP33(3,3),
     7B(3,3),B_INV(3,3)
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1ENUMAX=.4999D0, TOLER=1.0D-6)
C
C ----------------------------------------------------------------
C UMAT FOR ANISOTROPIC NONLINER KINEMATIC HARDENING PLASTICITY 
C USED FOR PLANE STRESS
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C PROPS(3) - SIGMA_0 
C PROPS(4) - R_0
C PROPS(5) - R_45
C PROPS(6) - R_90
C PROPS(7) - C
C PROPS(8) - GAMMA
C ----------------------------------------------------------------

C*******************JUST FOR TEST******************
      OPEN(unit=79,file='F:\JDGwork\ABAQUS\bugcheckoflinear.txt',
     1 status='unknown')
      !write(*,*)'PLZ INPUT AN INTEGER:'
      !READ(*,*) KTEMPREAD
      WRITE(79,*)'*****************************************************'
      WRITE(79,*)'******************* NEW UMAT ************************'
      WRITE(79,*)'*****************************************************'

C ----------------------------------------------------------------
C CHECK FOR ELEMENT TYPE
C
      IF (NDI.NE.2) THEN 
          WRITE (7, *) 'THIS UMAT MAY ONLY BE USED FOR ELEMENTS     
     1    WITH PLANE STRESS' 
          CALL XIT 
      ENDIF 
      IF (NTENS.NE.3) THEN 
          WRITE (7, *) 'THIS UMAT MAY ONLY BE USED FOR ELEMENTS     
     1    WITH PLANE STRESS' 
          CALL XIT 
      ENDIF 
C
C ELASTIC PROPERTIES
      EMOD=PROPS(1)
      ENU=PROPS(2)
C
C      WRITE(79,*)'***** KINC *****'
C      WRITE(79,*)KINC  
C      WRITE(79,"(I5)") KINC
C105   FORMAT(' ***APPEARED:', I5)
      WRITE(79,*)'***** OLD DDSDDE *****'
      WRITE(79,'(9(2X, ES13.6))')DDSDDE      
C ELASTIC STIFFNESS !!!HERE PLANE STRESS
      DDSDDE = 0
      DDSDDE(1,1) = EMOD / (ONE - ENU * ENU)
      DDSDDE(2,2) = DDSDDE(1,1)
      DDSDDE(1,2) = ENU * EMOD / (ONE - ENU * ENU)
      DDSDDE(2,1) = DDSDDE(1,2)
      DDSDDE(3,3) = EMOD / (ONE - ENU * ENU) * (ONE - ENU) / TWO
C !!!!!HERE WE USE GAMMA33 INSTEAD OF EPSILON33, IS IT RIGHT?
C
      WRITE(79,*)'***** DDSDDE ELASTIC*****'
      WRITE(79,'(9(2X, ES13.6))')DDSDDE
      WRITE(79,*)'***** OLD STRESS *****'
      WRITE(79,'(3(2X, ES13.6))')(STRESS(K1), K1 = 1,3)
      WRITE(79,*)'***** OLD STRAN *****'
      WRITE(79,'(3(2X, ES13.6))')(STRAN(K1), K1 = 1,3)

C RECOVER ELASTIC STRAINS, PLASTIC STRAINS, BACK STRESSES, EQUIVALENT PLASTIC STRAIN AND ROTATE
C NOTE: USE CODE 1 FOR (TENSOR) STRESS, CODE 2 FOR (ENGINEERING) STRAIN ???
C !!!!!HERE WE SHOULD CHECK THE USE OF 1 OR 2 AND NDI
      WRITE(79,*)'***** OLD STATEV *****'
      WRITE(79,'(11(2X, ES13.6))')(STATEV(K1), K1 = 1, 11)
      !write(*,*)'***1***'
      !CALL ROTSIG(STATEV(1), DROT, EELAS, 2, NDI, NSHR)
      !write(*,*)'***2***'
      !CALL ROTSIG(STATEV(4), DROT, EPLAS, 2, NDI, NSHR)
      !write(*,*)'***3***'
      !CALL ROTSIG(STATEV(7), DROT, ALPHA, 1, NDI, NSHR)
          DO K1=1,3
              EELAS(K1)=STATEV(K1)
              EPLAS(K1)=STATEV(K1+3)
              ALPHA(K1) =STATEV(K1+6)
          END DO
      EQEPLAS = STATEV(10)
C IF WE WANT TO CALCULATE PLASTIC STRAIN ENERGY
      SPD = STATEV(11)

C*******************JUST FOR TEST******************


C SAVE STRESS AND PLASTIC STRAINS AND
C CALCULATE PREDICTOR STRESS AND ELASTIC STRAIN
C
      DO K1 = 1, 3
          OLDS(K1) = STRESS(K1)
          OLDPL(K1) = EPLAS(K1)
          OLDALPHA(K1) = ALPHA(K1)
          EELAS(K1) = EELAS(K1) + DSTRAN(K1)
      ENDDO
      DO K1=1,3
          DO K2 = 1, 3
              STRESS(K2) = STRESS(K2) + DDSDDE(K2, K1) * DSTRAN(K1)
          END DO
      END DO 

      WRITE(79,*)'***** OLDS *****'
      WRITE(79,'(2X, ES13.6)') OLDS(1)
      WRITE(79,'(2X, ES13.6)') OLDS(2)
      WRITE(79,'(2X, ES13.6)') OLDS(3)

      WRITE(79,*)'***** DSTRAN *****'
      WRITE(79,'(3(2X, ES13.6))')
     1(DSTRAN(K1), K1 = 1, 3)
      WRITE(79,*)'***** PREDICT EELAS *****'
      WRITE(79,'(3(2X, ES13.6))')
     1(EELAS(K1), K1 = 1, 3)
      WRITE(79,*)'***** PREDICT STRESS *****'
      WRITE(79,'(3(2X, ES13.6))')
     1(STRESS(K1), K1 = 1, 3)

C THIS IS ALL OF THE ELASTIC PART, AND THE STRESS HAS BEEN UPDATES AS TEST STRESS
C ----------------------------------------------------------------      


      ONES=ZERO
      DO K1 = 1, 3
          ONES(K1,K1) = ONE
      END DO


C ----------------------------------------------------------------        
C PLASCIT PARAMETER      
      SIGMA_0 = PROPS(3) 
      R_0 = PROPS(4) 
      R_45 = PROPS(5) 
      R_90 = PROPS(6)
      C = PROPS(7)
      GAMMA = PROPS(8)

C ----------------------------------------------------------------      
C CALCULATE P USE (2) AND THEN CALCULATE ANISO USE (2)
      P_12 = R_0 / (1 + R_0)
      P_22 = R_0 * (1 + R_90) / R_90 / (1 + R_0)
      P_66 = (R_0 + R_90) * (1 + 2 * R_45) / R_90 / (1 + R_0)
      P = ZERO
      P(1,1) = ONE
      P(1,2) = -P_12
      P(2,1) = -P_12
      P(2,2) = P_22
      P(3,3) = P_66
      DO K1 = 1, 3
          DO K2 = 1, 3
              ANISO(K1,K2) = P(K1, K2)
          END DO
      END DO     
      
C----------------- JUST FOR TESTING THE ANISO POROPERTY
      !ANISO=0.D0
      !ANISO(1,1)=1.D0
      !ANISO(2,2)=13.D0/7.D0
      !ANISO(1,2)=-8.D0/7.D0
      !ANISO(2,1)=ANISO(1,2)
      !ANISO(3,3)=18.D0/7.D0
      
      
      

C ----------------------------------------------------------------      
C CALCULATE F USE (1)  --1*1, HERE AHPHA IS THE OLD ALPHA
      DO K1=1,3
          TEMP31(K1,1)=STRESS(K1)-ALPHA(K1)
      ENDDO
      TEMP1=MATMUL(MATMUL(TRANSPOSE(TEMP31),ANISO),TEMP31)
      F=SQRT(ONE/TWO*TEMP1(1,1))    

      write(*,*)'***4***'
      WRITE(79,*)'***** PREDICT F *****'
      write(79,'(2X, ES13.6)')F

C--------------------------INITIATE--------------------
      
C      DELT_LAMDA=STATEV(10) 
      DELT_EPLAS=0.0
      DELT_LAMDA=0.0

      L=0
      TICK=0
C ---------------------------CYCLE--------------------------------  
C CHECK FOR YIELD USE (16)
      DO WHILE ( (F - SIGMA_0) .GE. TOLER)
          TICK=TICK+1
          L=1
          write(79,*)'***************** NEW CYCLE STARTS **********
     1    **********'
C ----------------------------------------------------------------      
C CALCULATE A = DOT F / DOT STRESS  --3*1   //A CHANGES AS STRESS CHANGE
          write(*,*)'***5***'
          DO K1=1,3
              TEMP31(K1,1)=STRESS(K1)-ALPHA(K1)
          ENDDO
          A=MATMUL(ANISO,TEMP31)/F/TWO                    !具体公式推导见JDG附件
          WRITE(79,*)'***** A *****'
          WRITE(79,'(3(2X, ES13.6))')(A(K1,1),K1=1,3)

C CALCULATE DA = DOT A / DOT STRESS  --3*3   
          DA=(ANISO-2*MATMUL(A,TRANSPOSE(A)))/F/TWO          !具体公式推导见JDG附件
C SET THE VALUE OF PLASTIC STRAIN INCREMENT DELT_EPLAS(3) USE(5) (WE SHOULD FIRST HAVE A AND DELT_LAMDA FIRST)      
          DELT_EPLAS= DELT_LAMDA * A
C SET THE VALUE OF R_SIGMA USE (17)  --3*1
          DO K1=1,3
              TEMP31(K1,1)=DSTRAN(K1)
          ENDDO
          TEMP31 = MATMUL(DDSDDE , TEMP31-DELT_EPLAS)
          DO K1=1, 3
              R_SIGMA(K1,1) = STRESS(K1) - OLDS(K1) - TEMP31(K1,1)
          END DO

          WRITE(79,*)'????????????????'
          WRITE(79,'(3(2X, ES13.6))')(STRESS(K1),K1=1,3)
          WRITE(79,'(3(2X, ES13.6))')(OLDS(K1),K1=1,3)
          WRITE(79,'(3(2X, ES13.6))')(DELT_EPLAS(K1,1),K1=1,3)
          WRITE(79,'(3(2X, ES13.6))')(DSTRAN(K1),K1=1,3)
          WRITE(79,'(3(2X, ES13.6))')(TEMP31(K1,1),K1=1,3)
          WRITE(79,*)'????????????????'

          WRITE(79,*)'***** R_SIGMA *****'
          WRITE(79,'(3(2X, ES13.6))')(R_SIGMA(K1,1),K1=1,3)

C SET THE VALUE OF R_ALPHA USE (18)  --3*1
          DO K1=1, 3
              R_ALPHA(K1,1) = ALPHA(K1) - OLDALPHA(K1) 
     1        - TWO / THREE * C * DELT_EPLAS(K1,1) 
     2        + GAMMA * ALPHA(K1) * DELT_LAMDA
          END DO 
          WRITE(79,*)'*****R_ALPHA *****'
          WRITE(79,'(3(2X, ES13.6))')(R_ALPHA(K1,1),K1=1,3)


C SET THE VALUE OF R_F USE (19)  --SCALAR
          R_F = F - SIGMA_0
          WRITE(79,*)'***** R_F *****'
          WRITE(79,'(2X, ES13.6)')R_F

C ----------------------------------------------------------------        
C CALCULATE D1(3*3),D2(3*3),D3(3*3),BIG_FAI(3*3), BIG_GAMMA(3,1), BIG_PI(3,3), R_HAT(3,1)(r^), Q_WAVE(q~)
C-----------D1------3*3
          D1 = TWO / THREE * C * DELT_LAMDA * DA
          DO K1=1, NTENS
              D1(K1 ,K1) = D1(K1,K1) + ONE + GAMMA * DELT_LAMDA
          END DO
C-----------D2------3*3
      D2 = D1 * TWO / THREE * C
C-----------D3------3*3
      D3 = TWO / THREE * C * DELT_LAMDA * MATMUL(D1,DA)
C-----------BIG_FAI------3*3
      BIG_FAI = DELT_LAMDA * MATMUL(DDSDDE,DA)
     1-DELT_LAMDA* MATMUL(MATMUL(DDSDDE,DA),D3)
      DO K1=1,3
          BIG_FAI(K1,K1)=BIG_FAI(K1,K1)+ ONE
      ENDDO
      WRITE(79,*)'*****BIG_FAI*****'
      WRITE(79,'(9(2X, ES13.6))')BIG_FAI
C-----------BIG_GAMMA------3*1
      DO K1=1,3
          TEMP31(K1,1)=ALPHA(K1)
      ENDDO
      BIG_GAMMA = DELT_LAMDA * GAMMA * 
     1MATMUL(MATMUL(MATMUL(DDSDDE,DA),D1),TEMP31)
      WRITE(79,*)'*****BIG_GAMMA*****'
      WRITE(79,'(9(2X, ES13.6))')BIG_GAMMA
C-----------BIG_PI------3*3
      BIG_PI=DDSDDE-DELT_LAMDA*MATMUL(MATMUL(DDSDDE,DA),D2)
      WRITE(79,*)'*****BIG_PAI*****'
      WRITE(79,'(9(2X, ES13.6))')BIG_PI
C-----------R_HAT-------3*1
      R_HAT=R_SIGMA+
     1DELT_LAMDA*MATMUL(MATMUL(MATMUL(DDSDDE,DA),D1),R_ALPHA)
      WRITE(79,*)'***** R_HAT *****'
      WRITE(79,'(3(2X, ES13.6))')(R_HAT(K1,1),K1=1,3)
C-----------Q_WAVE------SCALAR
      
      B=BIG_FAI
            det_B = B(1,1)*(B(2,2)*B(3,3) - B(3,2)*B(2,3)) 
     1-        B(2,1)*(B(1,2)*B(3,3) - B(3,2)*B(1,3)) 
     2+        B(3,1)*(B(1,2)*B(2,3) - B(2,2)*B(1,3))
      det_B_inv = 1.D0/det_B
      B_inv(1,1) = det_B_inv*(B(2,2)*B(3,3)-B(3,2)*B(2,3))
      B_inv(1,2) = det_B_inv*(B(3,2)*B(1,3)-B(1,2)*B(3,3))
      B_inv(1,3) = det_B_inv*(B(1,2)*B(2,3)-B(2,2)*B(1,3))
      B_inv(2,1) = det_B_inv*(B(3,1)*B(2,3)-B(2,1)*B(3,3))
      B_inv(2,2) = det_B_inv*(B(1,1)*B(3,3)-B(3,1)*B(1,3))
      B_inv(2,3) = det_B_inv*(B(2,1)*B(1,3)-B(1,1)*B(2,3))
      B_inv(3,1) = det_B_inv*(B(2,1)*B(3,2)-B(3,1)*B(2,2))
      B_inv(3,2) = det_B_inv*(B(3,1)*B(1,2)-B(1,1)*B(3,2))
      B_inv(3,3) = det_B_inv*(B(1,1)*B(2,2)-B(2,1)*B(1,2))
      TEMP33=B_inv


      TEMP1=MATMUL(MATMUL(TRANSPOSE(A),TEMP33),R_HAT)
     1+MATMUL(MATMUL(MATMUL(TRANSPOSE(A),D3),TEMP33),R_HAT)
     2+MATMUL(MATMUL(TRANSPOSE(A),D1),R_ALPHA)
      Q_WAVE=TEMP1(1,1)
      WRITE(79,*)'***** Q_WAVE *****'
      WRITE(79,'(2X, ES13.6)')Q_WAVE

C ----------------------------------------------------------------        
C ITERATION FOR STRESS, ALPHA AND DELT_LAMDA  HERE TAKE CARE OF (3) OR (3,1)
C CALCULATE DOT_LAMDA(DOUBLE) USE (26) 
      B=BIG_FAI
            det_B = B(1,1)*(B(2,2)*B(3,3) - B(3,2)*B(2,3)) 
     1-        B(2,1)*(B(1,2)*B(3,3) - B(3,2)*B(1,3)) 
     2+        B(3,1)*(B(1,2)*B(2,3) - B(2,2)*B(1,3))
      det_B_inv = 1.D0/det_B
      B_inv(1,1) = det_B_inv*(B(2,2)*B(3,3)-B(3,2)*B(2,3))
      B_inv(1,2) = det_B_inv*(B(3,2)*B(1,3)-B(1,2)*B(3,3))
      B_inv(1,3) = det_B_inv*(B(1,2)*B(2,3)-B(2,2)*B(1,3))
      B_inv(2,1) = det_B_inv*(B(3,1)*B(2,3)-B(2,1)*B(3,3))
      B_inv(2,2) = det_B_inv*(B(1,1)*B(3,3)-B(3,1)*B(1,3))
      B_inv(2,3) = det_B_inv*(B(2,1)*B(1,3)-B(1,1)*B(2,3))
      B_inv(3,1) = det_B_inv*(B(2,1)*B(3,2)-B(3,1)*B(2,2))
      B_inv(3,2) = det_B_inv*(B(3,1)*B(1,2)-B(1,1)*B(3,2))
      B_inv(3,3) = det_B_inv*(B(1,1)*B(2,2)-B(2,1)*B(1,2))
      HIN_BIG_FAI=B_inv
      WRITE(79,*)'***** HIN_BIG_FAI*****'
      WRITE(79,'(9(2X, ES13.6))')HIN_BIG_FAI
      DO K1 = 1, 3
          TEMP31(K1,1) = ALPHA(K1)
      END DO  
      TEMP1=MATMUL(TRANSPOSE(A),MATMUL(HIN_BIG_FAI,
     1MATMUL(BIG_PI, A) + BIG_GAMMA)) 
     2+ MATMUL(TRANSPOSE(A), MATMUL(D2, A) - MATMUL(D1,GAMMA * TEMP31))
     3- MATMUL(TRANSPOSE(A), MATMUL(D3, MATMUL(HIN_BIG_FAI,
     4MATMUL(BIG_PI, A) + BIG_GAMMA)))
      DOT_LAMDA = (R_F - Q_WAVE) / TEMP1(1,1)
C CALCULATE DOT_STRESS(3,1) USE(25)  
      DOT_STRESS = - MATMUL(HIN_BIG_FAI , R_HAT)
     1- DOT_LAMDA * MATMUL(HIN_BIG_FAI , MATMUL(BIG_PI , A))
     2- DOT_LAMDA * MATMUL(HIN_BIG_FAI , BIG_GAMMA)

      WRITE(79,*)'***** BIG_PI *****'
      WRITE(79,'(9(2X, ES13.6))')BIG_PI
      WRITE(79,*)'***** DOT_LAMDA*****'
      WRITE(79,'(9(2X, ES13.6))')DOT_LAMDA
      WRITE(79,*)'***** DOT_STRESS*****'
      WRITE(79,'(9(2X, ES13.6))')DOT_STRESS


C CALCULATE DOT_ALPHA(3,1) USE(24)
      DOT_ALPHA = - MATMUL(D1 , R_ALPHA)
     1+ DOT_LAMDA * MATMUL(D2 , A)
     2+ MATMUL(D3 , DOT_STRESS)
     3- GAMMA * DOT_LAMDA * MATMUL(D1 , TEMP31)

C UPDATE STRESS, ALPHA AND LAMDA USE(27)
C UPDATE STRESS, UPDATE ALPHA
      DO K1=1, 3
          STRESS(K1) = STRESS(K1) + DOT_STRESS(K1,1)
          ALPHA(K1) = ALPHA(K1) + DOT_ALPHA(K1,1)
      END DO
C UPDATE DELT_LAMDA  
      DELT_LAMDA = DELT_LAMDA + DOT_LAMDA 
C UPDATE EQEPLAS
      EQEPLAS = EQEPLAS + DOT_LAMDA

      write(79,*)'***** STRESS_NEW *****'
      WRITE(79,'(3(2X, ES13.6))')
     1(STRESS(K1), K1 = 1, 3)
      WRITE(79,*)'***** ALPHA_NEW *****'
      WRITE(79,'(3(2X, ES13.6))')
     1(ALPHA(K1), K1 = 1, 3)
      WRITE(79,*)'***** DELT_LAMDA *****'
      WRITE(79,'(2X, ES13.6)')
     1DELT_LAMDA


C ----------------------------------------------------------------      
C CALCULATE F AGAIN USE (1)  --1*1
      DO K1=1,3
          TEMP31(K1,1)=STRESS(K1)-ALPHA(K1)
      ENDDO
      TEMP1=MATMUL(MATMUL(TRANSPOSE(TEMP31),ANISO),TEMP31)
      F=SQRT(ONE/TWO*TEMP1(1,1)) 


      WRITE(79,*)'***** NEW F *****'
      write(79,'(2X, ES13.6)')F

      
      IF(TICK .GE.100)THEN
          EXIT
      ENDIF

      END DO
C -------------------------END CYCLE------------------------------      
C ----------------------------------------------------------------        

      write(*,*)'***6***'
      IF (L.EQ. 1) THEN
C CALCULATE D_EP = DDSDDE USE (37)
C CALCULATE BIG_SIGMA(3,3) USE (36)

              BIG_SIGMA = ONES + DELT_LAMMA * MATMUL(DDSDDE , DA)
     1        - DELT_LAMMA * MATMUL(DDSDDE , MATMUL(DA , D3))
C CALCULATE BIG_R(3,3) USE (36)
                    B=BIG_SIGMA
            det_B = B(1,1)*(B(2,2)*B(3,3) - B(3,2)*B(2,3)) 
     1-        B(2,1)*(B(1,2)*B(3,3) - B(3,2)*B(1,3)) 
     2+        B(3,1)*(B(1,2)*B(2,3) - B(2,2)*B(1,3))
      det_B_inv = 1.D0/det_B
      B_inv(1,1) = det_B_inv*(B(2,2)*B(3,3)-B(3,2)*B(2,3))
      B_inv(1,2) = det_B_inv*(B(3,2)*B(1,3)-B(1,2)*B(3,3))
      B_inv(1,3) = det_B_inv*(B(1,2)*B(2,3)-B(2,2)*B(1,3))
      B_inv(2,1) = det_B_inv*(B(3,1)*B(2,3)-B(2,1)*B(3,3))
      B_inv(2,2) = det_B_inv*(B(1,1)*B(3,3)-B(3,1)*B(1,3))
      B_inv(2,3) = det_B_inv*(B(2,1)*B(1,3)-B(1,1)*B(2,3))
      B_inv(3,1) = det_B_inv*(B(2,1)*B(3,2)-B(3,1)*B(2,2))
      B_inv(3,2) = det_B_inv*(B(3,1)*B(1,2)-B(1,1)*B(3,2))
      B_inv(3,3) = det_B_inv*(B(1,1)*B(2,2)-B(2,1)*B(1,2))
      HIN_BIG_SIGMA=B_inv
              
              BIG_R = MATMUL(HIN_BIG_SIGMA , DDSDDE)
              
C CALCULATE BIG_LAMDA(3,1) USE (36)
              DO K1 = 1, 3
                  TEMP31(K1,1) = ALPHA(K1)
              END DO
              BIG_LAMDA = DELT_LAMMA * MATMUL(HIN_BIG_SIGMA, MATMUL(DDSD
     1        DE,MATMUL(DA, MATMUL(D2,A))))
     2        - DELT_LAMMA * GAMMA * MATMUL(HIN_BIG_SIGMA, MATMUL(DDSDDE
     3        ,MATMUL(DA, MATMUL(D1,TEMP31))))
C CALCULATE BIG_KESI(1,1) USE (37)
              BIG_KESI = MATMUL(TRANSPOSE(A), MATMUL(BIG_R,A))
     1        - MATMUL(TRANSPOSE(A), BIG_LAMDA)
     2        + MATMUL(TRANSPOSE(A), MATMUL(D2,A))
     3        - GAMMA * MATMUL(TRANSPOSE(A), MATMUL(D1,TEMP31))

C CALCULATE DDSDDE
              TEMP33 = BIG_KESI(1,1) * ONES + MATMUL(BIG_LAMDA, 
     1        MATMUL(TRANSPOSE(A), D3)) 
     2        - MATMUL(BIG_R, MATMUL(A, MATMUL(TRANSPOSE(A), D3)))
C !!!!!SHOULD CALCULATE HIN_TEMP33 FIRST
                    B=TEMP33
            det_B = B(1,1)*(B(2,2)*B(3,3) - B(3,2)*B(2,3)) 
     1-        B(2,1)*(B(1,2)*B(3,3) - B(3,2)*B(1,3)) 
     2+        B(3,1)*(B(1,2)*B(2,3) - B(2,2)*B(1,3))
      det_B_inv = 1.D0/det_B
      B_inv(1,1) = det_B_inv*(B(2,2)*B(3,3)-B(3,2)*B(2,3))
      B_inv(1,2) = det_B_inv*(B(3,2)*B(1,3)-B(1,2)*B(3,3))
      B_inv(1,3) = det_B_inv*(B(1,2)*B(2,3)-B(2,2)*B(1,3))
      B_inv(2,1) = det_B_inv*(B(3,1)*B(2,3)-B(2,1)*B(3,3))
      B_inv(2,2) = det_B_inv*(B(1,1)*B(3,3)-B(3,1)*B(1,3))
      B_inv(2,3) = det_B_inv*(B(2,1)*B(1,3)-B(1,1)*B(2,3))
      B_inv(3,1) = det_B_inv*(B(2,1)*B(3,2)-B(3,1)*B(2,2))
      B_inv(3,2) = det_B_inv*(B(3,1)*B(1,2)-B(1,1)*B(3,2))
      B_inv(3,3) = det_B_inv*(B(1,1)*B(2,2)-B(2,1)*B(1,2))
      HIN_TEMP33=B_inv
              
              DDSDDE = MATMUL(HIN_TEMP33, BIG_R * BIG_KESI(1,1)
     1        + MATMUL(BIG_LAMDA, MATMUL(TRANSPOSE(A), BIG_R))
     2        - MATMUL(BIG_R, MATMUL(A, MATMUL(TRANSPOSE(A), BIG_R))))

                    WRITE(79,*)'**ddsdde_NEW**'
      WRITE(79,'(9(2X, ES13.6))')DDSDDE
              
      ENDIF

C      DO K1=1,3
C          DO K2 = 1, 3
C              STRESS(K2) = STRESS(K2) + DDSDDE(K2, K1) * DSTRAN(K1)
C          END DO
C      END DO 




C STORE ELASTIC STRAINS, PLASTIC STRAINS, BACK STRESSES, EQUIVALENT PLASTIC STRAIN
C UPDATE EPLAS
      DELT_EPLAS= DELT_LAMDA * A
      DO K1=1,3
          EPLAS(K1)=EPLAS(K1)+DELT_EPLAS(K1,1)
      ENDDO
C UPDATE EELAS
      DO K1=1,3
          EELAS(K1)=EELAS(K1)-DELT_EPLAS(K1,1)
      ENDDO
      
C IN STATE VARIABLE ARRAY
      WRITE(79,*)'**EELAS_NEW**'
      WRITE(79,'(3(2X, ES13.6))')
     1(EELAS(K1), K1 = 1, 3)
      WRITE(79,*)'**EPLAS_NEW**'
      WRITE(79,'(3(2X, ES13.6))')
     1(EPLAS(K1), K1 = 1, 3)
      WRITE(79,*)'**ALPHA_NEW**'
      WRITE(79,'(3(2X, ES13.6))')
     1(ALPHA(K1), K1 = 1, 3)
          DO K1=1,3
              STATEV(K1)=EELAS(K1)
              STATEV(K1+3)=EPLAS(K1)
              STATEV(K1+6)=ALPHA(K1)   
          END DO
C UPDATE EQUIVALENT PLASTIC STRAIN
      STATEV(10)=EQEPLAS
C IF WE WANT TO CALCULATE PLASTIC STRAIN ENERGY
      DO K1=1,3
          SPD=SPD+(STRESS(K1)+OLDS(K1))*(EPLAS(K1)-OLDPL(K1))/TWO
      END DO
      STATEV(11)=SPD

      WRITE(79,*)'**STATEV_NEW**'
      WRITE(79,'(11(2X, ES13.6))')
     1(STATEV(K1), K1 = 1, 11)


C      WRITE(79,*)ONES
      WRITE(79,*)'**STRESS_NEW**'
      WRITE(79,'(3(2X, ES13.6))')
     1(STRESS(K1), K1 = 1, 3)



      RETURN
      END





!      subroutine MATRIX_INV(A,A_inv,ONES,CMNAME)
!C      !
!C      ! Returns A_inv, the inverse and det_A, the determinant
!C      ! Note that the det is of the original matrix, not the
!C      ! inverse
!C
!
!      INCLUDE 'ABA_PARAM.INC'
!      CHARACTER*80 CMNAME
!      DIMENSION  A(3,3),A_inv(3,3),ONES(3,3)
!
!      istat = 1
!      
!            !if (det_A .le. 0.d0) then
!      !    write(*,*) 'WARNING: subroutine MATRIX_INV:'
!      !    write(*,*) 'WARNING: det of mat=',det_A
!      !    istat = 0
!      !    return
!      !end if
!
!      det_B = B(1,1)*(B(2,2)*B(3,3) - B(3,2)*B(2,3)) 
!     1-        B(2,1)*(B(1,2)*B(3,3) - B(3,2)*B(1,3)) 
!     2+        B(3,1)*(B(1,2)*B(2,3) - B(2,2)*B(1,3))
!      det_B_inv = 1.D0/det_B
!      B_inv(1,1) = det_B_inv*(B(2,2)*B(3,3)-B(3,2)*B(2,3))
!      B_inv(1,2) = det_B_inv*(B(3,2)*B(1,3)-B(1,2)*B(3,3))
!      B_inv(1,3) = det_B_inv*(B(1,2)*B(2,3)-B(2,2)*B(1,3))
!      B_inv(2,1) = det_B_inv*(B(3,1)*B(2,3)-B(2,1)*B(3,3))
!      B_inv(2,2) = det_B_inv*(B(1,1)*B(3,3)-B(3,1)*B(1,3))
!      B_inv(2,3) = det_B_inv*(B(2,1)*B(1,3)-B(1,1)*B(2,3))
!      B_inv(3,1) = det_B_inv*(B(2,1)*B(3,2)-B(3,1)*B(2,2))
!      B_inv(3,2) = det_B_inv*(B(3,1)*B(1,2)-B(1,1)*B(3,2))
!      B_inv(3,3) = det_B_inv*(B(1,1)*B(2,2)-B(2,1)*B(1,2))
!
!      
!      ONES(3,3)=888888888
!
!      !OPEN(unit=78,file='F:\JDGwork\ABAQUS\jdgtest.txt',status='unknown')
!      !write(78,*)'>>>>>>>>>>>>>>>>>>>>>>>>>> inverse>>>>>>>>>>>>>>>>>>'
!      !WRITE(78,'(2(2X, ES13.6))')det_A,det_A_inv
!      !WRITE(78,'(9(2X, ES13.6))')A
!      !WRITE(78,'(9(2X, ES13.6))')A_inv
!      !WRITE(78,'(9(2X, ES13.6))')ONES
!
!
!
!      return
!      end subroutine MATRIX_INV