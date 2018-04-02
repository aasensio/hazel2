c rutina sublande2
c       code1 representa mediante letras minusculas los momentos angulares
c       orbitales semienteros, asi p=1/2, f=3/2, h=5/2,k=7/2,m=9/2,o=11/2
c                              r=13/2,t=15/2, u=17/2,v=19/2,w=21/2
      SUBROUTINE sublande2(MC,MULT,DESIGN,TAM,JI,JF,DL0,NP,NL,NR,DLP,DLL,
     *DLR,SP,SL,SR,g)

	integer mult(*),ji(*),jf(*),mc,np,nl,nr
	real tam(*),g(*)
	real dlp(*),dll(*),dlr(*),sp(*),sl(*),sr(*)

      CHARACTER*1 DESIGN(2),CODE(13),code1(21)
      data CODE/'S','P','D','F','G','H','I','K','L','M','N',
     *'O','Q'/
        data code1/'p','1','f','2','h','3','k','4','m','5','o'
     & ,'6','r','7','t','8','u','9','v','0','w'/


      JI1=JI(1)
      JI2=JI(2)
      JF1=JF(1) 
      JF2=JF(2) 
      DO 58 I=1,MC  
      DLL(I)=0. 
      SL(I)=0.  
      DLR(I)=0. 
  58  SR(I)=0.
      G(1)=0.   
      G(2)=0.   
      LEVEL=0   
      IF((JI1+JF1).EQ.0) LEVEL=1
      IF((JI2+JF2).EQ.0) LEVEL=-1   
      IF(LEVEL) 3,4,5   
3     I=1                    ! TRIPLETS WITH G=0
      GO TO 6   
5     I=2   
6     TAM(I)=1.              ! TOTAL ANGULAR MOMENTUM

      SPIN=0.5*FLOAT(MULT(I)-1)
      DO 7 J=1,13   
      IF(DESIGN(I).NE.CODE(J)) GO TO 7  
      OAM=J-1                ! ORBITAL ANGULAR MOMENTUM 
      GO TO 8   
7     CONTINUE  

      do 70 j=1,21

      if(design(i).ne.code1(j)) goto 70
      oam=float(j)/2.                ! ORBITAL ANGULAR MOMENTUM (semientero)
      go to 8 
70    continue

      print*,' '
      print*,'STOP: Check the transitions in the file containing the atomic parameters.'
      print*,'Some letter used to code the orbital angular momentum is erroneous.'      
      print*,' '
      print*,'__________________________________________________________________________________'
      stop


8     G(I)=1.5+(SPIN*(1.+SPIN)-OAM*(1.+OAM))/4. 
13    NP=1                   ! TRIPLETS WITH G=0 OR G1=G2   
      NL=NP
      NR=NP 
      SP(1)=1.  
      SL(1)=1.  
      SR(1)=1.  
      DLP(1)=0. 
      DLL(1)=G(I)   
      DLR(1)=-DLL(1)
      GO TO 9   
4     DO 10 I=1,2
      SPIN=0.5*FLOAT(MULT(I)-1) 
      DO 11 J=1,13  
      IF(DESIGN(I).NE.CODE(J)) GO TO 11 
      OAM=J-1   
      GO TO 10  
11    CONTINUE 

      do 110 j=1,21
      if(design(I).ne.code1(j)) goto 110  
      oam=float(j)/2.                ! ORBITAL ANGULAR MOMENTUM (semientero)
      go to 10
110   continue
 
      print*,' '
      print*,'STOP: Check the transitions in the file containing the atomic parameters.'
      print*,'Some letter used to code the orbital angular momentum is erroneous.'      
      print*,' '
      print*,'__________________________________________________________________________________'
      stop





C	LANDE FACTOR FOR EACH LABEL
10    G(I)=1.5+(SPIN*(1.+SPIN)-OAM*(1.+OAM))/(2.*TAM(I)*(1.+TAM(I)))
      IF(ABS(G(2)-G(1)).GT.5.E-6) GO TO 12  
      I=2
      GO TO 13  
12    LEVEL=JI2-JI1 
      IF(JF1.EQ.5) GO TO 14 
      IF(LEVEL) 15,16,17     ! INTEGRAL J'S 
15    NP=2*JI2+1
19    NL=NP 
      NR=NP 
      IF(NP.LE.MC) GO TO 18 
      STOP 'EXIT ZEEMAN'
16    NP=2*JI2  
      GO TO 19  
17    NP=2*JI2-1
      GO TO 19  
18    MUMIN=-JI2
      IF(JI1.LT.JI2) MUMIN=-JI1 
      MUMAX=-MUMIN  
      I=0   
      DO 20 MU=MUMIN,MUMAX   ! PI COMPONENTS
      IF(MU.EQ.0.AND.LEVEL.EQ.0) GO TO 20   
      I=I+1
      DLP(I)=FLOAT(MU)*(G(1)-G(2))  
      J=MU**2
      IF(LEVEL) 21,22,23
21    SP(I)=2*(JI1**2-J)
      GO TO 20  
22    SP(I)=2*J 
      GO TO 20  
23    SP(I)=2*(JI2**2-J)
20    CONTINUE  
      MUMIN=1-JI1   
      MUMAX=JI2 
      I=0   
      DO 24 MU=MUMIN,MUMAX   ! R-SIGMA COMPONENTS   
      I=I+1 
      DLR(I)=FLOAT(MU)*(G(1)-G(2))-G(1) 
      IF(LEVEL) 25,26,27
25    SR(I)=(JI1-MU)*(JI2-MU+2) 
      GO TO 24  
26    SR(I)=(JI2+MU)*(JI2-MU+1) 
      GO TO 24
27    SR(I)=(JI2+MU)*(JI1+MU)   
24    CONTINUE  
      MUMIN=-JI2
      MUMAX=JI1-1   
      I=0   
      DO 28 MU=MUMIN,MUMAX   ! L-SIGMA COMPONENTS   
      I=I+1 
      DLL(I)=FLOAT(MU)*(G(1)-G(2))+G(1) 
      IF(LEVEL) 29,30,31
29    SL(I)=(JI1+MU)*(JI2+MU+2) 
      GO TO 28  
30    SL(I)=(JI2-MU)*(JI2+MU+1) 
      GO TO 28  
31    SL(I)=(JI2-MU)*(JI1-MU)   
28    CONTINUE  
      GO TO 57  
14    I=2*JI1+1              ! HALF-INTEGRAL J'S

      IF(LEVEL) 32,33,34
32    NP=I-1
      NL=NP
      NR=NP 
36    IF(NP.LE.MC) GO TO 35 
      STOP 'EXIT ZEEMAN'
33    NP=I+1
      NL=I  
      NR=I  
      GO TO 36  
34    NP=I+1
      NL=NP 
      NR=NP 
      GO TO 36  
35    MUMIN=-JI2

      
      IF(JI1.LT.JI2) MUMIN=-JI1 
      MUMAX=1-MUMIN 
      I=0   
      DO 37 MU=MUMIN,MUMAX   ! PI COMPONENTS
      I=I+1 
      SPIN=FLOAT(MU)-0.5
      DLP(I)=(G(1)-G(2))*SPIN   
      SPIN=SPIN**2
      IF(LEVEL) 38,39,40
38    SP(I)=2.*((FLOAT(JI1)+0.5)**2-SPIN)   
      GO TO 37  
39    SP(I)=2.*SPIN 
      GO TO 37  
40    SP(I)=2.*((FLOAT(JI2)+0.5)**2-SPIN)   
37    CONTINUE  
      MUMIN=-JI1
      MUMAX=JI2 
      I=0   
      DO 41 MU=MUMIN,MUMAX   ! R-SIGMA COMPONENTS   
      I=I+1 
      DLR(I)=(FLOAT(MU)+0.5)*(G(1)-G(2))-G(1)   
      IF(LEVEL) 42,43,44
42    SR(I)=(JI1-MU)*(JI2-MU+2) 
      GO TO 41  
43    SR(I)=(JI2+MU+1)*(JI2-MU+1)   
      GO TO 41  
44    SR(I)=(JI2+MU+1)*(JI2+MU) 
41    CONTINUE
      MUMIN=-MUMAX  
      MUMAX=JI1 
      I=0   
      DO 45 MU=MUMIN,MUMAX   ! L-SIGMA COMPONENTS   
      I=I+1 
      DLL(I)=(FLOAT(MU)-0.5)*(G(1)-G(2))+G(1)   
      IF(LEVEL) 46,47,48
46    SL(I)=(JI1+MU)*(JI2+MU+2) 
      GO TO 45  
47    SL(I)=(JI2-MU+1)*(JI2+MU+1)   
      GO TO 45  
48    SL(I)=(JI2-MU+1)*(JI1-MU+1)
45    CONTINUE
57    SUM=0.
      DO 49 I=1,NP
49    SUM=SUM+SP(I)
      DO 50 I=1,NP
50    SP(I)=SP(I)/SUM
      SPIN=0.
      SUM=0.
      DO 51 I=1,NL           ! NL=NR ASSUMED.
      SPIN=SPIN+SL(I)
51    SUM=SUM+SR(I)
      DO 52 I=1,NL
      SL(I)=SL(I)/SPIN
52    SR(I)=SR(I)/SUM
   9  CONTINUE
      IF(MC.LE.18) GO TO 59 ! SKIP PRINTING IF SYNTHESIS
      WRITE(6,53) (MULT(I),DESIGN(I),JI(I),JF(I),I=1,2)
  53  FORMAT(3X,'Transition:',I3,A1,I2,'.',I1,' ----',I3,A1,I2,'.',I1/
     *3X,'Zeeman Pattern (Lorentz units) and normalized Intensities:'/
     *11X,'Pi',17X,'Sigma+',15X,'Sigma-')
 59   DO 54 I=1,NP
      IF(MC.LE.18) GO TO 60
      WRITE(6,55) DLP(I),SP(I),DLL(I),SL(I),DLR(I),SR(I)
  55  FORMAT(2X,3(F10.6,F11.8))
  60  DLP(I)=DLP(I)*DL0
      DLL(I)=DLL(I)*DL0
54    DLR(I)=DLR(I)*DL0
      IF(MC.LE.18) RETURN
      WRITE(6,56) G(1),G(2)
  56  FORMAT(3X,'Lande Factors: g(lower)=',F11.8,', g(upper)=',F11.8)
      RETURN
      END   
c________________________________________________________________________
