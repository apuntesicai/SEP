% TRABAJO DE SEP TEMA4: FLUJO DE CARGAS
% Fecha: 27/04/2016
% Autores: IGNACIO SANZ, GUILLERMO MAIRAL, JAVIER MIÑARRO,
% PEDRO PARDO, IGNACIO PASTOR, SANTIAGO LÓPEZ

clear all
format short
disp('  ')
disp('  ')
echo on

% EC 1
% Los datos de una red de 4 nudos son los siguientes:

% INI	FIN	Descripción
% 2	3	TRANSFORMADOR: 200 MVA, 410/20kV; ucc = 13%
% 1	4	LÍNEA: 400kV; zs = 0.03+j0.31 p.u., yp = j0.87 p.u.
% 1	4	LÍNEA: 400kV; zs = 0.03+j0.31 p.u., yp = j0.87 p.u.
% 1	2	LÍNEA: 400kV; zs = 0.02+j0.23 p.u., yp = j0.52 p.u.
% 2	4	LÍNEA: 400kV; zs = 0.01+j0.19 p.u., yp = j0.41 p.u.

% Los datos de los transformadores son relativos a sus propias bases.
% Los datos de las líneas son relativos a las bases del sistema.
% Obtener la matriz de admitancias nodales Ybus, considerando bases trifásicas del sistema de 400/220kV y 100 MVA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTES
J = sqrt(-1) ;

% Bases TRIFASICAS
Sbase = 100/3 ;
Ubase = [ 400 20 ]/sqrt(3) ;
Ibase = Sbase./Ubase 
Zbase = Ubase./Ibase 

% YBUS
YBUS = zeros(4,4)

%% Líneas
% Línea 1 - 4
zs_14 = 0.03 + J*0.31;
yp_14 = J*0.87;
YBUS(1,1) = YBUS(1,1) + 1/zs_14 + yp_14/2 ;
YBUS(1,4) = YBUS(1,4) - 1/zs_14 ; 
YBUS(4,1) = YBUS(4,1) - 1/zs_14 ; 
YBUS(4,4) = YBUS(4,4) + 1/zs_14 + yp_14/2 ; 
YBUS

% Línea 1 - 4 BIS
zs_14bis = 0.03 + J*0.31;
yp_14bis = J*0.87;
YBUS(1,1) = YBUS(1,1) + 1/zs_14 + yp_14/2 ;
YBUS(1,4) = YBUS(1,4) - 1/zs_14bis ; 
YBUS(4,1) = YBUS(4,1) - 1/zs_14bis ; 
YBUS(4,4) = YBUS(4,4) + 1/zs_14bis + yp_14bis/2 ; 
YBUS

% Línea 1 - 2
zs_12 = 0.02 + J*0.23;
yp_12 = J*0.52;
YBUS(1,1) = YBUS(1,1) + 1/zs_12 + yp_12/2 ;
YBUS(1,2) = YBUS(1,2) - 1/zs_12 ; 
YBUS(2,1) = YBUS(2,1) - 1/zs_12 ; 
YBUS(2,2) = YBUS(2,2) + 1/zs_12 + yp_12/2 ; 
YBUS

% Línea 2 - 4
zs_24 = 0.01 + J*0.19;
yp_24 = J*0.41;
YBUS(2,2) = YBUS(2,2) + 1/zs_24 + yp_24/2 ;
YBUS(2,4) = YBUS(2,4) - 1/zs_24 ; 
YBUS(4,2) = YBUS(4,2) - 1/zs_24 ; 
YBUS(4,4) = YBUS(4,4) + 1/zs_24 + yp_24/2 ; 
YBUS

%% Transformadores
% Transformador 2 - 3
trafos_Snom23 = 200/3 ;
trafos_VnomA23 = 410/sqrt(3) ;
trafos_VnomB23 = 20/sqrt(3) ;
trafos_ucc23 = 13/100 ;
trafos_zcc23_basesTRAFO = J*trafos_ucc23 ; 
trafos_zcc23 = trafos_zcc23_basesTRAFO*((trafos_VnomB23^2)/trafos_Snom23)/Zbase(2) ; % La ponemos en el lado de BAJA
trafos_t23 = (trafos_VnomA23/trafos_VnomB23)/(Ubase(1)/Ubase(2)) ;  % en el lado de ALTA
YBUS(2,2) = YBUS(2,2) + 1/((trafos_t23*trafos_t23)*trafos_zcc23) ;
YBUS(2,3) = YBUS(2,3) - 1/((trafos_t23)*trafos_zcc23) ; 
YBUS(3,2) = YBUS(3,2) - 1/((trafos_t23)*trafos_zcc23) ; 
YBUS(3,3) = YBUS(3,3) + 1/trafos_zcc23  ;  
YBUS

%%  Vectores de tensión del apartado anterior
V1 = 1.1005*exp(J*(-11.264)*pi/180);
V2 = 1.0737*exp(J*(-2.193)*pi/180);
V3 = 1.0300*exp(J*(-0.121)*pi/180);
V4 = 1.0100*exp(J*0.00*pi/180);

% Vector de tensiones (vector columna)
vecV = [V1; V2; V3; V4]

% Vector de corrientes netas inyectadas
vecI = YBUS*vecV

%% Potencias netas inyectadas, esto es, potencias calculadas
% Sc=Pc+J*Qc
Sc1 = vecV(1)*conj(vecI(1));
Sc2 = vecV(2)*conj(vecI(2));
Sc3 = vecV(3)*conj(vecI(3));
Sc4 = vecV(4)*conj(vecI(4));

% Dado que para las potencias que nos piden no tenemos mismatches:
% La potencia activa generada por el nudo 4 (SW) es:
Pc4 = real(Sc4);

% La potencia reactiva generada por el nudo 4 (SW) es:
Qc4 = imag(Sc4);

% La potencia reactiva generada por el nudo 3 es:
Qc3 = imag(Sc3);

%% Balance de potencias en lineas y trafos
% Línea 1 - 4
Ilinea_14 = V1*(yp_14/2) + (V1-V4)/(zs_14);
Ilinea_41 = V4*(yp_14/2) + (V4-V1)/(zs_14);
Plinea_14 = real(V1*conj(Ilinea_14));
Plinea_41 = real(V4*conj(Ilinea_41));

% Línea 1 - 4 bis
Ilinea_14bis = V1*(yp_14bis/2) + (V1-V4)/(zs_14bis);
Ilinea_41bis = V4*(yp_14bis/2) + (V4-V1)/(zs_14bis);
Plinea_14bis = real(V1*conj(Ilinea_14bis));
Plinea_41bis = real(V4*conj(Ilinea_41bis));

% Línea 1 - 2
Ilinea_12 = V1*(yp_12/2) + (V1-V2)/(zs_12);
Ilinea_21 = V2*(yp_12/2) + (V2-V1)/(zs_12);
Plinea_12 = real(V1*conj(Ilinea_12));
Plinea_21 = real(V2*conj(Ilinea_21));

% Línea 2 - 4
Ilinea_24 = V2*(yp_24/2) + (V2-V4)/(zs_24);
Ilinea_42 = V4*(yp_24/2) + (V4-V2)/(zs_24);
Plinea_24 = real(V2*conj(Ilinea_24));
Plinea_42 = real(V4*conj(Ilinea_42));

% Transformador 2 - 3
V2prima = V2/trafos_t23;
Itrafo_2prima3 = (V2prima-V3)/(trafos_zcc23);
Itrafo_23 = Itrafo_2prima3/trafos_t23;
Itrafo_32 = (V3-V2prima)/(trafos_zcc23);
Ptrafo_23 = real(V2*conj(Itrafo_23));
Ptrafo_32 = real(V3*conj(Itrafo_32));

%% RESULTADOS
echo off
disp('   ');
disp('Potencia activa generada por el nudo 4')
disp([ '    ' num2str(Pc4) ] )
disp('   ');
disp('Potencia reactiva generada por el nudo 4')
disp([ '    ' num2str(Qc4) ] )
disp('   ');
disp('Potencia reactiva generada por el nudo 3')
disp([ '    ' num2str(Qc3) ] )
disp('   ');
disp('Balance de potencia activa en lineas')
disp([ '· 1 a 4 : ' num2str(Plinea_14,'%7.3f') ' ; 4 a 1 : ' num2str(Plinea_41,'%7.3f') ' ; perdidas : ' num2str(Plinea_14+Plinea_41,'%7.3f') ] )
disp([ '· 1 a 4(bis) : ' num2str(Plinea_14bis,'%7.3f') ' ; 4 a 1(bis) : ' num2str(Plinea_41bis,'%7.3f') ' ; perdidas : ' num2str(Plinea_14bis+Plinea_41bis,'%7.3f') ] )
disp([ '· 1 a 2 : ' num2str(Plinea_12,'%7.3f') ' ; 2 a 1 : ' num2str(Plinea_21,'%7.3f') ' ; perdidas : ' num2str(Plinea_12+Plinea_21,'%7.3f') ] )
disp([ '· 2 a 4 : ' num2str(Plinea_24,'%7.3f') ' ; 4 a 2 : ' num2str(Plinea_42,'%7.3f') ' ; perdidas : ' num2str(Plinea_24+Plinea_42,'%7.3f') ] )
disp('Balance de potencia activa en trafos')
disp([ '· 2 a 3 : ' num2str(Ptrafo_23,'%7.3f') ' ; 3 a 2 : ' num2str(Ptrafo_32,'%7.3f') ' ; perdidas : ' num2str(Ptrafo_23+Ptrafo_32,'%7.3f') ] )




