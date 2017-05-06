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

%% RESULTADOS
echo off
disp('   ');
disp('Matriz de admitancias nodales YBUS');
disp(num2str(YBUS,'%8.3f'))