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

Ybus= YBUS;

%% Flujo de Cargas: Solución Completa.
% Se incluye también a su derecha la matriz de admitancias nodales Ybus de la misma,
% considerando como potencia base 100 MVA
% Tomando como punto de partida en perfil plano de tensiones,
% obtener la solución del problema de flujo de cargas planteado para esta red (considerar tolerancia en p.u. = 10-2),
% empleando para ello el flujo de cargas completo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATOS
x12 = 0.23 ;
x14bis = 0.31;
x14 = 0.31 ;
x23 = 0.065 ;
x24 = 0.19;

% Ybus DC
YbusDC = [
	[ (1/x12+1/x14+1/x14bis)       -1/x12               0        -(1/x14+1/x14bis)] ;
	[ -1/x12                 (1/x12+1/x23+1/x24)      -1/x23          -1/x24      ] ;
	[   0                    -1/x23          (1/x23)         0        ] ;
    [ -(1/x14+1/x14bis)            -1/x24                0      (1/x24+1/x14+1/x14bis)];                          
] 

% TIPO de NUDOS:
%   Nudo   Pe    Qe      tipo    mismP   mismQ   TH esp.   V esp.
%      1   NO    NO  -->   SW       NO      NO       SI       SI
%      2   SI    NO  -->   PV       SI      NO       NO       SI
%      3   SI    SI  -->   PQ       SI      SI       NO       NO

% El slack es el 1 -> recortar fila y columna
YbusDC_mismP = YbusDC(1:3,1:3)
invYbusDC_mismP = YbusDC_mismP\eye(size(YbusDC_mismP)) % A\B = inv(A)*B

% Potencias especificadas
Pe1 =  -210/Sbase %   Activa nudo 2, en pu
Pe2 =     0/Sbase
Pe3 =    60/Sbase %   Activa nudo 3, en pu

% Vector potencias inyectadas
Piny = [
	Pe1 ;
    Pe2 ;
	Pe3 ;
	]

% Ángulos de las tensiones
vecTH = invYbusDC_mismP*Piny
vecTH = [  vecTH ; 0.0 ] % Añadir el del slack
TH1 = vecTH(1)  
TH2 = vecTH(2) 
TH3 = vecTH(3) 
TH4 = vecTH(4) % referencia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTES
J = sqrt(-1) ;
toleranciaMISMATCHES = 1.0e-2 ;
numITERACIONES = 0 ;
Sbase = 100 ;
numeroNUDOS = 4 ;

% Ybus
Gbus = real(Ybus)
Bbus = imag(Ybus)

% TIPO de NUDOS:
%   Nudo   Pe    Qe      tipo    mismP   mismQ   TH esp.   V esp.
%      1   SI    SI  -->   PQ       SI      SI       NO       NO
%      2   SI    SI  -->   PQ       SI      SI       NO       NO
%      3   SI    SI  -->   PV       SI      NO       NO       SI
%      4   NO    NO  -->   SW       NO      NO       NO       NO

% Potencias especificadas
Pe1 =  -210/Sbase %   Activa nudo 1, en pu
Qe1 =   -20/Sbase %   Reactiva nudo 1, en pu
Pe2 =     0/Sbase %   Activa nudo 2, en pu
Qe2 =     0/Sbase %   Reactiva nudo 2, en pu
Pe3 =    60/Sbase %   Reactiva nudo 3, en pu

% Tensiones en modulo (V) argumento (TH) y forma fasorial (fasV)
% ??? -> Perfil plano de tensiones
V1 = 1.00 ;   fasV1 = V1*exp(J*TH1) ;
V2 = 1.00 ;   fasV2 = V2*exp(J*TH2) ;
V3 = 1.03 ;   fasV3 = V3*exp(J*TH3) ;
V4 = 1.01 ;   fasV4 = V3*exp(J*TH4) ;

% Vector de tensiones en forma fasorial (OJO: vector columna)
vecV = [ fasV1 ; fasV2 ; fasV3 ;fasV4]  

% Guardar para RESULTADOS
X_inicial = [ TH1 ; TH2 ; TH3 ; V1 ; V2 ]  ;

% COMPROBACION MISMATCHES
% Vector de corrientes netas inyectadas
vecI = Ybus*vecV

% Potencias netas inyectadas, esto es, potencias calculadas (Sc=Pc+J*Qc)
Sc1 = vecV(1)*conj(vecI(1))
Sc2 = vecV(2)*conj(vecI(2))
Sc3 = vecV(3)*conj(vecI(3))
Sc4 = vecV(4)*conj(vecI(4))

% En forma de vector
Sc = [ Sc1 ; Sc2 ; Sc3 ; Sc4] ;

% Mismatches, segun tabla:
MP1 = Pe1 - real(Sc1)
MQ1 = Qe1 - imag(Sc1)
MP2 = Pe2 - real(Sc2)
MQ2 = Qe2 - imag(Sc2)
MP3 = Pe3 - real(Sc3)

% Vector de Mismatches (OJO: vector columna)
vecMISMATCHES = [ MP1 ; MP2 ; MP3 ; MQ1 ; MQ2 ]  

% Guardar para RESULTADOS
vecMISMATCHES_inicial = vecMISMATCHES ;

% Error MAXIMO
errmaxMISMATCHES = max(abs(vecMISMATCHES)) 

% Proceso iteratico: se entrara y no se saldra hasta que se alcanze la tolerancia
while (errmaxMISMATCHES>toleranciaMISMATCHES)
	
	% Si estamos aqui, es necesario hacer una iteracion mas
	numITERACIONES = numITERACIONES + 1 
	
	% Matriz Jacobiana
	% IMPORTANTE: no confundir el NUMERO de nudo con POSICION en el subjacobiano correspondiente
	% subjacobiano Jpth (3x3)
	Jpth (1,1) = - imag(Sc1) - V1^2*Bbus(1,1) ;  % dPc1_dTH1
	Jpth (1,2) = V1*V2*( Gbus(1,2)*sin(TH1-TH2) - Bbus(1,2)*cos(TH1-TH2) ) ; % dPc1_dTH2
    Jpth (1,3) = V1*V3*( Gbus(1,3)*sin(TH1-TH3) - Bbus(1,3)*cos(TH1-TH3) ) ; % dPc1_dTH3
    
    Jpth (2,1) = V1*V2*( Gbus(1,2)*sin(TH1-TH2) - Bbus(1,2)*cos(TH1-TH2) ) ; % dPc1_dTH2
    Jpth (2,2) = - imag(Sc2) - V2^2*Bbus(2,2) ; % dPc2_dTH2
	Jpth (2,3) = V2*V3*( Gbus(2,3)*sin(TH2-TH3) - Bbus(2,3)*cos(TH2-TH3) ) ; % dPc2_dTH3

    Jpth (3,1) = V3*V1*( Gbus(3,1)*sin(TH3-TH1) - Bbus(3,1)*cos(TH3-TH1) ) ; % dPc3_dTH1
    Jpth (3,3) = - imag(Sc3) - V3^2*Bbus(3,3) ; % dPc3_dTH3
	Jpth (3,2) = V3*V2*( Gbus(3,2)*sin(TH3-TH2) - Bbus(3,2)*cos(TH3-TH2) ) ; % dPc3_dTH2
	Jpth
    
    % subjacobiano Jpv (3x2)
	Jpv (1,1) = real(Sc1) + V1^2*Gbus(1,1)  ; % V1*dPc1_dV1
    Jpv (1,2) = V1*V2*( Gbus(1,2)*cos(TH1-TH2) + Bbus(1,2)*sin(TH1-TH2) ) ;  % V2*dPc1_dV2
	
    Jpv (2,1) = V2*V1*( Gbus(2,1)*cos(TH2-TH1) + Bbus(2,1)*sin(TH2-TH1) ) ;  % V1*dPc2_dV1
    Jpv (2,2) = real(Sc2) + V2^2*Gbus(2,2)  ; % V2*dPc2_dV2
    
    Jpv (3,1) = V3*V1*( Gbus(3,1)*cos(TH3-TH1) + Bbus(3,1)*sin(TH3-TH1) ) ;  % V1*dPc3_dV1
    Jpv (3,2) = real(Sc3) + V3^2*Gbus(3,3)  ; % V2*dPc3_dV2
    Jpv
    
	% subjacobiano Jqth (2x3)
	Jqth (1,1) = real(Sc1) - V1^2*Gbus(1,1) ;  % dQc1_dTH1
    Jqth (1,2) = - V1*V2*( Gbus(1,2)*cos(TH1-TH2) + Bbus(1,2)*sin(TH1-TH2) ) ;  % dQc1_dTH2
	Jqth (1,3) = - V1*V3*( Gbus(1,3)*cos(TH1-TH3) + Bbus(1,3)*sin(TH1-TH3) ) ;  % dQc1_dTH3
	
    Jqth (2,2) = real(Sc2) - V2^2*Gbus(2,2) ;  % dQc2_dTH2
    Jqth (2,1) = - V2*V1*( Gbus(2,1)*cos(TH2-TH1) + Bbus(2,1)*sin(TH2-TH1) ) ;  % dQc2_dTH1
	Jqth (2,3) = - V2*V3*( Gbus(2,3)*cos(TH2-TH3) + Bbus(2,3)*sin(TH2-TH3) ) ;  % dQc2_dTH3
    Jqth
    
	% subjacobiano Jqv (2x2)
	Jqv (1,1) = imag(Sc1) - V1^2*Bbus(1,1) ; % V1*dQc1_dV1
	Jqv (1,2) = V1*V2*( Gbus(1,2)*sin(TH1-TH2) - Bbus(1,2)*cos(TH1-TH2) ) ; % V2*dQc1_dV2
    
    Jqv (2,1) = V1*V2*( Gbus(1,2)*sin(TH1-TH2) - Bbus(1,2)*cos(TH1-TH2) ) ; % V1*dQc2_dV1
	Jqv (2,2) = imag(Sc2) - V2^2*Bbus(2,2) ; % V2*dQc2_dV2
    
    % jacobiano
	JACOBIANO = [ 
		[  Jpth   Jpv  ] ;
		[  Jqth   Jqv  ] ;
		]
	% Guardar para RESULTADOS
	iteracionesJACOBIANO(:,:,numITERACIONES) = JACOBIANO ;
    
	% Vector de actualizaciones de las variables
	incX = inv(JACOBIANO)*vecMISMATCHES
    
	% Actualizaciones de las variables
	incTH1   = incX(1)
    incTH2   = incX(2)
	incTH3   = incX(3)
	incV1_V1 = incX(4)
    incV2_V2 = incX(5)
    
	% Guardar para RESULTADOS
	iteracionesINCX(:,numITERACIONES) = incX ;

	% Actualizar variables
    TH1 = TH1 + incTH1
	TH2 = TH2 + incTH2
	TH3 = TH3 + incTH3
	V1 = V1*( 1 + incV1_V1)
    V2 = V2*( 1 + incV2_V2)
	% Fasores y vector de tensiones en forma fasorial (OJO: vector columna)
	fasV1 = V1*exp(J*TH1) ;
    fasV2 = V2*exp(J*TH2) ;
	fasV3 = V3*exp(J*TH3) ;
    fasV4 = V4*exp(J*TH4) ;
	vecV = [ fasV1 ; fasV2 ; fasV3 ; fasV4 ]
    
	% Guardar para RESULTADOS
	iteracionesX(:,numITERACIONES) = [ TH1 ; TH2 ; TH3 ; V1 ; V2 ] ;

	% COMPROBACION MISMATCHES
	% Vector de corrientes netas inyectadas
	vecI = Ybus*vecV
    
	% Potencias netas inyectadas, esto es, potencias calculadas (Sc=Pc+J*Qc)
	Sc1 = vecV(1)*conj(vecI(1))
	Sc2 = vecV(2)*conj(vecI(2))
	Sc3 = vecV(3)*conj(vecI(3))
    Sc4 = vecV(4)*conj(vecI(4))
    
	% En forma de vector
	Sc = [ Sc1 ; Sc2 ; Sc3 ; Sc4] ;
    
	% Mismatches, segun tabla:
    MP1 = Pe1 - real(Sc1)
	MP2 = Pe2 - real(Sc2)
	MP3 = Pe3 - real(Sc3)
	MQ1 = Qe1 - imag(Sc1)
    MQ2 = Qe2 - imag(Sc2)
    
	% Vector de Mismatches (OJO: vector columna)
	vecMISMATCHES = [ MP1 ; MP2 ; MP3 ; MQ1 ; MQ2 ]
    
	% Guardar para RESULTADOS
	iteracionesMISMATCHES(:,numITERACIONES) = vecMISMATCHES ;
    
	% Error MAXIMO
	errmaxMISMATCHES = max(abs(vecMISMATCHES)) 

end

%% RESULTADOS
echo off
disp('   ');
disp('Evolucion de los MISMATCHES');
disp([ 'iter       ' num2str(0:numITERACIONES,'%8d') ])
disp([ ['MP1  ';'MP2  ';'MP3  ';'MQ1  ';'MQ2  '] num2str([ vecMISMATCHES_inicial iteracionesMISMATCHES ],'%8.4f')])
disp('   ');
disp('Evolucion del jacobiano');
for ii=1:numITERACIONES
	disp([ '*** iteracion  ' num2str(ii,'%2d') ])
	disp('        incTH1    incTH2    incTH3  incV1/V1  incV2/V2  ')
	disp([ ['MP1  ';'MP2  ';'MP3  ';'MQ1  ';'MQ2  '] num2str(iteracionesJACOBIANO(:,:,ii),'%10.3f')])
end
disp('   ');
disp('Evolucion del vector de actualizaciones');
disp([ 'iter            ' num2str(1:numITERACIONES,'%8d') ])
disp([ ['incTH1    ';'incTH2    ';'incTH3    ';'incV1/V1  ';'incV2/V2  '] num2str(iteracionesINCX,'%8.4f')])
disp('   ');
disp('Evolucion de las variables');
disp([ 'iter        ' num2str(0:numITERACIONES,'%8d') ])
disp([ ['TH1    ';'TH2    ';'TH3    ';'V1     ';'V2     '] num2str([ X_inicial iteracionesX ],'%8.4f')])
disp('   ');
disp('Resumen Solucion');
disp('n    V (pu)    TH (º)    Pc (MW)   Qc(Mvar) ')

for ii=1:numeroNUDOS
	disp(sprintf('%1d %9.4f %9.3f %10.2f %10.2f',ii,abs(vecV(ii)),angle(vecV(ii))*180/pi,Sbase*real(Sc(ii)),Sbase*imag(Sc(ii))))
end










