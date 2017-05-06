% TRABAJO DE SEP TEMA4: FLUJO DE CARGAS
% Fecha: 27/04/2016
% Autores: IGNACIO SANZ, GUILLERMO MAIRAL, JAVIER MIÑARRO,
% PEDRO PARDO, IGNACIO PASTOR, SANTIAGO LÓPEZ

clear all
format short
disp('  ')
disp('  ')
echo on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTES
J = sqrt(-1) ;
Sbase = 100 ;
numeroNUDOS = 4 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATOS
x12 = 0.23 ;
x14bis = 0.31;
x14 = 0.31 ;
x23 = 0.065 ;
x24 = 0.19;

% Matriz Ybus en DC
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
YbusDC_mismP = YbusDC(1:3,1:3); 
invYbusDC_mismP = YbusDC_mismP\eye(size(YbusDC_mismP)); % A\B = inv(A)*B

% Potencias especificadas
Pe1 =  -210/Sbase; %   Activa nudo 2, en pu
Pe2 =     0/Sbase;
Pe3 =    60/Sbase; %   Activa nudo 3, en pu

% Vector potencias inyectadas
Piny = [
	Pe1 ;
    Pe2 ;
	Pe3 ;
	];

% Ángulos de las tensiones
vecTH = invYbusDC_mismP*Piny;
vecTH = [  vecTH ; 0.0 ]; % Añadir el del slack
TH1 = vecTH(1);
TH2 = vecTH(2);
TH3 = vecTH(3);
TH4 = vecTH(4); % Referencia

disp('Vector de angulos(grados)')
vecTH = [ vecTH ]*180/pi
vecTH = vecTH*pi/180;

%% Análisis de contingencias
% Atención: con invYbusDC_mismP, puesto que posicion no corresponde con nudo

%% 1. Línea 1-2
% Flujo inicial
f12_0 = (TH1-TH2)/x12;

% Potencia ficticia (+ en 1, - en 2)
incP_ctg12 = f12_0/( 1 - (invYbusDC_mismP(1,1) + invYbusDC_mismP(2,2) - invYbusDC_mismP(1,2) - invYbusDC_mismP(2,1) )/x12); % los elementos correspondientes al nudo 1 son nulos
incTH_ctg12 = invYbusDC_mismP*[ incP_ctg12; -incP_ctg12 ; 0.0 ]; 
vecTH_ctg12 = vecTH + [  incTH_ctg12 ; 0.0 ]; % Añadir el del slack

%% 2.Línea 2-4
% Flujo inicial
f24_0 = (TH2-TH4)/x24;

% Potencia ficticia (+ en 2, - en 4)
incP_ctg24 = f24_0/( 1 - ( invYbusDC_mismP(2,2) + 0.0 - 0.0 - 0.0 )/x24); 
incTH_ctg24 = invYbusDC_mismP*[ 0.0 ; incP_ctg24 ; 0.0 ]; 
vecTH_ctg24 = vecTH + [  incTH_ctg24 ; 0.0]; % Añadir el del slack

%% 3. Línea 2-3 
% Flujo inicial
f23_0 = 2*(TH2-TH3)/x23;

% Potencia ficticia (+ en 1, - en 4)
incP_ctg23 = f23_0/( 1 - (invYbusDC_mismP(2,2) + invYbusDC_mismP(3,3) - invYbusDC_mismP(2,3) - invYbusDC_mismP(2,3) )/x14); % los elementos correspondientes al nudo 1 son nulos
incTH_ctg23 = invYbusDC_mismP*[ 0.0 ; incP_ctg23; -incP_ctg23 ]; 
vecTH_ctg23 = vecTH + [  incTH_ctg23 ; 0.0 ]; % Añadir el del slack

% Nota: Las dos líneas 1-4 son iguales, por tanto si se pierde sólo 1 de las dos,
% es indiferente cual se haya perdido

%% 4. Línea 1-4 (si se pierde sólo 1 de las dos)
% Flujo inicial
f14_0 = (TH1-TH4)/x14;

% Potencia ficticia (+ en 1, - en 4)
incP_ctg14 = f14_0/( 1 - (invYbusDC_mismP(1,1) + 0.0 - 0.0 - 0.0 )/x14); % los elementos correspondientes al nudo 1 son nulos
incTH_ctg14 = invYbusDC_mismP*[ incP_ctg14; 0.0 ; 0.0 ]; 
vecTH_ctg14 = vecTH + [  incTH_ctg14 ; 0.0 ]; % Añadir el del slack

%% 5. Línea 1-4 (si se pierden las dos líneas)
% Flujo inicial
f14_0_bis = 2*(TH1-TH4)/x14;

% Potencia ficticia (+ en 1, - en 4)
incP_ctg14_bis = f14_0/( 1 - (invYbusDC_mismP(1,1) + 0.0 - 0.0 - 0.0 )/2*x14); % los elementos correspondientes al nudo 1 son nulos
incTH_ctg14_bis = invYbusDC_mismP*[ +incP_ctg14_bis; 0.0 ; 0.0 ]; 
vecTH_ctg14_bis = vecTH + [  incTH_ctg14_bis ; 0.0 ]; % Añadir el del slack

%% RESULTADOS
echo off
disp('   ');
disp('Resumen Solucion');
disp('   TH [º] ')
disp('n    base    linea 12  linea 23  linea 24  linea 14  linea 14x2 ')
for ii=1:4
	disp(sprintf('%1d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f ',ii,vecTH(ii)*180/pi,vecTH_ctg12(ii)*180/pi,vecTH_ctg23(ii)*180/pi,vecTH_ctg24(ii)*180/pi,vecTH_ctg14(ii)*180/pi,vecTH_ctg14_bis(ii)*180/pi))
end










