% TRABAJO DE SEP TEMA1: PARÁMETROS Y MODELOS DE LÍNEAS
% Fecha: 07/02/2016
% Autores: IGNACIO SANZ, GUILLERMO MAIRAL, JAVIER MIÑARRO,
% PEDRO PARDO, IGNACIO PASTOR, SANTIAGO LÓPEZ

%%% CONSTANTES
mu0  = 4*pi*1.0E-7;
eps0 = 8.854187817E-12;
J = sqrt(-1);

%%% DATOS
RROterreno = 110.0;
frecuencia = 50;
L = 9200; % Longitud de la línea (m)

% TERRENO
rg = 9.869E-7*frecuencia;
De  =658.368*sqrt(RROterreno/frecuencia);

% Coordenadas (conductor y pantalla tienen las mismas)
R_cubierta = 32.6e-3; % Radio de la cubierta (m)
Distancia_entre_cubiertas = 30e-3; % (m)
D_ejes = 2*R_cubierta + Distancia_entre_cubiertas; % separacion entre ejes de las fases en m
x = +0.0 + (D_ejes/sqrt(3))*[ cos(+90*pi/180) cos(-30*pi/180) cos(-150*pi/180) ] ;
y = +0.0 + (D_ejes/sqrt(3))*[ sin(+90*pi/180) sin(-30*pi/180) sin(-150*pi/180) ] ;
RMGs_cond   = [ 12.4 12.4 12.4 ]*1.0E-3 ; % A metros
RADIOs_pant = [ 29.5 29.5 29.5 ]*1.0E-3 ; % A metros
RESISTENCIAs_cond = [  0.045 0.045 0.045 ]*1.0E-3 ; % A OHM/m
RESISTENCIAs_pant = [ 0.136 0.136 0.136 ]*1.0E-3 ; % A OHM/m
epsr = 2.6 ;
RADIO_cond =  15.6e-3 ;
RADIO_aisl =  25.8e-3 ;

%% 1. ZSERIE: separamos cond-cond, pant-pant y cond-pant
for fila=1:3
	for colum=1:3
		if fila==colum 
			Zserie_cc(fila,colum) = RESISTENCIAs_cond(fila) + rg + J*(mu0/(2*pi))*2*pi*frecuencia*log(De/RMGs_cond(fila)) ;
			Zserie_cp(fila,colum) =                         + rg + J*(mu0/(2*pi))*2*pi*frecuencia*log(De/RADIOs_pant(fila)) ;
			Zserie_pp(fila,colum) = RESISTENCIAs_pant(fila) + rg + J*(mu0/(2*pi))*2*pi*frecuencia*log(De/RADIOs_pant(fila)) ;
		else 
			distancia = sqrt( (x(fila)-x(colum))^2 + (y(fila)-y(colum))^2 ) ;
			Zserie_cc(fila,colum) = rg + J*(mu0/(2*pi))*2*pi*frecuencia*log(De/distancia) ;
			Zserie_cp(fila,colum) = rg + J*(mu0/(2*pi))*2*pi*frecuencia*log(De/distancia) ;
			Zserie_pp(fila,colum) = rg + J*(mu0/(2*pi))*2*pi*frecuencia*log(De/distancia) ;
		end
	end
end

Zserie_pc = transpose(Zserie_cp); 
Zserie = [ 
	Zserie_cc    Zserie_cp ;
	Zserie_pc    Zserie_pp ;
	]

%% 2. Eliminación de conductores pasivos
% Condición de contorno ambas extremos a tierra
disp('Zserie sin Conductores Pasivos en OHM/Km')
Zserie_sinP = (Zserie_cc - Zserie_cp*(Zserie_pp\Zserie_pc))*1e3

%% 3. Matriz de Impedancias de Conductores Equivalentes
% Monofasico equivalente
ROT = [ 0 1 0 ; 0 0 1 ; 1 0 0 ]; 
invROT = ROT';  % La inversa es la traspuesta
Zserie_sinP_equilibrada = (1/3)*( Zserie_sinP + ROT*Zserie_sinP*invROT + invROT*Zserie_sinP*ROT );
ZD = Zserie_sinP_equilibrada(1,1);  % Tambien valen (2,2) o (3,3)
ZM = Zserie_sinP_equilibrada(1,2);  % Tambien valen (1,3), (2,3), (3,2), (2,1) o (3,1)

disp('Zserie Monofasico Equivalente en OHM/Km')
Zmonofasicoequivalente = ZD - ZM

%% 4. Matriz de Admitancias Paralelo
Ccable = 2*pi*eps0*epsr/log(RADIO_aisl/RADIO_cond); 
Cparalelo_cc = [ Ccable 0 0 ; 0 Ccable 0 ; 0 0 Ccable ]; 
Cparalelo_pp = +Cparalelo_cc; 
Cparalelo_cp = -Cparalelo_cc; 
Cparalelo_pc = -Cparalelo_cc; 
Cparalelo = [
	Cparalelo_cc Cparalelo_cp ;
	Cparalelo_pc Cparalelo_pp ;
];

disp('Y paralelo de conductores, en micorS/km')
Yparalelo = J*Cparalelo*2*pi*frecuencia*1e9

%% 5. Matriz de Admitancias de Conductores equivalentes
% Matriz de conductores
Yparalelo_sinP = J*Cparalelo_cc*2*pi*frecuencia*1e9;

% Monofasico equivalente
ROT = [ 0 1 0 ; 0 0 1 ; 1 0 0 ]; 
invROT = ROT';  % La inversa es la traspuesta
Yparalelo_sinP_equilibrada = (1/3)*( Yparalelo_sinP + ROT*Yparalelo_sinP*invROT + invROT*Yparalelo_sinP*ROT );
YD = Yparalelo_sinP_equilibrada(1,1)  % Tambien valen (2,2) o (3,3)
YM = Yparalelo_sinP_equilibrada(1,2)  % Tambien valen (1,3), (2,3), (3,2), (2,1) o (3,1)
disp('Yparalelo Monofasico Equivalente en microS/Km')
Ymonofasicoequivalente = YD - YM

%% 6. Monofasicos equivalentes simplificados
D_ejes = 2*R_cubierta + Distancia_entre_cubiertas; % separacion entre ejes de las fases en m
xfases = +0.0 + (D_ejes/sqrt(3))*[ cos(+90*pi/180) cos(-30*pi/180) cos(-150*pi/180) ]; 
yfases = +0.0 + (D_ejes/sqrt(3))*[ sin(+90*pi/180) sin(-30*pi/180) sin(-150*pi/180) ]; 
RMG   = 12.4*1.0E-3;  % A metros
RESISTENCIA = 0.045*1.0E-3;  % A OHM/m
DRS = sqrt( (xfases(1)-xfases(2))^2 + (yfases(1)-yfases(2))^2 ); 
DST = sqrt( (xfases(2)-xfases(3))^2 + (yfases(2)-yfases(3))^2 ); 
DTR = sqrt( (xfases(3)-xfases(1))^2 + (yfases(3)-yfases(1))^2 ); 
DMG = (DRS*DST*DTR)^(1/3); 
RMGEQ = RMG; 

% Resistencia
RESISTENCIA_EQ = RESISTENCIA/1;  % 1 conductores por fase

% Reactancia
REACTANCIA_EQ = (mu0/(2*pi))*2*pi*frecuencia*log(DMG/RMGEQ); 

% Susceptancia: cable aislado
SUSCEPTANCIA_EQ = 2*pi*frecuencia*Ccable; 

% Z e Y
Zmonofasicoequivalente_simp = RESISTENCIA_EQ + J*REACTANCIA_EQ;
Ymonofasicoequivalente_simp =                + J*SUSCEPTANCIA_EQ;

disp('Zserie de conductores Metodo Simplificado, en OHM/km')
Zmonofasicoequivalente_simp = (RESISTENCIA_EQ + j*REACTANCIA_EQ)*1.0e3
disp('Yparalelo de conductores Metodo Simplificado, en microS/km')
Ymonofasicoequivalente_simp =                 + j*SUSCEPTANCIA_EQ*1.0e9


%%%%%%%%%%%% PARÁMETROS DISTRIBUIDOS %%%%%%%%%%%%
%% 7.1. Modelo equivalente en pi de la linea: Real y Simplificado.
% a) Impedancia Caracteristica
Zs = Zmonofasicoequivalente*1e-3; % Ohm/m
Yp = Ymonofasicoequivalente*1e-9; % S/m
Zc = sqrt(Zs/Yp);

% b)Constante de Propagacion
gamma = sqrt(Zs*Yp);

% c) Modelo Real de parametros distribuidos
% Elementos de la matriz de transferencia (A,B,C,D)
A = cosh(gamma*L);
B = -Zc*sinh(gamma*L);
C = -sinh(gamma*L)/Zc;
D = cosh(gamma*L);
M = [A B;
    C D];

% d) Modelo exacto en Pi de la línea
Zpi = Zc*sinh(gamma*L);
Ypi = 2*tanh(gamma*L/2)/Zc;

% e) Elementos de la matriz de transferencia (A,B,C,D)
A = 1 + Zpi*Ypi/2;
B = -Zpi;
C = -Ypi/2*(2 + Zpi*Ypi/2);
D = 1 + Zpi*Ypi/2;
disp('Matriz de transferencia del modelo en PI completo')
M_PI_comp = [A B;
    C D]

%% 7.2. Modelo equivalente en pi de la linea: Simplificado.
Zpi = Zmonofasicoequivalente*1e-3*L;
Ypi = Ymonofasicoequivalente*1e-9*L;

% e) Elementos de la matriz de transferencia (A,B,C,D)
A = 1 + Zpi*Ypi/2;
B = -Zpi;
C = -Ypi/2*(2+Zpi*Ypi/2);
D = 1 + Zpi*Ypi/2;
disp('Matriz de transferencia del modelo en PI simplificado')
M_PI_simp = [A B;
    C D]

%% 8.1.Con 40 MW en la carga, Tensión necesaria en el generador y la potencia producida por este. Obtener el balance de potencias en la línea.
% Matriz de transferencia. Línea corta, sin despreciar la admitancia paralelo
Zs_prima = Zmonofasicoequivalente*1e-3; % Ohm/m
Yp_prima = Ymonofasicoequivalente*1e-9; % S/m

% Impedancia Caracteristica
ZC = sqrt(Zs_prima/Yp_prima);

% Constante de propagación
gamma = sqrt(Zs_prima*Yp_prima);
Zs_a = Zs_prima*L;
Yp_a = Yp_prima*L;

% Matriz de transferencia
ABCD_A = 1 + Zs_a*Yp_a/2;
ABCD_B = -Zs_a;
ABCD_C = -(Yp_a/2)*(2 + Zs_a*Yp_a/2);
ABCD_D = 1 + Zs_a*Yp_a/2;
ABCD_a2 = [
	ABCD_A  ABCD_B ;
	ABCD_C  ABCD_D ;
    ];

% Condiciones de contorno de la linea en la carga
P_load = 40e6;                                    % Potencia en la carga
V_load = 66000/sqrt(3);                           % Tensión monofásica en la carga
S_load= (P_load/0.95)*(0.95+j*sin(acos(0.95)))/3; % Potencia aparente monofásica en la carga
I_load = (S_load/V_load)';                        % Intensidad en la carga
v_load = [V_load I_load];                         % Vector v_linea
v_load_hip_1 = v_load;                            % Guardamos el vector para después

% Tensión e intensidad en el generador
disp('Tensión e intensidad en bornes del generador (40 MW)')
v_0 = inv(ABCD_a2)*v_load'; % Vector V_0 en el generador

% Resultados
fprintf('Tensión en bornes del generador (40 MW)\n %.2f /_%.2fº kV\n',abs(v_0(1))*1e-3,angle(v_0(1))*180/pi)
fprintf('Intensidad en bornes del generador (40 MW)\n %.2f /_%.2fº A\n',abs(v_0(2)),angle(v_0(2))*180/pi)

%% 8.2. Potencia producida por el generador
% Generador
V_generador = v_0(1);                     % Tensión Monofasica
I_generador = v_0(2);                     % Intensidad de Linea
S_generador = 3*V_generador*I_generador'; % Potencia trifásica
fprintf('Potencia en bornes del generador\n %.2f /_%.2fº MVA \n',abs(S_generador)*1e-6,angle(S_generador)*180/pi)

% Potencia trifasica de la carga
S_3load = 3*S_load;                       

% Balance de Potencias
disp('Balance de Potencias. Potencia consumida en la linea')
S_linea = S_generador - S_3load;             
fprintf('Potencia Activa\n %.2f MW\n',real(S_linea)*1e-6); 
fprintf('Potencia Reactiva\n %.2f MVar\n',imag(S_linea)*1e-6);

%% 9.1.Con 140 MW en la carga, Tensión necesaria en el generador y la potencia producida por este. Obtener el balance de potencias en la línea.
% Condiciones de contorno de la linea en la carga
P_load = 140e6;                                   % Potencia en la carga
V_load = 66000/sqrt(3);                           % Tensión monofásica en la carga
S_load= (P_load/0.95)*(0.95+j*sin(acos(0.95)))/3; % Potencia aparente monofásica en la carga
I_load = (S_load/V_load)';                        % Intensidad en la carga
v_load = [ V_load I_load ];                       % Vector v_linea
v_load_hip_2 = v_load;                            % Guardamos el vector para después

% Tensión e intensidad en el generador
disp('Tensión e intensidad en bornes del generador (140 MW)')
v_0 = inv(ABCD_a2)*v_load'; 

% Resultados
fprintf('Tensión en bornes del generador (140 MW)\n %.2f /_%.2fº kV\n',abs(v_0(1))*1e-3,angle(v_0(1))*180/pi)
fprintf('Intensidad en bornes del generador (140 MW)\n %.2f /_%.2fº A\n',abs(v_0(2)),angle(v_0(2))*180/pi)

%% 9.2. Potencia producida por el generador
% Generador
V_generador = v_0(1);                     % Tensión Monofasica
I_generador = v_0(2);                     % Intensidad de Linea
S_generador = 3*V_generador*I_generador'; % Potencia trifásica
fprintf('Potencia en bornes del generador\n %.2f /_%.2fº MVA \n',abs(S_generador)*1e-6,angle(S_generador)*180/pi)

% Potencia trifasica de la carga
S_3load = 3*S_load; 

% Balance de Potencias
disp('Balance de Potencias. Potencia consumida en la linea')
S_linea = S_generador - S_3load; 
fprintf('Potencia Activa\n %.2f MW\n',real(S_linea)*1e-6); 
fprintf('Potencia Reactiva\n %.2f MVar\n',imag(S_linea)*1e-6);

%% 10. La intensidad por los conductores pasivos en las hipotesis de carga 4) y 5).
% Intensidad por la linea ddel apartado 4
% Para construir el vector de intensidades por los conductores, utilizamos
% la intensidad del monofasico equivalente y las desfasamos +-120º en el
% tiempo para obtener las tres fases R,S,T. Ídem para tensiones.

%% 10.1. Hipotesis apartado 4)
a = 1*exp(J*2*pi/3);
I_load_hip_1 = [ v_load_hip_1(2); v_load_hip_1(2)*a ; v_load_hip_1(2)*a^(2) ];
I_serie_hip_1 = -Zserie_pp\Zserie_cp*I_load_hip_1; 
fprintf('Intensidad por las pantallas hipótesis 1\n')
for i=1:3
fprintf('%.2f /_ %.2fº A\n',abs(I_serie_hip_1(i)),angle(I_serie_hip_1(i))*180/pi);
end

%% 10.2. Hipotesis apartado 5)
a = exp(J*2*pi/3);
I_load_hip_2 = [ v_load_hip_2(2); v_load_hip_2(2)*a ; v_load_hip_2(2)*a^(2) ];
I_serie_hip_2 = -Zserie_pp\Zserie_cp*I_load_hip_2;
fprintf('Intensidad por las pantallas hipótesis 2\n')
for i=1:3
fprintf('%.2f /_ %.2fº A\n',abs(I_serie_hip_2(i)),angle(I_serie_hip_2(i))*180/pi);
end

