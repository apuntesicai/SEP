% TRABAJO DE SEP TEMA 3: REDES TRIFÁSICAS DESEQUILIBRADAS
% Fecha: 9/04/2016
% Autores: IGNACIO SANZ, GUILLERMO MAIRAL, JAVIER MIÑARRO,
% PEDRO PARDO, IGNACIO PASTOR, SANTIAGO LÓPEZ

clear all
format short
disp('  ')
echo off

% La red de la figura presenta el esquema de evacuación de una central que cuenta con dos grupos de
% generación. En el momento anterior a la falta, ambas máquinas presentan tensión de vacío nominal
% y la red infinita se encuentra a 390kV. Se considera además que todos los grupos presentan ángulo
% cero (con la conveniente corrección debida a su transformador). 
% Los datos de la red son los siguientes:
% ALTERNADOR a: 300 MVA, 35 kV; zg1 = j0.3 p.u., zg2 = j0.2 p.u., zg0 = j0.06 p.u.; ZN = j1.5 ?
% TRANSFORMADOR a: 300 MVA, 395/35 kV, YNd5; zcc = j0.13 p.u.
% ALTERNADOR b: 200 MVA, 20 kV; zg1 = j0.2 p.u., zg2 = j0.1 p.u., zg0 = j0.04 p.u.; ZN = 0.0 ?
% TRANSFORMADOR b: 200 MVA, 410/20 kV, YNyn0; zcc = j0.15 p.u.
% LÍNEA: X1=190, X0=543.
% Falta franca fase-tierra 
% a) Intensidades de avería (kA)
% b) Intensidades que SALEN de la red infinita y de ambos transformadores en alta (kA)
% c) Intensidades de los alternadores (kA)
% d) Intensidades del neutro hacia tierra en alternadores, transformadores y la red infinita (kA)
% e) Tensiones en el punto de la avería (kV)
% f) Tensiones en bornes de los alternadores (kV) 

j = sqrt(-1) ;
a = exp(j*2*pi/3) ;
A = [
	[ 1   1   1  ] ;
	[ 1  a*a  a  ] ;
	[ 1   a  a*a ] ;
	] ;

invA = inv(A) ;

%% 1. Magnitudes Base del Problema
SB = 100/3; % MVA
UB = [400 35 20]/sqrt(3); % kV
IB = SB./UB; % KA
ZB = UB./IB; % Ohm

tB_1 = UB(1)/UB(2);
tB_2 = UB(1)/UB(3);

%% 2. Alternador 1
Sg1    = 300/3;
Ug1    = UB(2)/sqrt(3); 

zg1_G1 = j*0.3;   zg1_1 = zg1_G1*(Ug1*Ug1/Sg1)/ZB(2);
zg2_G1 = j*0.2;   zg2_1 = zg2_G1*(Ug1*Ug1/Sg1)/ZB(2);
zg0_G1 = j*0.06;  zg0_1 = zg0_G1*(Ug1*Ug1/Sg1)/ZB(2);
zN_G1  = j*1.5;   zgN_1 = zN_G1*(Ug1*Ug1/Sg1)/ZB(2);

%% 3. Transformador 1 (YNd5)
St1     = 300/3;
Ut1     = [35 395]/sqrt(3);
zcc1_t1 = j*0.13; zcc_1 = zcc1_t1*(Ut1(2)*Ut1(2)/St1)/ZB(1); %lado de Alta

phB_t1 = [ 0 150 ]*pi/180;
t_1    = (Ut1(2)/Ut1(1))/tB_1;  % lado de Alta, Alta/Baja
t0_1   = t_1 ;
t1_1   = t_1*exp(j*phB_t1(2));
t2_1   = t_1*exp(-j*phB_t1(2));
O12_t1 = [ t0_1 ; t1_1 ; t2_1 ];
ug1    = (Ug1/UB(2)); % tension en el generador 1

%% 4. Alternador 2
Sg2    = 200/3;
Ug2    = UB(3)/sqrt(3); 

zg1_G2 = j*0.2;  zg1_2 = zg1_G2*(Ug2*Ug2/Sg2)/ZB(3);
zg2_G2 = j*0.1;  zg2_2 = zg2_G2*(Ug2*Ug2/Sg2)/ZB(3);
zg0_G2 = j*0.04; zg0_2 = zg0_G2*(Ug2*Ug2/Sg2)/ZB(3);
zN_G2  = j*2.7;  zgN_2 = zN_G2*(Ug2*Ug2/Sg2)/ZB(3);

%% 5. Transformador 2 (YNy0)
St2     = 200/3;
Ut2     = [20 410]/sqrt(3);
zcc2_t2 = j*0.15; zcc_2 = zcc2_t2*(Ut2(2)*Ut2(2)/St2)/ZB(1); % lado de Alta

phB_t2  = [ 0 0 ]*pi/180;
t_2     = (Ut2(2)/Ut2(1))/tB_2;  %  lado de Alta
t0_2    = t_2;
t1_2    = t_2*exp(j*phB_t2(2));
t2_2    = t_2*exp(-j*phB_t2(2));
O12_t2  = [ t0_2 ; t1_2 ; t2_2 ];
ug2     = (Ug2/UB(3)); % tension en el generador 2

%% 6. Linea:
X1 = 190; % Ohm
X0 = 543; % Ohm

z1_linea = j*X1/ZB(1);
z2_linea = z1_linea;
z0_linea = j*X0/ZB(1);

%% 7. Red:
Uinf = 390/sqrt(3);
uinf = Uinf/UB(1);

%% 0. Circuitos equivalentes de Thevenin en el punto F.
% Secuencia Directa (1) 
% Generador1 en alta
Ug1_alta   = ug1*O12_t1(2);
zg1_1_alta = zg1_1/O12_t1(1)^2;

% Generador 2 en alta
Ug2_alta   = ug2*O12_t2(2);
zg1_2_alta = zg1_2/O12_t2(1)^2;

% Z_thevenin
zth1_1_izq = zcc_1 + zg1_1_alta;
zth1_2_izq = zcc_2 + zg1_2_alta;
zth1_der   = z1_linea;
zth1_v     = [zth1_1_izq zth1_2_izq zth1_der];

% Ethevenin
eth1_1_izq = Ug1_alta;
eth1_2_izq = Ug2_alta;
eth1_der   = uinf;
eth1       = [eth1_1_izq eth1_2_izq eth1_der];

% Vthevenin
vth1 = eth1(1)*(paralelo(zth1_v(2),zth1_v(3)))/(paralelo(zth1_v(2),zth1_v(3)) + zth1_v(1));
vth2 = eth1(2)*(paralelo(zth1_v(1),zth1_v(3)))/(paralelo(zth1_v(1),zth1_v(3)) + zth1_v(2));
vth3 = eth1(3)*(paralelo(zth1_v(2),zth1_v(1)))/(paralelo(zth1_v(2),zth1_v(1)) + zth1_v(3));
vth_total = vth1 + vth2 + vth3;

% Equivalente thevenin del conjunto
zth1_total = paralelo(paralelo(zth1_v(1),zth1_v(2)),zth1_v(3));
eth1_total = vth_total;

%%% Secuencia Inversa (2)
% Generador1 + Trafo1
zth2_1_izq = zcc_1 + zg2_1*O12_t1(1)^2;
zth2_der   = z1_linea;

% Generador2 + Trafo2
zth2_2_izq = zcc_2 + zg2_2*O12_t2(1)^2;

% Equivalente thevenin del conjunto
zth2_v     = [zth2_1_izq zth2_2_izq zth2_der];
zth2_total = paralelo(paralelo(zth2_v(1),zth2_v(2)),zth2_v(3));

%%% Secuencia Homopolar (0)
% Generador1 + Trafo1
zth0_1_izq = zcc_1;
zth0_der   = z0_linea;

% Generador2 + Trafo2
zth0_2_izq = zcc_2 + (zg0_2 + zgN_2)*O12_t2(1)^2;

% Equivalente thevenin del conjunto
zth0_v      = [zth0_1_izq zth0_2_izq zth0_der];
zth0_total1 = paralelo(zth0_v(1),zth0_v(2));
zth0_total  = paralelo(zth0_total1,zth0_v(3));

%% 1. Falta Franca Fase-Tierra 
disp('%% 1. Falta Franca Fase-Tierra');
% a) Intensidad de Averia
i_averia     = eth1_total/(zth0_total + zth1_total + zth2_total);
i_av_O12     = [ i_averia; i_averia; i_averia ];
rst_i_averia = A*i_av_O12*IB(1);

disp('a) Intensidades de averia (kA y º)');
disp( [ abs(rst_i_averia)  180*angle(rst_i_averia)/pi ] ) 

%% b) Intensidades que SALEN de la red infinita y de ambos trafos en alta.
%%% Secuencia Homopolar
Vfalta_O = -i_averia*zth0_total;
i_red_O  = -Vfalta_O/zth0_der;

i_alta_O   = i_averia-i_red_O;
i_alta_O_1 = -Vfalta_O/zth0_1_izq;
i_alta_O_2 = -Vfalta_O/zth0_2_izq;

%%% Secuencia Inversa
Vfalta_2 = -i_averia*zth2_total;
i_red_2  = -Vfalta_2/zth2_der;

i_alta_2   = i_averia-i_red_2;
i_alta_2_1 = -Vfalta_2/zth2_1_izq;
i_alta_2_2 = -Vfalta_2/zth2_1_izq;

%%% Directa
Vfalta_1 = -Vfalta_2-Vfalta_O;
i_red_1  = (uinf - Vfalta_1)/zth1_der;

i_alta_1   = i_averia-i_red_1;
i_alta_1_1 = (eth1_1_izq - Vfalta_1)/zth1_1_izq;
i_alta_1_2 = (eth1_2_izq - Vfalta_1)/zth1_2_izq;;

%% b.1) Intensidades de Red
O12_i_red = [i_red_O; i_red_1; i_red_2];
rst_i_red = A*O12_i_red*IB(1);

disp('b.1) Intensidades que SALEN de la Red (kA y º)');
disp( [ abs(rst_i_red) 180*angle(rst_i_red)/pi ] )

%% b.2) Intensidades que salen del trafo de alta en los trafos t1 y t2.
O12_i_alta_t1 = [ i_alta_O_1; i_alta_1_1; i_alta_2_1 ];
O12_i_alta_t2 = [ i_alta_O_2; i_alta_1_2; i_alta_2_2 ];
rst_i_alta_t1 = A*O12_i_alta_t1*IB(1);
rst_i_alta_t2 = A*O12_i_alta_t2*IB(1);

disp('b.2.1) Intensidades de alta que salen del Transformador 1 (kA y º)');
disp( [ abs(rst_i_alta_t1) 180*angle(rst_i_alta_t1)/pi ] )
disp('b.2.2) Intensidades de alta que salen del Transformador 2 (kA y º)');
disp( [ abs(rst_i_alta_t2) 180*angle(rst_i_alta_t2)/pi ] )

%% c) Intensidades de los alternadores.
O12_i_generador_1    = O12_i_alta_t1.*conj(O12_t1);
O12_i_generador_1(1) = 0;
rst_i_generador_1    = A*O12_i_generador_1*IB(2);

O12_i_generador_2 = O12_i_alta_t2.*conj(O12_t2);
rst_i_generador_2 = A*O12_i_generador_2*IB(3);

disp('c.1) Intensidades que salen del generador 1 (kA y º)');
disp( [ abs(rst_i_generador_1) 180*angle(rst_i_generador_1)/pi ] )
disp('c.2) Intensidades que salen del generador 2 (kA y º)');
disp( [ abs(rst_i_generador_2) 180*angle(rst_i_generador_2)/pi ] )

%% d) Intensidades del neutro hacia tierra en alternadores, transformadores y la red infinita (kA).
rst_i_neutro_red     = 3*i_red_O*IB(1);
rst_i_neutro_alta_t1 = 0;                            % Transformador en triángulo IO = 0
rst_i_neutro_alta_t2 = 3*O12_i_generador_2(1)*IB(3); % IN = 3*IO
rst_i_neutro_gen1    = 3*O12_i_generador_1(1)*IB(2); % Transformador en triángulo IO = 0
rst_i_neutro_gen2    = 3*O12_i_generador_2(1)*IB(3);

disp('d.1) Intensidades del neutro hacia tierra en  la red infinita (kA).');
disp( [ abs(rst_i_neutro_red) 180*angle(rst_i_neutro_red)/pi ] )
disp('d.2) Intensidades del neutro hacia tierra en  el transformador 1 (kA).');
disp( [ abs(rst_i_neutro_alta_t1) 180*angle(rst_i_neutro_alta_t1)/pi ] )
disp('d.3) Intensidades del neutro hacia tierra en  el transformador 2(kA).');
disp( [ abs(rst_i_neutro_alta_t2) 180*angle(rst_i_neutro_alta_t2)/pi ] )
disp('d.4) Intensidades del neutro hacia tierra en  el generador 1 (kA).');
disp( [ abs(rst_i_neutro_gen1) 180*angle(rst_i_neutro_gen1)/pi ] )
disp('d.5) Intensidades del neutro hacia tierra en  el generaodr 2 (kA).');
disp( [ abs(rst_i_neutro_gen2) 180*angle(rst_i_neutro_gen2)/pi ] )

%% e) Tension en el punto de averia.
O12_Vfalta = [Vfalta_O; Vfalta_1; Vfalta_2];
rst_Vfalta = A*O12_Vfalta*UB(1);

disp('e) Tension en el punto de averia(kV y º).');
disp( [ abs(rst_Vfalta) 180*angle(rst_Vfalta)/pi ] )

%% f) Tensiones en Bornas del Alternador.
O12_Vgen1_bornas = (O12_Vfalta + zcc_1*O12_i_alta_t1)./O12_t1;
rst_Vgen1_bornas = A*(O12_Vgen1_bornas)*UB(2);
O12_Vgen2_bornas = (O12_Vfalta + zcc_2*O12_i_alta_t2)./O12_t2;
rst_Vgen2_bornas = A*(O12_Vgen2_bornas)*UB(3);

disp('f.1) Tensiones en bornas del generador 1 (kV y º)');
disp( [ abs(rst_Vgen1_bornas) 180*angle(rst_Vgen1_bornas)/pi ] )
disp('f.2) Tensiones en bornas del generador 2 (kV y º)');
disp( [ abs(rst_Vgen2_bornas) 180*angle(rst_Vgen2_bornas)/pi ] )


%% 2. Falta Franca Fase-Fase
disp('%% 2. Falta Franca Fase-Fase');

%% a)Intensidades de averia 
% Secuencia Directa (1)
i_averia_ff_1 = eth1_total/(zth1_total + zth2_total);
% Secuencia Inversa (2)
i_averia_ff_2 = -i_averia_ff_1; 
% Secuencia Homopolar (O)
i_averia_ff_O = 0;

O12_i_averia_ff = [ i_averia_ff_O;i_averia_ff_1;i_averia_ff_2 ];
rst_i_averia_ff = A*O12_i_averia_ff*IB(1);

disp('a) Intensidades de averia (kA y º)');
disp( [ abs(rst_i_averia_ff) 180*angle(rst_i_averia_ff)/pi ] ) 

% Tension de Falta Secuencia Directa (1)
Vfalta1_ff  = eth1_total - i_averia_ff_1*zth1_total;
i_red_ff_1  = (uinf - Vfalta1_ff)/zth1_der;
i_alta_ff_1 = i_averia_ff_1 - i_red_ff_1;

% Tension de Falta Secuencia Inversa (2)
Vfalta2_ff  = i_averia_ff_1*zth2_total;
i_red_ff_2  = -Vfalta2_ff/zth2_der;
i_alta_ff_2 = i_averia_ff_2 - i_red_ff_2;

% Tensión de Falta Homopolar (O)
VfaltaO_ff  = 0;
i_red_ff_O  = 0;
i_alta_ff_O = 0;

%% b) Intensidades que SALEN de la red infinita y de ambos trafos en alta.
% b.1) Intensidades de Red
O12_i_red_ff = [0; i_red_ff_1; i_red_ff_2];
rst_i_red_ff = A*O12_i_red_ff*IB(1);

disp('b.1) Intensidades que SALEN de la Red (kA y º)');
disp( [ abs(rst_i_red_ff) 180*angle(rst_i_red_ff)/pi ] )

% b.2) Intensidades que Salen del trafo de alta en los trafos t1 y t2.
% Secuencia Directa (1)
i_alta_ff_1_1 = (eth1_1_izq - Vfalta1_ff)/zth1_1_izq;
i_alta_ff_1_2 = (eth1_2_izq - Vfalta1_ff)/zth1_2_izq;

% Secuencia Inversa (2)
i_alta_ff_2_1 = -Vfalta2_ff/(zth2_1_izq);
i_alta_ff_2_2 = -Vfalta2_ff/(zth2_2_izq);

% Secuencia Homopolar (O)
i_alta_ff_O_1 = 0;
i_alta_ff_O_2 = 0;

% Trafo 1
O12_i_alta_ff_t1 = [ i_alta_ff_O_1; i_alta_ff_1_1; i_alta_ff_2_1 ]; %h,D,I

% Trafo 2
O12_i_alta_ff_t2 = [ i_alta_ff_O_2; i_alta_ff_1_2; i_alta_ff_2_2 ]; %h,D,I

rst_i_alta_ff_t1 = A*O12_i_alta_ff_t1*IB(1);
rst_i_alta_ff_t2 = A*O12_i_alta_ff_t2*IB(1);

disp('b.2.1) Intensidades de alta que salen del Transformador 1 (kA y º)');
disp( [ abs(rst_i_alta_ff_t1) 180*angle(rst_i_alta_ff_t1)/pi ] )
disp('b.2.2) Intensidades de alta que salen del Transformador 2 (kA y º)');
disp( [ abs(rst_i_alta_ff_t2) 180*angle(rst_i_alta_ff_t2)/pi ] )

%% c) Intensidades de los alternadores.
O12_i_generador_ff_1 = O12_i_alta_ff_t1.*conj(O12_t1);
rst_i_generador_ff_1 = A*O12_i_generador_ff_1*IB(2);
O12_i_generador_ff_2 = O12_i_alta_ff_t2.*conj(O12_t2);
rst_i_generador_ff_2 = A*O12_i_generador_ff_2*IB(3);

disp('c.1) Intensidades que salen del generador 1 (kA y º)');
disp( [ abs(rst_i_generador_ff_1) 180*angle(rst_i_generador_ff_1)/pi ] )
disp('c.2) Intensidades que salen del generador 2 (kA y º)');
disp( [ abs(rst_i_generador_ff_2) 180*angle(rst_i_generador_ff_2)/pi ] )

%% d) Intensidades del neutro hacia tierra en alternadores, transformadores y la red infinita (kA).
rst_i_neutro_red_ff     = 3*i_red_ff_O;
rst_i_neutro_alta_ff_t1 = 3*i_alta_ff_O;
rst_i_neutro_alta_ff_t2 = 3*i_alta_ff_O;
rst_i_neutro_gen1_ff    = 3*O12_i_generador_ff_1(1);
rst_i_neutro_gen2_ff    = 3*O12_i_generador_ff_2(1);

disp('d.1) Intensidades del neutro hacia tierra en  la red infinita (kA).');
disp( [ abs(rst_i_neutro_red_ff) 180*angle(rst_i_neutro_red_ff)/pi ] )
disp('d.2) Intensidades del neutro hacia tierra en  el transformador 1 (kA).');
disp( [ abs(rst_i_neutro_alta_ff_t1) 180*angle(rst_i_neutro_alta_ff_t1)/pi ] )
disp('d.3) Intensidades del neutro hacia tierra en  el transformador 2(kA).');
disp( [ abs(rst_i_neutro_alta_ff_t2) 180*angle(rst_i_neutro_alta_ff_t2)/pi ] )
disp('d.4) Intensidades del neutro hacia tierra en  el generador 1 (kA).');
disp( [ abs(rst_i_neutro_gen1_ff) 180*angle(rst_i_neutro_gen1_ff)/pi ] )
disp('d.5) Intensidades del neutro hacia tierra en  el generaodr 2 (kA).');
disp( [ abs(rst_i_neutro_gen2_ff) 180*angle(rst_i_neutro_gen2_ff)/pi ] )

%% e) Tension en el punto de averia.
O12_Vfalta_ff = [ VfaltaO_ff; Vfalta1_ff; Vfalta2_ff];
rst_Vfalta_ff = A*O12_Vfalta_ff*UB(1)*sqrt(3);

disp('e) Tension en el punto de averia(kV y º).');
disp( [ abs(rst_Vfalta_ff) 180*angle(rst_Vfalta_ff)/pi ] )

%% f) Tensiones en Bornas del Alternador.
O12_Vgen1_bornas_ff = (O12_Vfalta_ff + zcc_1*O12_i_alta_ff_t1)./O12_t1;
rst_Vgen1_bornas_ff = A*(O12_Vgen1_bornas_ff.*O12_t1)*UB(2);

O12_Vgen2_bornas_ff = (O12_Vfalta_ff + zcc_2*O12_i_alta_ff_t2)./O12_t2;
rst_Vgen2_bornas_ff = A*(O12_Vgen2_bornas_ff.*O12_t2)*UB(3);

disp('f.1) Tensiones en bornas del generador 1 (kV y º)');
disp( [ abs(rst_Vgen1_bornas_ff) 180*angle(rst_Vgen1_bornas_ff)/pi ] )
disp('f.2) Tensiones en bornas del generador 2 (kV y º)');
disp( [ abs(rst_Vgen2_bornas_ff) 180*angle(rst_Vgen2_bornas_ff)/pi ] )


%% 3. Falta Franca Fase-Fase-Tierra
disp('%% 3. Falta Franca Fase-Fase-Tierra')

%% a) Intensidades de averia 
% Secuencia Directa (1)
i_averia_fft_1 = eth1_total/(zth1_total+paralelo(zth2_total,zth0_total)); 
Vfalta_fft_1   = eth1_total - i_averia_fft_1*zth1_total;

% Secuencia Inversa (2)
Vfalta_fft_2   = Vfalta_fft_1;
i_averia_fft_2 = -Vfalta_fft_2/zth2_total;

% Secuencia Homopolar (O)
Vfalta_fft_O   =  Vfalta_fft_1;
i_averia_fft_O = -Vfalta_fft_O/zth0_total;

O12_i_averia_fft = [ i_averia_fft_O; i_averia_fft_1; i_averia_fft_2 ];
rst_i_averia_fft = A*O12_i_averia_fft*IB(1);

disp('a) Intensidades de averia (kA y º)');
disp( [ abs(rst_i_averia_fft) 180*angle(rst_i_averia_fft)/pi ] ) 

%% b) Intensidades que SALEN de la red infinita y de ambos trafos en alta.
% b.1) Intensidades de Red
i_red_fft_1 = -Vfalta_fft_1/zth1_der;
i_red_fft_2 = -Vfalta_fft_2/zth2_der;
i_red_fft_O = -Vfalta_fft_O/zth0_der;

i_alta_fft_1 = i_averia_fft_1 - i_red_fft_1;
i_alta_fft_2 = i_averia_fft_2 - i_red_fft_2;
i_alta_fft_O = i_averia_fft_O - i_red_fft_O;

O12_i_red_fft = [ i_red_fft_O; i_red_fft_1;i_red_fft_2 ];
rst_i_red_fft = A*O12_i_red_fft*IB(1);

disp('b.1) Intensidades que SALEN de la Red (kA y º)');
disp( [ abs(rst_i_red_fft) 180*angle(rst_i_red_fft)/pi ] )

% Secuencia Directa (1)
i_alta_fft_1_1 = (eth1_1_izq - Vfalta_fft_1)/(zth1_1_izq);
i_alta_fft_1_2 = (eth1_2_izq - Vfalta_fft_1)/(zth1_2_izq);

% Secuencia Inversa (2)
i_alta_fft_2_1 = -Vfalta_fft_2/(zth2_1_izq);
i_alta_fft_2_2 = -Vfalta_fft_2/(zth2_2_izq);

% Secuencia Homopolar (O)
i_alta_fft_O_1 = -Vfalta_fft_2/(zth0_1_izq);
i_alta_fft_O_2 = -Vfalta_fft_2/(zth0_2_izq);

% Trafo 1
O12_i_alta_fft_t1 = [ i_alta_fft_O_1; i_alta_fft_1_1; i_alta_fft_2_1 ]; %h,D,I

% Trafo 2
O12_i_alta_fft_t2 = [ i_alta_fft_O_2; i_alta_fft_1_2; i_alta_fft_2_2 ]; %h,D,I

rst_i_alta_fft_t1 = A*O12_i_alta_fft_t1*IB(1);
rst_i_alta_fft_t2 = A*O12_i_alta_fft_t2*IB(1);

disp('b.2.1) Intensidades de alta que salen del Transformador 1 (kA y º)');
disp( [ abs(rst_i_alta_fft_t1) 180*angle(rst_i_alta_fft_t1)/pi ] )
disp('b.2.2) Intensidades de alta que salen del Transformador 2 (kA y º)');
disp( [ abs(rst_i_alta_fft_t2) 180*angle(rst_i_alta_fft_t2)/pi ] )

%% c) Intensidades de los alternadores.
O12_i_generador_fft_1    = O12_i_alta_fft_t1.*conj(O12_t1);
O12_i_generador_fft_1(1) = 0;
rst_i_generador_fft_1    = A*O12_i_generador_fft_1*IB(2);

O12_i_generador_fft_2 = O12_i_alta_fft_t2.*conj(O12_t2);
rst_i_generador_fft_2 = A*O12_i_generador_fft_2*IB(3);

disp('c.1) Intensidades que salen del generador 1 (kA y º)');
disp( [ abs(rst_i_generador_fft_1) 180*angle(rst_i_generador_fft_1)/pi ] )
disp('c.2) Intensidades que salen del generador 2 (kA y º)');
disp( [ abs(rst_i_generador_fft_2) 180*angle(rst_i_generador_fft_2)/pi ] )

%% d) Intensidades del neutro hacia tierra en alternadores, transformadores y la red infinita (kA).
rst_i_neutro_red_fft     = 3*i_red_fft_O;
rst_i_neutro_alta_fft_t1 = 3*i_alta_fft_O_1;
rst_i_neutro_alta_fft_t2 = 3*i_alta_fft_O_2;
rst_i_neutro_gen1_fft    = 3*O12_i_generador_fft_1(1);
rst_i_neutro_gen2_fft    = 3*O12_i_generador_fft_2(1);

disp('d.1) Intensidades del neutro hacia tierra en  la red infinita (kA).');
disp( [ abs(rst_i_neutro_red_fft) 180*angle(rst_i_neutro_red_fft)/pi ] )
disp('d.2) Intensidades del neutro hacia tierra en  el transformador 1 (kA).');
disp( [ abs(rst_i_neutro_alta_fft_t1) 180*angle(rst_i_neutro_alta_fft_t1)/pi ] )
disp('d.3) Intensidades del neutro hacia tierra en  el transformador 2(kA).');
disp( [ abs(rst_i_neutro_alta_fft_t2) 180*angle(rst_i_neutro_alta_fft_t2)/pi ] )
disp('d.4) Intensidades del neutro hacia tierra en  el generador 1 (kA).');
disp( [ abs(rst_i_neutro_gen1_fft) 180*angle(rst_i_neutro_gen1_fft)/pi ] )
disp('d.5) Intensidades del neutro hacia tierra en  el generaodr 2 (kA).');
disp( [ abs(rst_i_neutro_gen2_fft) 180*angle(rst_i_neutro_gen2_fft)/pi ] )

%% e) Tension en el punto de averia.
O12_Vfalta_fft = [ Vfalta_fft_O; Vfalta_fft_1; Vfalta_fft_2 ];
rst_Vfalta_fft = A*O12_Vfalta_fft*UB(1);

disp('e) Tension en el punto de averia(kV y º).');
disp( [ abs(rst_Vfalta_fft) 180*angle(rst_Vfalta_fft)/pi ] )

%% f) Tensiones en Bornas del Alternador.
O12_Vgen1_bornas_fft = (O12_Vfalta_fft + zcc_1*O12_i_alta_fft_t1)./O12_t1;
rst_Vgen1_bornas_fft = A*(O12_Vgen1_bornas_fft)*UB(2);

O12_Vgen2_bornas_fft = (O12_Vfalta_fft+zcc_2*O12_i_alta_fft_t2)./O12_t2;
rst_Vgen2_bornas_fft = A*(O12_Vgen2_bornas_fft)*UB(3);

disp('f.1) Tensiones en bornas del generador 1 (kV y º)');
disp( [ abs(rst_Vgen1_bornas_fft) 180*angle(rst_Vgen1_bornas_fft)/pi ] )
disp('f.2) Tensiones en bornas del generador 2 (kV y º)');
disp( [ abs(rst_Vgen2_bornas_fft) 180*angle(rst_Vgen2_bornas_fft)/pi ] )



