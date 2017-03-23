%Programa final do projeto 1 parte 1B

clc;
clear;

%%Dados de projeto

V0 = 220; %Tens�o rms da fonte
V0_min = 220*0.9; %Tens�o rms m�nima da fonte
V0_max = 220*1.1; %Tens�o rms m�xima da fonte
f = 60; %Frequ�ncia da rede el�trica em Hertz
Vto = 0.85; %Queda de tens�o no diodo
I_FRMS = 150; %Corrente eficaz m�xima que o diodo suporta
I_FAV = 95; %Corrente m�dia m�xima que o diodo suporta
R_diodo = 0.003; %Resist�ncia do diodo
R_alimentador_entrada = 0.064; %Resist�ncia dos cabos que saem da fonte e chegam no conversor (ida e volta)
R_alimentador_saida = 0.065; %Resist�ncia dos cabos que saem do conversor e chegam na carga (ida e volta)
Rjc = 0.55; %Resist�ncia t�rmica jun��o-capsula
Rcd = 0.2; %Resist�ncia t�rmica capsula-dissipador
Tj = 180; %Temperatura da jun��o
Ta = 60; %Temperatura ambiente
m = 2.05; %Rela��o de transforma��o do trafo
n = 0;

w = 2*pi*f; %Frequ�ncia da rede el�trica em rad/s
T = 1/f; %Per�odo da rede el�trica em segundos
Beta = 272; %�ngulo de condu��o do diodo em graus
t_Beta = roundn((Beta*T/360),-5); %Tempo de condu��o do diodo
syms t %Declarando o tempo como vari�vel simb�lica



%%Tens�es
V = @(t)(sqrt(2)*V0*m*sin(w*t));
V_quadrado = @(t)V(t).^2;

V_min = @(t)(sqrt(2)*V0_min*m*sin(w*t));
V_quadrado_min = @(t)V_min(t).^2;

V_max = @(t)(sqrt(2)*V0_max*m*sin(w*t));
V_quadrado_max = @(t) V_max(t).^2;

V_rms_fonte = roundn(sqrt((1/T)*integral(V_quadrado,0,T)),-5);
V_rms_carga = roundn(sqrt((1/T)*integral(V_quadrado,0,t_Beta)),-5);
V_med_fonte = roundn((1/T)*integral(V,0,T),-5);
V_med_carga = roundn((1/T)*integral(V,0,t_Beta),-5);

V_rms_fonte_min = roundn(sqrt((1/T)*integral(V_quadrado_min,0,T)),-5);
V_rms_carga_min = roundn(sqrt((1/T)*integral(V_quadrado_min,0,t_Beta)),-5);
V_med_fonte_min = roundn((1/T)*integral(V_min,0,T),-5);
V_med_carga_min = roundn((1/T)*integral(V_min,0,t_Beta),-5);

V_rms_fonte_max = roundn(sqrt((1/T)*integral(V_quadrado_max,0,T)),-5);
V_rms_carga_max = roundn(sqrt((1/T)*integral(V_quadrado_max,0,t_Beta)),-5);
V_med_fonte_max = roundn((1/T)*integral(V_max,0,T),-5);
V_med_carga_max = roundn((1/T)*integral(V_max,0,t_Beta),-5);

%%Imped�ncia
R_carga_min = roundn(V_med_carga_max/I_FAV,-5); %%Resist�ncia de carga m�nima
R_carga_comercial = roundn(R_carga_min/0.8,-5); %%Resist�ncia comercial utilizada
R_carga_min = 1.15; %%Arredondamento
R = R_alimentador_entrada + R_diodo + R_alimentador_saida + R_carga_min; %Resist�ncia total vista da fonte
L = 0.0116; %Indut�ncia calculada no simulador
X = w*L; %Reat�ncia do circuito
Phi = atan(X/R); %�ngulo de carga em rad
Phi_graus = radtodeg(Phi); %Angulo de carga em graus

%%Correntes
I = @(t) (sqrt(2)*V0*m/sqrt(R^2+X^2))*(sin(w*t-Phi)-(sin(-Phi)*exp(-(t*R/L))));
I_quadrado = @(t) I(t).^2;

I_min = @(t) (sqrt(2)*V0_min*m/sqrt(R^2+X^2))*(sin(w*t-Phi)-(sin(-Phi)*exp(-(t*R/L))));
I_quadrado_min = @(t) I_min(t).^2;

I_max = @(t) (sqrt(2)*V0_max*m/sqrt(R^2+X^2))*(sin(w*t-Phi)-(sin(-Phi)*exp(-(t*R/L))));
I_quadrado_max = @(t) I_max(t).^2;

I_rms_fonte = roundn(sqrt((1/T)*integral(I_quadrado,0,t_Beta)),-5);
I_rms_carga = roundn(sqrt((1/T)*integral(I_quadrado,0,t_Beta)),-5);
I_med_fonte = roundn((1/T)*integral(I,0,t_Beta),-5);
I_med_carga = roundn((1/T)*integral(I,0,t_Beta),-5);

I_rms_fonte_min = roundn(sqrt((1/T)*integral(I_quadrado_min,0,t_Beta)),-5);
I_rms_carga_min = roundn(sqrt((1/T)*integral(I_quadrado_min,0,t_Beta)),-5);
I_med_fonte_min = roundn((1/T)*integral(I_min,0,t_Beta),-5);
I_med_carga_min = roundn((1/T)*integral(I_min,0,t_Beta),-5);

I_rms_fonte_max = roundn(sqrt((1/T)*integral(I_quadrado_max,0,t_Beta)),-5);
I_rms_carga_max = roundn(sqrt((1/T)*integral(I_quadrado_max,0,t_Beta)),-5);
I_med_fonte_max = roundn((1/T)*integral(I_max,0,t_Beta),-5);
I_med_carga_max = roundn((1/T)*integral(I_max,0,t_Beta),-5);

%%Informa��es na fonte/prim�rio
V_rms_primario = V0;
V_rms_primario_max = V0_max;
V_rms_primario_min = V0_min;
V_med_primario = 0;
V_med_primario_max = 0;
V_med_primario_min = 0;

I_rms_primario = I_rms_fonte*m;
I_rms_primario_max = I_rms_fonte_max*m;
I_rms_primario_min = I_rms_fonte_min*m;
I_med_primario = I_med_fonte*m;
I_med_primario_max = I_med_fonte_max*m;
I_med_primario_min = I_med_fonte_min*m;


%%S�rie de Fourier da corrente na carga


IL = @(t) (sqrt(2)*V0_max*m/sqrt(R^2+X^2))*(sin(w*t-Phi)-(sin(-Phi)*exp(-(t*R/L)))); %Equa��o de condu��o do diodo para carga puramente resistiva

a0 = (2/T)*integral(IL,0,t_Beta); %Coeficiente da componente CC
IL_DC = a0/2; %Corrente CC na carga

a = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena os coeficientes a_n
b = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena os coeficientes b_n
c = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena os coeficientes c_n
theta = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena as defasagens theta_n


for n = 1:25
    ILa = @(t) IL(t).*cos(n*w*t);
    ILb = @(t) IL(t).*sin(n*w*t);
    a(n) = roundn((2/T)*integral(ILa,0,t_Beta),-5); %Coeficiente a_n do harm�nico n
    b(n) = roundn((2/T)*integral(ILb,0,t_Beta),-5); %Coeficiente b_n do harm�nico n
    if(b(n) == 0)
        b(n) = 0;
    end
    c(n) = roundn(sqrt((a(n)^2)+(b(n)^2)),-5); %Coeficiente c_n do harm�nico n
    theta(n) = roundn(atan(a(n)/b(n)),-5); %Defasagem theta_n do harm�nico n
    if(a(n) == 0)
        theta(n) = 0;
    end
end

%%S�rie de Fourier da tens�o na carga

%VL = @(t) (sqrt(2)*V0_max*m)*(sin(w*t-Phi)-(sin(-Phi)*exp(-(t*R/L)))); %Equa��o da tens�o na carga
VL = @(t) (sqrt(2)*V0_max*m)*(sin(w*t)); %Equa��o da tens�o na carga

a0_V = (2/T)*integral(VL,0,t_Beta); %Coeficiente da componente CC
VL_DC = a0_V/2; %Tens�o CC na carga

a_V = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena os coeficientes a_n
b_V = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena os coeficientes b_n
c_V = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena os coeficientes c_n
theta_V = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Vetor que armazena as defasagens theta_n


for n = 1:25
    VLa = @(t) VL(t).*cos(n*w*t);
    VLb = @(t) VL(t).*sin(n*w*t);
    a_V(n) = roundn((2/T)*integral(VLa,0,t_Beta),-5); %Coeficiente a_n do harm�nico n
    b_V(n) = roundn((2/T)*integral(VLb,0,t_Beta),-5); %Coeficiente b_n do harm�nico n
    if(b_V(n) == 0)
        b_V(n) = 0;
    end
    c_V(n) = roundn(sqrt((a_V(n)^2)+(b_V(n)^2)),-5); %Coeficiente c_n do harm�nico n
    theta_V(n) = roundn(atan(a_V(n)/b_V(n)),-5); %Defasagem theta_n do harm�nico n
    if(a_V(n) == 0)
        theta_V(n) = 0;
    end
end

%%Pot�ncia ativa
somatorio_potencia = 0;
for n = 1:25
    somatorio_potencia = somatorio_potencia + ((c_V(n)/sqrt(2))*(c(n)/sqrt(2))*cos(theta_V(n)-theta(n)));
end
P_carga = VL_DC*IL_DC + somatorio_potencia;
P_alimentador_entrada = R_alimentador_entrada*(I_rms_fonte_max^2);
P_alimentador_saida = R_alimentador_saida*(I_rms_carga_max^2);
P_diodo = Vto*I_med_fonte + R_diodo*(I_rms_fonte^2);
P_fonte = P_carga + P_alimentador_saida + P_diodo + P_alimentador_entrada;

%%Pot�ncia aparente
S_carga = V_rms_carga_max*I_rms_carga_max;
S_fonte = V_rms_fonte_max*I_rms_fonte_max;

%%Rendimento
n_sistema = P_carga/P_fonte;
n_conversor = (P_carga + P_alimentador_saida)/(P_fonte-P_alimentador_entrada);

%%Fator de pot�ncia
FP = P_fonte/S_fonte;
FP_fundamental = (c_V(1)*c(1)*cos(theta_V(1)-theta(1))/(V_rms_fonte_max*I_rms_fonte_max));


%%THD de corrente

I1 = c(1); %Valor da componente fundamental
Somatorio_corrente = 0; %Somatorio que guardar� a soma dos harm�nicos
for n = 2:25
    Somatorio_corrente = Somatorio_corrente + (c(n)^2);
end

THD_corrente = sqrt(Somatorio_corrente)/I1; %Fator de distor��o harm�nica de corrente

%%Calculando os ham�nicos pares em rela��o ao fundamental

Percentual_harmonico = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
for n = 2:25
    Percentual_harmonico(n) = (c(n)/c(1))*100;
end

%%Calculando o fator de ondula��o de corrente na carga

Somatorio_ripple = 0; %Somat�rio que guardar� os harmonicos ao quadrado
for n = 1:25
    Somatorio_ripple = Somatorio_ripple + ((c(n)/sqrt(2))^2);
end
I_ac = sqrt(Somatorio_ripple);
Ripple_corrente = (I_ac/(a0/2))*100; %Ripple de corrente


%%Calculando o THD de tens�o na carga

V1 = c_V(1); %Valor da componente fundamental
Somatorio_tensao = 0; %Somatorio que guardar� a soma dos harm�nicos
for n = 2:25
    Somatorio_tensao = Somatorio_tensao + (c_V(n)^2);
end

THD_tensao = sqrt(Somatorio_tensao)/V1; %Fator de distor��o harm�nica de tens�o

%%Calculando os ham�nicos pares em rela��o ao fundamental

Percentual_harmonico_V = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
for n = 2:25
    Percentual_harmonico_V(n) = (c_V(n)/c_V(1))*100;
end

%%Calculando o fator de ondula��o de tens�o na carga

Somatorio_ripple_V = 0; %Somat�rio que guardar� os harmonicos ao quadrado
for n = 1:25
    Somatorio_ripple_V = Somatorio_ripple_V + ((c_V(n)/sqrt(2))^2);
end
V_ac = sqrt(Somatorio_ripple_V);
Ripple_tensao = (V_ac/(a0_V/2))*100; %Ripple de corrente

%%Calculando o dissipador
P_diodo_nominal = Vto*I_med_carga + R_diodo*(I_rms_carga^2);
Rja = (Tj - Ta)/P_diodo_nominal;
Rda = Rja - Rjc - Rcd;
T_atingida = P_diodo_nominal*(Rcd + Rda) + Ta;

%%Gerando os gr�ficos

tempo = 0:0.00001:T; %Amostrando o eixo do tempo

 %Gr�fico da tens�o na fonte
figure(1)
hold on
for n = 1:25
    plot(tempo,(c_V(n).*sin(n*w*tempo+theta_V(n))));
end
title('25 primeiros harm�nicos da tens�o na carga');
xlabel('Tempo(s)');
ylabel('Tens�o(V)');
hold off

%%Gerando os gr�ficos dos harm�nicos de corrente na carga
figure(2)
hold on
for n = 1:25
    plot(tempo,(c(n).*sin(n*w*tempo+theta(n))));
end
title('25 primeiros harm�nicos da corrente do circuito');
xlabel('Tempo(s)');
ylabel('Corrente(A)');
hold off

%%Gerando arquivo de sa�da txt com as informa��es

fid = fopen('Trabalho 1B.txt','wt');
fprintf(fid,'             TRABALHO 1B DE ELETR�NICA DE POT�NCIA\n');
fprintf(fid,'          RETIFICADOR DE MEIA ONDA COM CARGA RESISTIVA\n\n');
fprintf(fid,'             CARACTERIZA��O DA FONTE \n\n');
fprintf(fid,'Par�metro     M�nimo     Nominal     M�ximo     Unidade \n');
fprintf(fid,'  Vrms       %.4f    %.4f    %.4f    Vrms \n',V_rms_primario_min,V_rms_primario,V_rms_primario_max);
fprintf(fid,'  Vmed       %.4f      %.4f      %.4f       V \n',V_med_primario_min,V_med_primario,V_med_primario_max);
fprintf(fid,'  Irms       %.4f    %.4f    %.4f    Arms \n',I_rms_primario_min,I_rms_primario,I_rms_primario_max);
fprintf(fid,'  Imed       %.4f    %.4f    %.4f     A \n\n',I_med_primario_min,I_med_primario,I_med_primario_max);
fprintf(fid,'             CARACTERIZA��O DO SECUND�RIO \n\n');
fprintf(fid,'Par�metro     M�nimo     Nominal     M�ximo     Unidade \n');
fprintf(fid,'  Vrms       %.4f    %.4f    %.4f    Vrms \n',V_rms_fonte_min,V_rms_fonte,V_rms_fonte_max);
fprintf(fid,'  Vmed      %.4f     %.4f     %.4f       V \n',V_med_fonte_min,V_med_fonte,V_med_fonte_max);
fprintf(fid,'  Irms       %.4f    %.4f    %.4f    Arms \n',I_rms_fonte_min,I_rms_fonte,I_rms_fonte_max);
fprintf(fid,'  Imed       %.4f     %.4f     %.4f      A \n\n',I_med_fonte_min,I_med_fonte,I_med_fonte_max);

fprintf(fid,'               CARACTERIZA��O DA CARGA \n\n');
fprintf(fid,'Par�metro     M�nimo     Nominal     M�ximo     Unidade \n');
fprintf(fid,'  Vrms       %.4f    %.4f   %.4f     Vrms \n',V_rms_carga_min,V_rms_carga,V_rms_carga_max);
fprintf(fid,'  Vmed       %.4f     %.4f    %.4f      V \n',V_med_carga_min,V_med_carga,V_med_carga_max);
fprintf(fid,'  Irms       %.4f    %.4f   %.4f     Arms \n',I_rms_carga_min,I_rms_carga,I_rms_carga_max);
fprintf(fid,'  Imed       %.4f     %.4f    %.4f       A \n\n',I_med_carga_min,I_med_carga,I_med_carga_max);


fprintf(fid,'Pot�ncia ativa m�xima na fonte = %.4f W\n',P_fonte);
fprintf(fid,'Pot�ncia ativa m�xima na carga = %.4f W\n',P_carga);
fprintf(fid,'Pot�ncia aparente m�xima na fonte = %.4f VA\n',S_fonte);
fprintf(fid,'Pot�ncia aparente m�xima na carga = %.4f VA\n\n',S_carga);
fprintf(fid,'Rendimento do sistema = %.4f\n',n_sistema);
fprintf(fid,'Rendimento do conversor = %.4f\n\n',n_conversor);
fprintf(fid,'Fator de pot�ncia = %.4f\n',FP);
fprintf(fid,'Fator de pot�ncia em rela��o ao primeiro harm�nico = %.4f\n\n',FP_fundamental);
fprintf(fid,'THD de tens�o = %.4f\n',THD_tensao);
fprintf(fid,'THD de corrente = %.4f\n\n',THD_corrente);
fprintf(fid,'Componente CC da tensao = a0/2 = %.4f V\n',a0_V/2);
fprintf(fid,'Componente CC da corrente na carga = a0/2 = %.4f A\n\n',a0/2);
fprintf(fid,'Ripple de corrente = %.4f\n',Ripple_corrente);
fprintf(fid,'Ripple de tens�o = %.4f\n\n',Ripple_tensao);
fprintf(fid,'Pot�ncia nominal dissipada pelo diodo = %.4f W\n',P_diodo_nominal);
fprintf(fid,'Rja = %.4f �C/W\n',Rja);
fprintf(fid,'Rda = %.4f �C/W\n\n\n',Rda);

fprintf(fid,'             S�RIE DE FOURIER DA CORRENTE\n\n');
fprintf(fid,'                    COEFICIENTES\n\n');
fprintf(fid,'  n          a               b             c           theta\n');
fprintf(fid,'  1       %.4f        %.4f       %.4f     %.4f\n',a(1),b(1),c(1),theta(1));
fprintf(fid,'  2        %.4f        %.4f       %.4f      %.4f\n',a(2),b(2),c(2),theta(2));
fprintf(fid,'  3       %.4f         %.4f        %.4f        %.4f\n',a(3),b(3),c(3),theta(3));
fprintf(fid,'  4       %.4f         %.4f        %.4f        %.4f\n',a(4),b(4),c(4),theta(4));
fprintf(fid,'  5        %.4f          %.4f        %.4f        %.4f\n',a(5),b(5),c(5),theta(5));
fprintf(fid,'  6        %.4f         %.4f        %.4f       %.4f\n',a(6),b(6),c(6),theta(6));
fprintf(fid,'  7       %.4f         %.4f        %.4f        %.4f\n',a(7),b(7),c(7),theta(7));
fprintf(fid,'  8       %.4f          %.4f        %.4f       %.4f\n',a(8),b(8),c(8),theta(8));
fprintf(fid,'  9        %.4f          %.4f        %.4f        %.4f\n',a(9),b(9),c(9),theta(9));
fprintf(fid,'  10       %.4f         %.4f        %.4f       %.4f\n',a(10),b(10),c(10),theta(10));
fprintf(fid,'  11      %.4f         %.4f        %.4f        %.4f\n',a(11),b(11),c(11),theta(11));
fprintf(fid,'  12      %.4f          %.4f        %.4f       %.4f\n',a(12),b(12),c(12),theta(12));
fprintf(fid,'  13       %.4f          %.4f        %.4f        %.4f\n',a(13),b(13),c(13),theta(13));
fprintf(fid,'  14       %.4f         %.4f        %.4f       %.4f\n',a(14),b(14),c(14),theta(14));
fprintf(fid,'  15      %.4f         %.4f        %.4f        %.4f\n',a(15),b(15),c(15),theta(15));
fprintf(fid,'  16      %.4f          %.4f        %.4f       %.4f\n',a(16),b(16),c(16),theta(16));
fprintf(fid,'  17       %.4f          %.4f        %.4f        %.4f\n',a(17),b(17),c(17),theta(17));
fprintf(fid,'  18       %.4f         %.4f        %.4f       %.4f\n',a(18),b(18),c(18),theta(18));
fprintf(fid,'  19      %.4f         %.4f        %.4f        %.4f\n',a(19),b(19),c(19),theta(19));
fprintf(fid,'  20      %.4f          %.4f        %.4f       %.4f\n',a(20),b(20),c(20),theta(20));
fprintf(fid,'  21       %.4f          %.4f        %.4f        %.4f\n',a(21),b(21),c(21),theta(21));
fprintf(fid,'  22       %.4f         %.4f        %.4f       %.4f\n',a(22),b(22),c(22),theta(22));
fprintf(fid,'  23      %.4f         %.4f        %.4f        %.4f\n',a(23),b(23),c(23),theta(23));
fprintf(fid,'  24      %.4f          %.4f        %.4f       %.4f\n',a(24),b(24),c(24),theta(24));
fprintf(fid,'  25       %.4f          %.4f        %.4f        %.4f\n\n\n',a(25),b(25),c(25),theta(25));

fprintf(fid,'PERCENTUAL DOS HARM�NICOS PARES EM RELA��O AO FUNDAMENTAL\n\n');
fprintf(fid,'        Harm�nico           Rela��o com o fundamental\n');
fprintf(fid,'            2                        %.4f\n',Percentual_harmonico(2));
fprintf(fid,'            4                        %.4f\n',Percentual_harmonico(4));
fprintf(fid,'            6                        %.4f\n',Percentual_harmonico(6));
fprintf(fid,'            8                        %.4f\n',Percentual_harmonico(8));
fprintf(fid,'            10                       %.4f\n',Percentual_harmonico(10));
fprintf(fid,'            12                       %.4f\n',Percentual_harmonico(12));
fprintf(fid,'            14                       %.4f\n',Percentual_harmonico(14));
fprintf(fid,'            16                       %.4f\n',Percentual_harmonico(16));
fprintf(fid,'            18                       %.4f\n',Percentual_harmonico(18));
fprintf(fid,'            20                       %.4f\n',Percentual_harmonico(20));
fprintf(fid,'            22                       %.4f\n',Percentual_harmonico(22));
fprintf(fid,'            24                       %.4f\n\n\n',Percentual_harmonico(24));

fprintf(fid,'          S�RIE DE FOURIER DA TENS�O NA CARGA\n\n');
fprintf(fid,'                    COEFICIENTES\n\n');
fprintf(fid,'  n          a               b             c           theta\n');
fprintf(fid,'  1       %.4f         %.4f      %.4f      %.4f\n',a_V(1),b_V(1),c_V(1),theta_V(1));
fprintf(fid,'  2      %.4f         %.4f      %.4f      %.4f\n',a_V(2),b_V(2),c_V(2),theta_V(2));
fprintf(fid,'  3      %.4f        %.4f        %.4f      %.4f\n',a_V(3),b_V(3),c_V(3),theta_V(3));
fprintf(fid,'  4      %.4f          %.4f       %.4f      %.4f\n',a_V(4),b_V(4),c_V(4),theta_V(4));
fprintf(fid,'  5       %.4f          %.4f        %.4f       %.4f\n',a_V(5),b_V(5),c_V(5),theta_V(5));
fprintf(fid,'  6       %.4f          %.4f       %.4f      %.4f\n',a_V(6),b_V(6),c_V(6),theta_V(6));
fprintf(fid,'  7      %.4f         %.4f        %.4f       %.4f\n',a_V(7),b_V(7),c_V(7),theta_V(7));
fprintf(fid,'  8      %.4f          %.4f       %.4f      %.4f\n',a_V(8),b_V(8),c_V(8),theta_V(8));
fprintf(fid,'  9       %.4f          %.4f        %.4f       %.4f\n',a_V(9),b_V(9),c_V(9),theta_V(9));
fprintf(fid,'  10      %.4f          %.4f       %.4f      %.4f\n',a_V(10),b_V(10),c_V(10),theta_V(10));
fprintf(fid,'  11     %.4f         %.4f        %.4f       %.4f\n',a_V(11),b_V(11),c_V(11),theta_V(11));
fprintf(fid,'  12     %.4f           %.4f       %.4f      %.4f\n',a_V(12),b_V(12),c_V(12),theta_V(12));
fprintf(fid,'  13      %.4f          %.4f        %.4f       %.4f\n',a_V(13),b_V(13),c_V(13),theta_V(13));
fprintf(fid,'  14      %.4f          %.4f       %.4f      %.4f\n',a_V(14),b_V(14),c_V(14),theta_V(14));
fprintf(fid,'  15     %.4f         %.4f        %.4f       %.4f\n',a_V(15),b_V(15),c_V(15),theta_V(15));
fprintf(fid,'  16     %.4f           %.4f       %.4f      %.4f\n',a_V(16),b_V(16),c_V(16),theta_V(16));
fprintf(fid,'  17      %.4f          %.4f        %.4f       %.4f\n',a_V(17),b_V(17),c_V(17),theta_V(17));
fprintf(fid,'  18      %.4f          %.4f       %.4f      %.4f\n',a_V(18),b_V(18),c_V(18),theta_V(18));
fprintf(fid,'  19     %.4f         %.4f        %.4f       %.4f\n',a_V(19),b_V(19),c_V(19),theta_V(19));
fprintf(fid,'  20     %.4f           %.4f        %.4f      %.4f\n',a_V(20),b_V(20),c_V(20),theta_V(20));
fprintf(fid,'  21      %.4f           %.4f        %.4f       %.4f\n',a_V(21),b_V(21),c_V(21),theta_V(21));
fprintf(fid,'  22      %.4f          %.4f        %.4f       %.4f\n',a_V(22),b_V(22),c_V(22),theta_V(22));
fprintf(fid,'  23     %.4f          %.4f        %.4f       %.4f\n',a_V(23),b_V(23),c_V(23),theta_V(23));
fprintf(fid,'  24     %.4f           %.4f        %.4f       %.4f\n',a_V(24),b_V(24),c_V(24),theta_V(24));
fprintf(fid,'  25      %.4f           %.4f        %.4f        %.4f\n\n\n',a_V(25),b_V(25),c_V(25),theta_V(25));

fprintf(fid,'PERCENTUAL DOS HARM�NICOS PARES EM RELA��O AO FUNDAMENTAL\n\n');
fprintf(fid,'        Harm�nico           Rela��o com o fundamental\n');
fprintf(fid,'            2                        %.4f\n',Percentual_harmonico_V(2));
fprintf(fid,'            4                        %.4f\n',Percentual_harmonico_V(4));
fprintf(fid,'            6                        %.4f\n',Percentual_harmonico_V(6));
fprintf(fid,'            8                        %.4f\n',Percentual_harmonico_V(8));
fprintf(fid,'            10                       %.4f\n',Percentual_harmonico_V(10));
fprintf(fid,'            12                       %.4f\n',Percentual_harmonico_V(12));
fprintf(fid,'            14                       %.4f\n',Percentual_harmonico_V(14));
fprintf(fid,'            16                       %.4f\n',Percentual_harmonico_V(16));
fprintf(fid,'            18                       %.4f\n',Percentual_harmonico_V(18));
fprintf(fid,'            20                       %.4f\n',Percentual_harmonico_V(20));
fprintf(fid,'            22                       %.4f\n',Percentual_harmonico_V(22));
fprintf(fid,'            24                       %.4f\n',Percentual_harmonico_V(24));