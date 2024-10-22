clc; close all; clear 
%% Simulación de electrodeposicioón a partir de automatas celculares 
%% Espacio
x=100;   y=100;
interaciones=200;
Space=zeros(x+1,y+1,interaciones);
Space_o=zeros(x+1,y+1);
n_Cu=100; n_CN=100;

%% Posiciones Posibles 
k=0;   Posicion=zeros(x*y,2);  
for i=1:x
    for j=1:y
        k=k+1;
        Posicion(k,:)=[i,j];
    end
end 

%% Vecindario
veci=cell(x,y);
for f=1:length(Posicion)
    i=Posicion(f,1);  j=Posicion(f,2);
    switch i && j
        case i==1 && j~=y && j~=1                        % Vecindario frontera inferior
            veci{i,j}=[i+1,j;i,j+1;x,j;i,j-1];
        case j==1 && i~=x && i~=1                        % Vecindario frontera izquierda
            veci{i,j}=[i+1,j;i,j+1;i-1,j;i,y];
        case i==x && j~=1 && j~=y                        % Vecindario frontera superior
            veci{i,j}=[1,j;i,j+1;i-1,j;i,j-1];
        case j==y && i~=1 && i~=x                        % Vecindario frontera derecha
            veci{i,j}=[i+1,j;i,1;i-1,j;i,j-1];
        case j==1 && i==1                                % Vecindario punta interior-izquierda
            veci{i,j}=[i+1,j;i,j+1;x,j;i,y];
        case j==y && i==x                                % Vecindario punta superior-derecha
            veci{i,j}=[1,j;i,1;i-1,j;i,j-1];
        case j==y && i==1                                % Vecindario punta inferior-derecha
            veci{i,j}=[i+1,j;i,1;x,j;i,j-1];
        case j==1 && i==x                                % Vecindario punta superior-izquierda
            veci{i,j}=[1,j;i,j+1;i-1,j;i,y];
        otherwise
            veci{i,j}=[i+1,j;i,j+1;i-1,j;i,j-1];         % Vecindario de puntos internos
    end
end

%% Condiciones iniciales 
P_ac=[0 0.25 0.5 0.75 1];       % Probabilidad acumulada de realizar un moviento en cualquier direacción
mo=cell(x,y);
pos=zeros(n_Cu+n_CN,2,interaciones);
Q=zeros(x,y);

for i=1:n_Cu
    while 1
        p=[randi(x,1,1),randi(y,1,1)];
        if Space_o(p(1),p(2))~=1&&Space_o(p(1),p(2))~=2
            Space_o(p(1),p(2))=1;
            pos(i,:,1)=[p(1) p(2)];
            break;
        end
    end
end
for i=1:n_CN
    while 1
        p=[randi(x,1,1),randi(y,1,1)];
        [m1,m2,r]=m(p(1),p(2),P_ac,veci);
        if Space_o(p(1),p(2))==0&&Space_o(m1,m2)==0
            Space_o(p(1),p(2))=2;
             pos(n_Cu+i,:,1)=[p(1) p(2)];
            break;
        end
    end
end

Space(:,:,1)=Space_o;
%% Movimiento
P=0.5;                             % Probabilidad de no realizar un movimiento
D=zeros(x+1,y+1); D1=zeros(x+1,y+1); R=zeros(x+1,y+1);
N=zeros(interaciones,6); xt=zeros(interaciones,6);  pos2=zeros(n_Cu+n_CN,2);

for k=1:interaciones
    pos1=Pos(pos(:,:,k),n_Cu,n_CN);
    pos3=pos1(:,1);  pos4=pos1(:,2);
    pos3(pos3==0)=[];  pos4(pos4==0)=[];
    pos1=[pos3,pos4];
    for h=1:length(pos1(:,1))
        i=pos1(h,1);   j=pos1(h,2);
        %% Ion Cobre
        if Space(i,j,k)==1
                [D,R,pos2]=Reaccion(i,j,k,h,R,D,Space,pos2,veci,2);
            if R(i,j)~=1
                [D,pos2]=Movimiento(i,j,k,h,x,y,pos2,D,P,pos1,Space,veci,1);
            end
            %% Ion Cianuro
        elseif Space(i,j,k)==2
            if R(i,j)~=1
                [D,pos2]=Movimiento(i,j,k,h,x,y,pos2,D,P,pos1,Space,veci,2);
            end
            %% Ion Cianuro de Cobre I
            
        elseif Space(i,j,k)==3
            [D,R,pos2]=Reaccion(i,j,k,h,R,D,Space,pos2,veci,1);
            if R(i,j)~=1
                [D,pos2]=Movimiento(i,j,k,h,x,y,pos2,D,P,pos1,Space,veci,3);
            end
            %% Ion Cianuro de Cobre II
        elseif Space(i,j,k)==4
            [D,R,pos2]=Reaccion(i,j,k,h,R,D,Space,pos2,veci,1);
            if R(i,j)~=1
                [D,pos2]=Movimiento(i,j,k,h,x,y,pos2,D,P,pos1,Space,veci,4);
            end
            %% Ion Cianuro de Cobre III
        elseif Space(i,j,k)==5
            [D,R,pos2]=Reaccion(i,j,k,h,R,D,Space,pos2,veci,1);
           if R(i,j)~=1
                [D,pos2]=Movimiento(i,j,k,h,x,y,pos2,D,P,pos1,Space,veci,5);
           end
            %% Ion Cianuro de Cobre IV
        elseif Space(i,j,k)==6
           if R(i,j)~=1
                [D,pos2]=Movimiento(i,j,k,h,x,y,pos2,D,P,pos1,Space,veci,6);
           end
        end
    end
    %% Se guardan los cambios en cada iteración y se borra la matirz de cambio
    R=zeros(x+1,y+1);
    Space(:,:,k+1)=Space(:,:,k)+D;
    D1(:,:,k)=D(:,:);
    D(:,:)=0;
    %% Concentraciones
    N(k,:)=[length(find(Space(:,:,k)==1)),length(find(Space(:,:,k)==2)),length(find(Space(:,:,k)==3)),length(find(Space(:,:,k)==4)),length(find(Space(:,:,k)==5)),length(find(Space(:,:,k)==6))];
    suma_N=sum(N(k,:));
    xt(k,:)=N(k,:)/suma_N;
    pos(:,:,k+1)=pos2;
    pos2=zeros(n_Cu+n_CN,2);
end

%% Simulación
for l=1:interaciones 
    data=Space(:,:,l);
    hAxes=gca;
    imagesc(hAxes,data);
    colormap(hAxes, [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1; 1 1 0] )
    colorbar
    pause(0)
end

%% Concentraciones 
figure (1)
plot(xt);

%% Funciones 
%%  Direccion del movimiento
function [m1,m2,r]=m(i,j,P_ac,veci)
e=rand();                               % 0<e<0.25 (vecino de arriba) 0.25<e<0.5 (vecino derecho)
r=find(P_ac>e);   r=r(1)-1;             % 0.5<e<0.75 (vecino de abajo) 0.75<e<1 (vecino izquierdo)
mol=veci{i,j}(r,:);
m1=mol(1);  m2=mol(2);
end

%% Posiciones aleatorias al momento del decidir el movimiento
function [Posicion1]=Pos(Posicion,n_Cu,n_CN)
Posicion1=zeros(n_Cu+n_CN,2);
for h=1:length(Posicion)
    e=randi(n_Cu+n_CN-h+1);
    while e==0
        e=randi(n_Cu+n_CN-h+1);
    end
    Posicion1(h,:)=Posicion(e,:);
    Posicion(e,:)=[];
end
end

%% Función de movimiento
function [D,pos]=Movimiento(i,j,k,h,x,y,pos,D,P,pos1,Space,veci,n)
vecin=veci{i,j};  o=0;
r1i=vecin(1,1); r1j=vecin(1,2); r2i=vecin(2,1); r2j=vecin(2,2);
r3i=vecin(3,1); r3j=vecin(3,2); r4i=vecin(4,1); r4j=vecin(4,2);

%% Se pregunta si el automata realizará el movimiento
e=rand();
if e<P
    pos(h,:)=[i j];
    o=1;
end
while o==0
    %% Se elige otra dirrección y se guardan los cambios
    P_ac=interacion(i,j,k,pos1,x,y,Space);
    [m1,m2,~]=m(i,j,P_ac,veci);
    %% Se guarda el moviemiento en la matiz de cambio
    D(m1,m2)=D(m1,m2)+n;    D(i,j)=D(i,j)-n;           % Matiz de cambio "D"
    if D(m1,m2)==n && Space(m1,m2,k)==0
        pos(h,:)=[m1 m2];
        break;
    end
    %% De estar ocupada la posicion, se borra el movimiento y se repite el proceso
    D(m1,m2)=D(m1,m2)-n;    D(i,j)=D(i,j)+n;           % Matiz de cambio "D"
    if Space(r1i,r1j,k)==0&&Space(r2i,r2j,k)==0&&Space(r3i,r3j,k)==0&&Space(r4i,r4j,k)==0
        pos(h,:)=[i j];
        break;
    end
end
end

%% Funcion Reaccion
function [D,R,pos]=Reaccion(i,j,k,h,R,D,Space,pos,veci,n)
%% Reacción
vec_ij=veci{i,j};
if Space(vec_ij(1,1),vec_ij(1,2),k)==2&&D(vec_ij(1,1),vec_ij(1,2))~=-2
    R(i,j)=1; R(vec_ij(1,1),vec_ij(1,2))=1;
    D(i,j)=n;   D(vec_ij(1,1),vec_ij(1,2))=-2;
    pos(h,:)=[i;j];   
elseif Space(vec_ij(2,1),vec_ij(2,2),k)==2&&D(vec_ij(2,1),vec_ij(2,2))~=-2
    R(i,j)=1; R(vec_ij(2,1),vec_ij(2,2))=1;
    D(i,j)=n;   D(vec_ij(2,1),vec_ij(2,2))=-2;
    pos(h,:)=[i;j];  
elseif Space(vec_ij(3,1),vec_ij(3,2),k)==2&&D(vec_ij(3,1),vec_ij(3,2))~=-2
    R(i,j)=1; R(vec_ij(3,1),vec_ij(3,2))=1;
    D(i,j)=n;   D(vec_ij(3,1),vec_ij(3,2))=-2;
    pos(h,:)=[i;j];  
elseif Space(vec_ij(4,1),vec_ij(4,2),k)==2&&D(vec_ij(4,1),vec_ij(4,2))~=-2
    R(i,j)=1; R(vec_ij(4,1),vec_ij(4,2))=1;
    D(i,j)=n;   D(vec_ij(4,1),vec_ij(4,2))=-2;
    pos(h,:)=[i;j];
end
end

%% Calculo de las fuerzas Intermoleculares
function [P_ac]=interacion(i,j,k,pos1,x,y,Space)
m_ij=Space(i,j,k);
Q_ij=round(-0.075*m_ij^5+1.4167*m_ij^4-10.125*m_ij^3+33.583*m_ij^2-50.8*m_ij+27);
F1=0; F2=0; F3=0; F4=0;
if Q_ij~=0
    l=[x,y]; K=8987551788;      
    for d=1:length(pos1(:,1))
        rk=pos1(d,:);
        m_r=Space(rk(1),rk(2),k);
        Q_r=round(-0.075*m_r^5+1.4167*m_r^4-10.125*m_r^3+33.583*m_r^2-50.8*m_r+27);
        if Q_ij*Q_r>0
            r=[i,j]-rk;   
        else
            r=rk-[i,j];
        end
        if r(1)==0 && r(2)==0
        else 
            %% Condicion de Imagen Minima
            for h=1:2
                if r(h)>l(h)/2
                    r(h)=r(h)-l(h);
                elseif r(h)<=-l(h)/2
                    r(h)=r(h)+l(h);
                end
            end
            rd=sqrt(r(2)^2+r(1)^2);
            u=r/rd;            % [i,j] Vector unitario
            F=K*(abs(Q_ij)*abs(Q_r))/rd^2;
            if u(2)>0
                F2=F2+abs(F*u(2));
            else
                F4=F4+abs(F*u(2));
            end
            if u(1)<0
                F3=F3+abs(F*u(1));
            else
                F1= F1+abs(F*u(1));
            end
        end
    end
    Ft=F1+F2+F3+F4;
    P1=F1/Ft;  P2=F2/Ft; P3=F3/Ft; P4=F4/Ft;
    P_ac=[0,P1,P1+P2,P1+P2+P3,P1+P2+P3+P4];
else
    P_ac=[0 0.25 0.5 0.75 1];
end
end
