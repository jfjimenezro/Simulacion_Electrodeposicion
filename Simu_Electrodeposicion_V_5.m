clc; close all; clear
%% Simulación de electrodeposicioón a partir de automatas celculares
%% Espacio
x=90;   y=90;  z=90;
interaciones=500;
Space=zeros(x+1,y+1,z+1,interaciones);
Space_o=zeros(x+1,y+1,z+1);
n_Cu=111; n_CN=1000;

%% Posiciones Posibles
k=0;   Posicion=zeros(x*y*z,3);
for i=1:x
    for j=1:y
        for l=1:z
            k=k+1;
            Posicion(k,:)=[i,j,l];
        end
    end
end

%% Vecindario
veci=cell(x,y,z);
for f=1:length(Posicion)
    i=Posicion(f,1);  j=Posicion(f,2); l=Posicion(f,3);
    switch i && j && l
        case i==1 && j==1 && l==1
            veci{i,j,l}=[i+1,j,l;i,j+1,l;x,j,l;i,y,l;i,j,z;i,j,l+1];
        case i==x && j==1 && l==1
            veci{i,j,l}=[1,j,l;i,j+1,l;i-1,j,l;i,y,l;i,l,z;i,j,l+1];
        case i==1 && j==y && l==1
            veci{i,j,l}=[i+1,j,l;i,1,l;x,j,l;i,j-1,l;i,j,z;i,j,l+1];
        case i==x && j==y && l==1
            veci{i,j,l}=[1,j,l;i,1,l;i-1,j,l;i,j-1,l;i,j,z;i,j,l+1];
        case i==1 && j==1 && l==z
            veci{i,j,l}=[i+1,j,l;i,j+1,l;x,j,l;i,y,l;i,j,l-1;i,j,1];
        case i==x && j==1 && l==z
            veci{i,j,l}=[1,j,l;i,j+1,l;i-1,j,l;i,y,l;i,j,l-1;i,j,1];
        case i==1 && j==y && l==z
            veci{i,j,l}=[i+1,j,l;i,1,l;x,j,l;i,j-1,l;i,j,l-1;i,j,1];
        case i==x && j==y && l==z
            veci{i,j,l}=[1,j,l;i,1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,1];
        case i~=x && i~=1 && j==1 && l==1
            veci{i,j,l}=[i+1,j,l;i,j+1,l;i-1,j,l;i,y,l;i,j,z;i,j,l+1];
        case i==1 && j~=1 && j~=y && l==1
            veci{i,j,l}=[i+1,j,l;i,j+1,l;x,j,l;i,j-1,l;i,j,z;i,j,l+1];
        case i~=1 && i~=x && j==y && l==1
            veci{i,j,l}=[i+1,j,l;i,1,l;i-1,j,l;i,j-1,l;i,j,z;i,j,l+1];
        case i==x && j~=1 && j~=y && l==1
            veci{i,j,l}=[1,j,l;i,j+1,l;i-1,j,l;i,j-1,l;i,j,z;i,j,l+1];
        case i~=1 && i~=x && j==1 && l==z
            veci{i,j,l}=[i+1,j,l;i,j+1,l;i-1,j,l;i,y,l;i,j,l-1;i,j,1];
        case i==1 && j~=1 && j~=y && l==z
            veci{i,j,l}=[i+1,j,l;i,j+1,l;x,j,l;i,j-1,l;i,j,l-1;i,j,1];
        case i~=1 && i~=x && j==y && l==z
            veci{i,j,l}=[i+1,j,l;i,1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,1];
        case i==x && j~=1 && j~=y && l==z
            veci{i,j,l}=[1,j,l;i,j+1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,1];
        case i==1 && j==1 && l~=1 && l~=z
            veci{i,j,l}=[i+1,j,l;i,j+1,l;x,j,l;i,y,l;i,j,l-1;i,j,l+1];
        case i==1 && j==y && l~=1 && l~=z
            veci{i,j,l}=[i+1,j,l;i,1,l;x,j,l;i,j-1,l;i,j,l-1;i,j,l+1];
        case i==x && j==1 && l~=1 && l~=z
            veci{i,j,l}=[1,j,l;i,j+1,l;i-1,j,l;i,y,l;i,l,l-1;i,j,l+1];
        case i==x && j==y && l~=1 && l~=z
            veci{i,j,l}=[1,j,l;i,1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,l+1];
        case i~=x && i~=1 && j~=x && j~=1 && l==1
            veci{i,j,l}=[i+1,j,l;i,j+1,l;i-1,j,l;i,j-1,l;i,j,z;i,j,l+1];
        case i~=x && i~=1 && j~=x && j~=1 && l==z
            veci{i,j,l}=[i+1,j,l;i,j+1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,1];
        case i==1 && j~=x && j~=1 && l~=x && l~=1
            veci{i,j,l}=[i+1,j,l;i,j+1,l;x,j,l;i,j-1,l;i,j,l-1;i,j,l+1];
        case i==x && j~=x && j~=1 && l~=x && l~=1
            veci{i,j,l}=[1,j,l;i,j+1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,l+1];
        case i~=x && i~=1 && j==1 && l~=x && l~=1
            veci{i,j,l}=[i+1,j,l;i,j+1,l;i-1,j,l;i,y,l;i,j,l-1;i,j,l+1];
        case i~=x && i~=1 && j==y && l~=x && l~=1
            veci{i,j,l}=[i+1,j,l;i,1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,l+1];
        otherwise
            veci{i,j,l}=[i+1,j,l;i,j+1,l;i-1,j,l;i,j-1,l;i,j,l-1;i,j,l+1];         % Vecindario de puntos internos
    end
end

%% Condiciones iniciales
% organizar las posiciones aleatorias de los automatas 
P_ac=[0 0.25 0.5 0.75 1];       % Probabilidad acumulada de realizar un moviento en cualquier dirección
pos=zeros(n_Cu+n_CN,4,interaciones);

for i=1:n_Cu
    while 1
        p=[randi(x,1,1),randi(y,1,1),randi(z,1,1)];% vector posicion
        if Space_o(p(1),p(2),p(3))~=1&&Space_o(p(1),p(2),p(3))~=2
            Space_o(p(1),p(2),p(3))=1;
            pos(i,:,1)=[p(1) p(2) p(3) 1];
            break;
        end
    end
end
for i=1:n_CN
    while 1
        p=[randi(x,1,1),randi(y,1,1),randi(z,1,1)];
        if Space_o(p(1),p(2),p(3))==0
            Space_o(p(1),p(2),p(3))=2;
            pos(n_Cu+i,:,1)=[p(1) p(2) p(3) 2];
            break;
        end
    end
end

Space(:,:,:,1)=Space_o;

% figure (1)
% plot3(pos(:,1),pos(:,2),pos(:,3),"gs");
% grid on
%% Movimiento
P=0.5;                             % Probabilidad de no realizar un movimiento
D=zeros(x+1,y+1,z+1); D1=zeros(x+1,y+1,z+1); R=zeros(x+1,y+1,z+1);
N=zeros(interaciones,6); xt=zeros(interaciones,6);  pos2=zeros(n_Cu+n_CN,4);

for k=1:interaciones
    pos1=Pos(pos(:,:,k),n_Cu,n_CN);
    pos3=pos1(:,1);  pos4=pos1(:,2);  pos5=pos1(:,3); pos6=pos1(:,4);
    pos3(pos3==0)=[];  pos4(pos4==0)=[]; pos5(pos5==0)=[];  pos6(pos6==0)=[];
    pos1=[pos3,pos4,pos5,pos6];
    for h=1:length(pos1(:,1))
        i=pos1(h,1);   j=pos1(h,2); l=pos1(h,3); T=pos1(h,4);
        %% Ion Cobre
        if Space(i,j,l,k)==1
            [D,R,pos2]=Reaccion(i,j,l,T,k,h,R,D,Space,pos2,veci,2);
            if R(i,j,l)~=1
                [D,pos2]=Movimiento(i,j,l,T,k,h,x,y,z,pos2,D,P,pos1,Space,veci,1);
            end
            %% Ion Cianuro
        elseif Space(i,j,l,k)==2
            if R(i,j,l)~=1
                [D,pos2]=Movimiento(i,j,l,T,k,h,x,y,z,pos2,D,P,pos1,Space,veci,2);
            end
            %% Ion Cianuro de Cobre I
        elseif Space(i,j,l,k)==3
            [D,R,pos2]=Reaccion(i,j,l,T,k,h,R,D,Space,pos2,veci,1);
            if R(i,j,l)~=1
                [D,pos2]=Movimiento(i,j,l,T,k,h,x,y,z,pos2,D,P,pos1,Space,veci,3);
            end
            %% Ion Cianuro de Cobre II
        elseif Space(i,j,l,k)==4
            [D,R,pos2]=Reaccion(i,j,l,T,k,h,R,D,Space,pos2,veci,1);
            if R(i,j,l)~=1
                [D,pos2]=Movimiento(i,j,l,T,k,h,x,y,z,pos2,D,P,pos1,Space,veci,4);
            end
            %% Ion Cianuro de Cobre III
        elseif Space(i,j,l,k)==5
            [D,R,pos2]=Reaccion(i,j,l,T,k,h,R,D,Space,pos2,veci,1);
            if R(i,j,l)~=1
                [D,pos2]=Movimiento(i,j,l,T,k,h,x,y,z,pos2,D,P,pos1,Space,veci,5);
            end
            %% Ion Cianuro de Cobre IV
        elseif Space(i,j,l,k)==6
            if R(i,j,l)~=1
                [D,pos2]=Movimiento(i,j,l,T,k,h,x,y,z,pos2,D,P,pos1,Space,veci,6);
            end
        end
    end
    %% Se guardan los cambios en cada iteración y se borra la matirz de cambio
    R=zeros(x+1,y+1,z+1);
    Space(:,:,:,k+1)=Space(:,:,:,k)+D;
    D1(:,:,:,k)=D(:,:,:);
    D(:,:,:)=0;
    %% Concentraciones
    N(k,:)=[length(find(Space(:,:,:,k)==1)),length(find(Space(:,:,:,k)==2)),length(find(Space(:,:,:,k)==3)),length(find(Space(:,:,:,k)==4)),length(find(Space(:,:,:,k)==5)),length(find(Space(:,:,:,k)==6))];
    suma_N=sum(N(k,:));
    xt(k,:)=N(k,:)/suma_N;
    pos(:,:,k+1)=pos2;
    pos2=zeros(n_Cu+n_CN,4);
end

xx=pos(:,1,:); yy=pos(:,2,:); zz=pos(:,3,:); TT=pos(:,4,:);
%% Simulación
for g=1:interaciones
    voxelSurf(Space(:,:,:,g),false);
    cb=colorbar;
    colormap ([0 1 1; 1 0 0; 0 0 1; 0 0 0; 1 1 0; 0 1 0]);
    grid on
    F=getframe;
end

%% Concentraciones
figure (1)
plot(xt);

%% Funciones
%%  Direccion del movimiento
function [m1,m2,m3,r]=m(i,j,l,P_ac,veci)
e=rand();                               % 0<e<0.25 (vecino de arriba) 0.25<e<0.5 (vecino derecho)
r=find(P_ac>e);   r=r(1)-1;             % 0.5<e<0.75 (vecino de abajo) 0.75<e<1 (vecino izquierdo)
mol=veci{i,j,l}(r,:);
m1=mol(1);  m2=mol(2); m3=mol(3);
end

%% Posiciones aleatorias al momento del decidir el movimiento
function [Posicion1]=Pos(Posicion,n_Cu,n_CN)
Posicion1=zeros(n_Cu+n_CN,4);
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
function [D,pos]=Movimiento(i,j,l,T,k,h,x,y,z,pos,D,P,pos1,Space,veci,n)
vecin=veci{i,j,l};  o=0;
r1i=vecin(1,1); r1j=vecin(1,2); r1l=vecin(1,3); r2i=vecin(2,1); r2j=vecin(2,2); r2l=vecin(2,3);
r3i=vecin(3,1); r3j=vecin(3,2); r3l=vecin(3,3); r4i=vecin(4,1); r4j=vecin(4,2); r4l=vecin(4,3);
r5i=vecin(5,1); r5j=vecin(5,2); r5l=vecin(5,3); r6i=vecin(6,1); r6j=vecin(6,2); r6l=vecin(6,3);

%% Se pregunta si el automata realizará el movimiento
e=rand();
if e<P
    pos(h,:)=[i j l T];
    o=1;
end
while o==0
    %% Se elige otra dirrección y se guardan los cambios
    P_ac=interacion(i,j,l,k,pos1,x,y,z,Space);
    [m1,m2,m3,~]=m(i,j,l,P_ac,veci);
    %% Se guarda el moviemiento en la matiz de cambio
    D(m1,m2,m3)=D(m1,m2,m3)+n;    D(i,j,l)=D(i,j,l)-n;           % Matiz de cambio "D"
    if D(m1,m2,m3)==n && Space(m1,m2,m3,k)==0
        pos(h,:)=[m1 m2 m3 T];
        break;
    end
    %% De estar ocupada la posicion, se borra el movimiento y se repite el proceso
    D(m1,m2,m3)=D(m1,m2,m3)-n;    D(i,j,l)=D(i,j,l)+n;           % Matiz de cambio "D"
    if Space(r1i,r1j,r1l,k)~=0&&Space(r2i,r2j,r2l,k)~=0&&Space(r3i,r3j,r3l,k)~=0&&Space(r4i,r4j,r4l,k)~=0&&Space(r5i,r5j,r5l,k)~=0&&Space(r6i,r6j,r6l,k)~=0
        pos(h,:)=[i j l T];
        break;
    end
end
end

%% Funcion Reaccion
function [D,R,pos]=Reaccion(i,j,l,T,k,h,R,D,Space,pos,veci,n)
%% Reacción
vec_ij=veci{i,j,l};
if Space(vec_ij(1,1),vec_ij(1,2),vec_ij(1,3),k)==2&&D(vec_ij(1,1),vec_ij(1,2),vec_ij(1,3))~=-2
    R(i,j,l)=1;   R(vec_ij(1,1),vec_ij(1,2),vec_ij(1,3))=1;
    D(i,j,l)=n;   D(vec_ij(1,1),vec_ij(1,2),vec_ij(1,3))=-2;
    pos(h,:)=[i;j;l;T+n];
elseif Space(vec_ij(2,1),vec_ij(2,2),vec_ij(2,3),k)==2&&D(vec_ij(2,1),vec_ij(2,2),vec_ij(2,3))~=-2
    R(i,j,l)=1;   R(vec_ij(2,1),vec_ij(2,2),vec_ij(2,3))=1;
    D(i,j,l)=n;   D(vec_ij(2,1),vec_ij(2,2),vec_ij(2,3))=-2;
    pos(h,:)=[i;j;l;T+n];
elseif Space(vec_ij(3,1),vec_ij(3,2),vec_ij(3,3),k)==2&&D(vec_ij(3,1),vec_ij(3,2),vec_ij(3,3))~=-2
    R(i,j,l)=1;   R(vec_ij(3,1),vec_ij(3,2),vec_ij(3,3))=1;
    D(i,j,l)=n;   D(vec_ij(3,1),vec_ij(3,2),vec_ij(3,3))=-2;
    pos(h,:)=[i;j;l;T+n];
elseif Space(vec_ij(4,1),vec_ij(4,2),vec_ij(4,3),k)==2&&D(vec_ij(4,1),vec_ij(4,2),vec_ij(4,3))~=-2
    R(i,j,l)=1;   R(vec_ij(4,1),vec_ij(4,2),vec_ij(4,3))=1;
    D(i,j,l)=n;   D(vec_ij(4,1),vec_ij(4,2),vec_ij(4,3))=-2;
    pos(h,:)=[i;j;l;T+n];
elseif Space(vec_ij(5,1),vec_ij(5,2),vec_ij(5,3),k)==2&&D(vec_ij(5,1),vec_ij(5,2),vec_ij(5,3))~=-2
    R(i,j,l)=1;   R(vec_ij(5,1),vec_ij(5,2),vec_ij(5,3))=1;
    D(i,j,l)=n;   D(vec_ij(5,1),vec_ij(5,2),vec_ij(5,3))=-2;
    pos(h,:)=[i;j;l;T+n];
elseif Space(vec_ij(6,1),vec_ij(6,2),vec_ij(6,3),k)==2&&D(vec_ij(6,1),vec_ij(6,2),vec_ij(6,3))~=-2
    R(i,j,l)=1;   R(vec_ij(6,1),vec_ij(6,2),vec_ij(6,3))=1;
    D(i,j,l)=n;   D(vec_ij(6,1),vec_ij(6,2),vec_ij(6,3))=-2;
    pos(h,:)=[i;j;l;T+n];
end
end

%% Calculo de las fuerzas Intermoleculares
function [P_ac]=interacion(i,j,l,k,pos1,x,y,z,Space)
m_ij=Space(i,j,l,k);
Q_ij=round(-0.075*m_ij^5+1.4167*m_ij^4-10.125*m_ij^3+33.583*m_ij^2-50.8*m_ij+27);
F1=0; F2=0; F3=0; F4=0; F5=0; F6=0;
if Q_ij~=0
    lo=[x,y,z]; K=8987551788;
    for d=1:length(pos1(:,1))
        rk=pos1(d,1:3);
        m_r=Space(rk(1),rk(2),rk(3),k);
        Q_r=round(-0.075*m_r^5+1.4167*m_r^4-10.125*m_r^3+33.583*m_r^2-50.8*m_r+27);
        if Q_ij*Q_r>0
            r=[i,j,l]-rk;
        else
            r=rk-[i,j,l];
        end
        if r(1)==0 && r(2)==0 && r(3)==0
        else
            %% Condicion de Imagen Minima
            for h=1:3
                if r(h)>lo(h)/2
                    r(h)=r(h)-lo(h);
                elseif r(h)<=-lo(h)/2
                    r(h)=r(h)+lo(h);
                end
            end
            rd=sqrt(r(2)^2+r(1)^2+r(3)^2);
            u=r/rd;            % [i,j,l] Vector unitario
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
            if u(3)>0
                F6=F6+abs(F*u(3));
            else
                F5=F5+abs(F*u(3));
            end
        end
    end
    Ft=F1+F2+F3+F4+F5+F6;
    P1=F1/Ft;  P2=F2/Ft; P3=F3/Ft; P4=F4/Ft; P5=F5/Ft;  P6=F6/Ft;
    P_ac=[0,P1,P1+P2,P1+P2+P3,P1+P2+P3+P4,P1+P2+P3+P4+P5,P1+P2+P3+P4+P5+P6];
else
    P_ac=[0 1/6 2/6 3/6 4/6 5/6 1];
end
end