% --------------------------------------------------------------------------
%% Post-Processing of NUMERICAL RESULTS
% --------------------------------------------------------------------------

%Vertical displacement and rotation
fileName=sprintf('DISPV1.inp');[ Vx1] = ExtinftxtMF( fileName, N, nw); fileName=sprintf('DISPV2.inp');
[ Vx2] = ExtinftxtMF( fileName, N, nw); fileName=sprintf('DISPV3.inp');
[ Vx3] = ExtinftxtMF( fileName, N, nw); 
LastData1=(nw-1)*N+nw;
LastData2=nw*N+nw;

alphan1=(Vx1(1:(N+1),3)-Vx1(LastData1:LastData2,3))/sum(lf) ;
alphan2=(Vx2(1:(N+1),3)-Vx2(LastData1:LastData2,3))/sum(lf); 
alphan3=(Vx3(1:(N+1),3)-Vx3(LastData1:LastData2,3))/sum(lf); 
alphan4=(Vx3(1:(N+1),3)-Vx3(LastData1:LastData2,3))/sum(lf) ;
alphan=[alphan1 alphan2 alphan3 alphan4];
alphanum={};
for j=1:4
    for i=1:nw
alphanum{j}(:,i)=alphan(:,j);
    end
end


%Transverse displacement
Unum={};
fileName=sprintf('DISPT1.inp');[Unum1] = ExtinftxtMF( fileName, N, nw); fileName=sprintf('DISPT2.inp');[Unum2] = ExtinftxtMF( fileName, N, nw);
fileName=sprintf('DISPT3.inp');[Unum3] = ExtinftxtMF( fileName, N, nw); 
    for i=1:nw
        Unum{1}(:,i)=Unum1((i-1)*N+i:i*N+i,3);
        Unum{2}(:,i)=Unum2((i-1)*N+i:i*N+i,3);
        Unum{3}(:,i)=Unum3((i-1)*N+i:i*N+i,3);
        Unum{4}(:,i)=Unum3((i-1)*N+i:i*N+i,3);
    end
    
    
figure     
plot(Unum{1,1}(:,1),x,Ux(:,1),x ) ; grid on; hold on;
 
xlabel('Displacement U*','FontSize',12); ylabel('x*','FontSize',12);
legend(': Mode 1: castem ', ': Mode 1: analytical');
title(' Modal Shapes ','FontSize',12);
hold off;

Uxunum1={};
Uxunum2={};
Uxunum3={};
alphanum1={};
%Derivatives of transverse displacement
h=x(2)-x(1);
for k=1: length(omega)
    for i=1:nw
Uxunum1{k}(:,i)=diff(Unum{k}(:,i))/h;
Uxunum2{k}(:,i)=diff(Uxunum1{k}(:,i))/h;
Uxunum3{k}(:,i)=diff(Uxunum2{k}(:,i))/h;
alphanum1{k}(:,i)=diff(alphanum{k}(:,i))/h;
    end
end

%Nodal rotation
fileName=sprintf('ROT1.inp');
[ Thetanum1] = ExtinftxtMF( fileName, N, nw); 
fileName=sprintf('ROT2.inp');
[ Thetanum2] = ExtinftxtMF( fileName, N, nw); 
fileName=sprintf('ROT3.inp');
[ Thetanum3] = ExtinftxtMF( fileName, N, nw); 
thetanum={};

    for i=1:nw
        thetanum{1}(:,i)=Thetanum1((i-1)*N+i:i*N+i,3);
        thetanum{2}(:,i)=Thetanum2((i-1)*N+i:i*N+i,3);
        thetanum{3}(:,i)=Thetanum3((i-1)*N+i:i*N+i,3);
        thetanum{4}(:,i)=Thetanum3((i-1)*N+i:i*N+i,3);
    end
%Internal forces BASED ON ANALYTICAL FORMULATION
[Tnum, Mnum, M_innum, M_in1num, M_globalnum, Trnum]=forces_M2(af,aw,hf,hw,lf,lw,Ew,Unum,Uxunum1,Uxunum2,Uxunum3,xd,pe,omega,L,thetanum,alphanum,alphanum1,length(aw),length(af));
%Internal forces From CASTEM
   
fileName=sprintf('IFN1MOD1.inp'); [FCastemW1_1]=ExtinftxtForces( fileName, N);fileName=sprintf('IFN1MOD2.inp'); [FCastemW1_2]=ExtinftxtForces( fileName, N);
fileName=sprintf('IFN1MOD3.inp'); [FCastemW1_3]=ExtinftxtForces( fileName, N);fileName=sprintf('IFN2MOD1.inp'); [FCastemW2_1]=ExtinftxtForces( fileName, N);
fileName=sprintf('IFN2MOD2.inp'); [FCastemW2_2]=ExtinftxtForces( fileName, N);fileName=sprintf('IFN2MOD3.inp'); [FCastemW2_3]=ExtinftxtForces( fileName, N);
j=0;
TCastem={};
for i=1:nw
    if i==1 || i==nw
    TCastem{1}(:,i)=FCastemW1_1{1}(:,4)*-1;
    TCastem{2}(:,i)=FCastemW1_2{1}(:,4)*-1;
    TCastem{3}(:,i)=FCastemW1_3{1}(:,4)*-1;
    else
    j=j+1;    
    TCastem{1}(:,i)=FCastemW2_1{j}(:,4)*-1;
    TCastem{2}(:,i)=FCastemW2_2{j}(:,4)*-1;
    TCastem{3}(:,i)=FCastemW2_3{j}(:,4)*-1;
    
    end
end

j=0;
AxCastem={};
for i=1:nw
    if i==1 || i==nw
    AxCastem{1}(:,2)=FCastemW1_1{1}(:,3)*-1;
    AxCastem{2}(:,2)=FCastemW1_2{1}(:,3)*-1;
    AxCastem{3}(:,2)=FCastemW1_3{1}(:,3)*-1;
    else
        j=j+1;
    AxCastem{1}(:,i)=FCastemW2_1{j}(:,3)*-1;
    AxCastem{2}(:,i)=FCastemW2_2{j}(:,3)*-1;
    AxCastem{3}(:,i)=FCastemW2_3{j}(:,3)*-1;
    end
end