%% --------------------------------------------------------------------------
% Internal Forces (ANALYTICAL RESULTS)
% --------------------------------------------------------------------------
Omega_graph=omega; %Saving omega for graph in plot function
omega=sqrt(wn); % natural frequencies
xd=x;
 for k=1:length(E(1,:))
     for g=1:length(x)
     for i=1:nw
     Uy{k}(g,i)=alpha(g,k)*UY_numuA(i);%REVISION IS REQUIRED HERE
     end
     end
 end
lfac=zeros(nw,1);
 for i=1:nw
     if i==1
     lfac(i)=0;
     else
         lfac(i)=lf(i-1)+lfac(i-1);
     end
        
 end
[thetaa, thetan, nodal_rot]=rotation(af,aw,hf,hw,lf,lw,Ew,Uy,Ux1, xd,alpha,pe,omega,Theta_numuA,Theta_numuU,L ,length(aw),length(af));
[Na, T, M]=Internalforces(af,aw,hf,hw,lf,lw,Ew,Ux,Uy,Ux1,Ux2,Ux3, xd,alpha,alpha1,alpha2,pe,omega,nodal_rot,L,length(aw),length(af));
ForcesHPDM
  for h=1:length(omega)
     for i=1:nw
     for g=1:length(x)
         M_glb{h}(:,i)=(M_global{h}(:,i)*(lfac(i)-lfac(end-(i-1)))/(sum(lf)));
     end
     end
 end
 
 Nam={};
 for h=1:length(omega)
     for i=1:nw
     for g=1:length(x)
         Nam{h}(:,i)=M_glb{h}(:,i)/(sum(lf));
     end
     end
 end
 