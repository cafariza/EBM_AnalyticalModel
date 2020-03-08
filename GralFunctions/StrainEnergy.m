%Strain Energy Calculation
%SE=1/2*int_xi_xf {T(x)^2/K+Mi(x)^2/Ki+Mg(x)^2/Kg}dx 

dx=xd(2)-xd(1);
SE_T=zeros(length(xd)-1,length(omega));
SE_Mi=zeros(length(xd)-1,length(omega));
SE_Mg=zeros(length(xd)-1,length(omega));
SE=zeros(length(xd)-1,length(omega));
for j=1: length(omega)
for i=1: length(xd)-1
    SE_T(i,j)=1/2*(dx*((Tr{j}(i,1))^2+(Tr{j}(i+1,1))^2)/2)/K;           %Strain Energy due to SHEAR
    SE_Mi(i,j)=1/2*(dx*((M_in{j}(i,1))^2+(M_in{j}(i+1,1))^2)/2)/Ki;     %Strain Energy due to Inner Bending
    SE_Mg(i,j)=1/2*(dx*((M_glb{j}(i,1))^2+(M_glb{j}(i+1,1))^2)/2)/Kgb;  %Strain Energy due to GlobalBending
    SE(i,j)= SE_T(i,j)+SE_Mi(i,j)+SE_Mg(i,j);                           % Total Strain Energy 
end
Total_SE(j)=sum(SE(:,j));
end

%Plot Parameters 
ImageFontSizeM=11;
ImageFontSize=9;
LegendFontSize=7;
MarkSize=6;
FontName='Garamond';
AxisFontName='Garamond';
axes1 = axes('FontSize',ImageFontSize,'FontName',AxisFontName);
box(axes1,'on');
hold(axes1,'all');
figure 
subplot(1,3,1)
 fig1=barh(xd(1:(end)), [SE(:,1); SE(end,1)]);
 set (fig1(1), 'FaceColor', [.7 .7 .7]);
 
 ylabel('x [mm]','FontSize',ImageFontSize,'FontName',AxisFontName);

  subplot(1,3,2)
 fig1=barh(xd(1:(end)), [SE(:,2); SE(end,2)]);
 set (fig1(1), 'FaceColor', [.7 .7 .7]);
 title(' Strain energy distribution along the structure','FontSize',ImageFontSizeM,'FontName',FontName);
xlabel('Strain energy [MN mm]','FontSize',ImageFontSize,'FontName',AxisFontName); 
 subplot(1,3,3)
 fig1=barh(xd(1:(end)), [SE(:,3); SE(end,3)]);
 set (fig1(1), 'FaceColor', [.7 .7 .7]);
