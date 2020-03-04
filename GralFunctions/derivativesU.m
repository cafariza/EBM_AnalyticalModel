function [Ux1D,Ux2D,Ux3D,Ux4D,Ux5D]=approxderivU(lw,L,E,H)
A=lw/10000;
Dx=A/L;
x=0:Dx:H/L;
Ux=zeros(length(x),length(E(1,:)));
for i=1:length(E(1,:))
UxUlti=E(4,i)*cos(b1(i)*x(end))+E(3,i)*sin(b1(i)*x(end))+E(1,i)*cosh(b2(i)*x(end))+E(5,i)*sinh(b2(i)*x(end))+E(2,i)*cosh(b3(i)*x(end))+E(6,i)*sinh(b3(i)*x(end));
UxCorr=(E(4,i)*cos(b1(i)*x(1))+E(3,i)*sin(b1(i)*x(1))+E(1,i)*cosh(b2(i)*x(1))+E(5,i)*sinh(b2(i)*x(1))+E(2,i)*cosh(b3(i)*x(1))+E(6,i)*sinh(b3(i)*x(1)))/UxUlti;
for g=1:length(x)
Ux(g,i)=(E(4,i)*cos(b1(i)*x(g))+E(3,i)*sin(b1(i)*x(g))+E(1,i)*cosh(b2(i)*x(g))+E(5,i)*sinh(b2(i)*x(g))+E(2,i)*cosh(b3(i)*x(g))+E(6,i)*sinh(b3(i)*x(g)))/UxUlti-UxCorr;
end
end
d1U=diff(Ux)/Dx; d2U=diff(d1U)/Dx; d3U=diff(d2U)/Dx; d4U=diff(d3U)/Dx;
d5U=diff(d4U)/Dx;

d1U(length(x),:)=d1U(end,:);
d2U(length(x)-1,:)=d2U(end,:);
d2U(length(x),:)=d2U(end,:);
d3U(length(x)-2,:)=d3U(end,:);
d3U(length(x)-1,:)=d3U(end,:);
d3U(length(x),:)=d3U(end,:);
d4U(length(x)-3,:)=d4U(end,:);
d4U(length(x)-2,:)=d4U(end,:);
d4U(length(x)-1,:)=d4U(end,:);
d4U(length(x),:)=d4U(end,:);
d5U(length(x)-4,:)=d5U(end,:);
d5U(length(x)-3,:)=d5U(end,:);
d5U(length(x)-2,:)=d5U(end,:);
d5U(length(x)-1,:)=d5U(end,:);
d5U(length(x),:)=d5U(end,:);

xB=0:lw/L:H/L;
for i=1:length(x)
    for j=1: length(xB)
if (x(i)-xB(j))<Dx
Ux1D(j,:)=d1U(i,:);
Ux2D(j,:)=d2U(i,:);
Ux3D(j,:)=d3U(i,:);
Ux4D(j,:)=d4U(i,:);
Ux5D(j,:)=d5U(i,:);
end

end
end

end