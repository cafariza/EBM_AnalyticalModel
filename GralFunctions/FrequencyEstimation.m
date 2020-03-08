%% HPDM: NATURAL FREQUENCY COMPUTATIONS

omegar=zeros(d,1);
% --------------------------------------------------------------------------
% Global structural behavior
% --------------------------------------------------------------------------
x=zeros(d,1);
y=zeros(d,1);
C=zeros(d,1);
L=zeros(d,1);
epsilon=zeros(d,1);

for k=1:d
% Scale Factor
L(k)=2*N*lw/(pi*(2*k-1));
epsilon(k)=pi*(2*k-1)/(2*N);
% Macroscopic Constants
C(k)=Kgb/(K*(L(k))^2);
gama=Ki/Kgb;
xE(k)=log(C(k))/log(epsilon(k));
yE(k)=log(gama)/log(epsilon(k));
end

% --------------------------------------------------------------------------
% Frequency estimation
% --------------------------------------------------------------------------

omega=(2*pi()*(0.8/N+0.006)):(2*pi()*(0.8/N+0.006)):30*100*pi()/N;
deter=zeros(length(omega),1); 
frange=cell(d,1,1);
% Scale Factor
L=2*N(1)*lw/(pi);
% Macroscopic Constants
C=Kgb/(K*(L)^2);
beta=C*(L/lw)^2;
gama=Ki/Kgb;

for i=1:length(omega)
omegamsq1=Lam*omega(i)^2*(L)^2/K;

p=-((1+gama)^2/(3*C^2*gama^2)+omegamsq1/(C*gama));
q=2*(1+gama)^3/(27*C^3*gama^3)-omegamsq1/(C^2*gama)*(1-(1+gama)/(3*gama));
eqn1=1/3*acos((q/2)*sqrt(-27/p^3));

B1=2*sqrt(-p/3)*cos(eqn1+2*pi/3)+(1+gama)/(3*C*gama);
B2=2*sqrt(-p/3)*cos(eqn1+4*pi/3)+(1+gama)/(3*C*gama);
B3=2*sqrt(-p/3)*cos(eqn1)+(1+gama)/(3*C*gama);

b1=sqrt(-B1);
b2=sqrt(B2);
b3=sqrt(B3);
deter(i)=(cos(pi*b1/2)*(B1^2-B2^2)*(B1^2-B3^2)*(2*B2*B3/(cosh(pi*b2/2)*cosh(pi*b3/2))+tanh(pi*b2/2)*tanh(pi*b3/2)*b2*b3*(B2+B3))+(B2^2-B1^2)*(B2^2-B3^2)*(2*B1*B3/(cosh(pi*b3/2))-sin(pi*b1/2)*tanh(pi*b3/2)*b1*b3*(B1+B3))+(B3^2-B1^2)*(B3^2-B2^2)*(2*B1*B2/(cosh(pi*b2/2))-sin(pi*b1/2)*tanh(pi*b2/2)*b1*b2*(B1+B2))-cos(pi*b1/2)*(B1^2*(B2^2-B3^2)^2+B3^2*(B2^2-B1^2)^2+B2^2*(B1^2-B3^2)^2));
end
idx = imag(double(deter)) == 0;
deter = deter(idx);
deter=real(deter); %Taking the real number
%Looking for proper natural frequency range 
n=0;

for i=1:length(deter)-1
   if deter(i)*deter(i+1)<0
        frange{i}=[omega(i)^2 omega(i+1)^2];
        n=n+1;
        
        switch n
            case 10
                break
        end
    end
end
frange= frange(~cellfun('isempty',frange));

syms omegasq 
omegamsq=sqrt(Lam*omegasq*(L)^2/K);

p=-((1+gama)^2/(3*C^2*gama^2)+omegamsq^2/(C*gama));
q=2*(1+gama)^3/(27*C^3*gama^3)-omegamsq^2/(C^2*gama)*(1-(1+gama)/(3*gama));
eqn1=1/3*acos((q/2)*sqrt(-27/p^3));

B1=2*sqrt(-p/3)*cos(eqn1+2*pi/3)+(1+gama)/(3*C*gama);
B2=2*sqrt(-p/3)*cos(eqn1+4*pi/3)+(1+gama)/(3*C*gama);
B3=2*sqrt(-p/3)*cos(eqn1)+(1+gama)/(3*C*gama);

b1=sqrt(-B1);
b2=sqrt(B2);
b3=sqrt(B3);

eqn=(cos(pi*b1/2)*(B1^2-B2^2)*(B1^2-B3^2)*(2*B2*B3/(cosh(pi*b2/2)*cosh(pi*b3/2))+tanh(pi*b2/2)*tanh(pi*b3/2)*b2*b3*(B2+B3))+(B2^2-B1^2)*(B2^2-B3^2)*(2*B1*B3/(cosh(pi*b3/2))-sin(pi*b1/2)*tanh(pi*b3/2)*b1*b3*(B1+B3))+(B3^2-B1^2)*(B3^2-B2^2)*(2*B1*B2/(cosh(pi*b2/2))-sin(pi*b1/2)*tanh(pi*b2/2)*b1*b2*(B1+B2))-cos(pi*b1/2)*(B1^2*(B2^2-B3^2)^2+B3^2*(B2^2-B1^2)^2+B2^2*(B1^2-B3^2)^2))==0;

for n=1:length(frange)
omegastor=vpasolve(eqn,omegasq,[frange{n}(1),frange{n}(2)],'random',true);
if isempty(omegastor)
omegastor=0;
end
omegar(n)=omegastor;
end

wn=sort(omegar);
Ndecimals = 5 ;
fround = 10.^Ndecimals; 
wn = round(fround*wn)/fround;
inisize_wn=length(wn);
omesqu=wn;
for i=1:length(wn)
if wn(i)==0
omesqu(i) = [];
end
end
wn=omesqu(1:4);

    for s=1:d
    f_natural(s)=sqrt(wn(s))/(2*pi);
    f_ratio(s)=f_natural(s)/f_natural(1);
    end
    
f_natural    %Hz
f_ratio