function mAc=Mac(Phi1,Phi2)
% This function calculates mac between phi1 and phi2
% The MAC value between two modes is essentially the normalized dot product of the complex modal vector at each common nodes (i.e., points), as shown in Equation 1.  It can also be thought of as the square of correlation between two modal vectors ?r and ?s. 
%  Equation 1: Modal Assurance Criterion equation for comparing two mode shapes
% If a linear relationship exists (i.e., the vectors move the same way) between the two complex vectors, the MAC value will be near to one. If they are linearly independent, the MAC value will be small (near zero).
 
mAc= (abs(Phi1'*Phi2))^2/((Phi1'*Phi1)*(Phi2'*Phi2));
end