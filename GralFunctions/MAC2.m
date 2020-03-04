function mac=MAC2(phi1,phi2)
% This function is to compute and plot Modal Assurance Criterion (MAC) matrix between identified mode shapes
% rectangle around the peaks.
% phi: matrix of the identified mode shapes
% mac: MAC matrix
% Example: load ModeShapes.mat;MAC_Matrix=MAC(phi)
for I=1:size(phi1,2)
    for J=1:size(phi2,2)
        mac(I,J)=Mac(phi1(:,I),phi2(:,J));
    end
end
% plot mac matrix
figure
h=bar3(mac); grid on; hold on;
set(h(1),'facecolor','red');
set(h(2),'facecolor','blue');
set(h(3),'facecolor','green');
title('MAC criterion:','FontSize',15)
xlabel('CASTEM mode shapes','FontSize',11); ylabel('FEM Mode shapes','FontSize',11);
legend(': Data mode 1', ': Data mode 2',': Data mode 3' );
end