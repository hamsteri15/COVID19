function [dU,dV,dW] = project(dU,dV,dW,KX,KY,KZ,AA,OnePerK,KXXP,KYYP,KZZP,KXYP,KXZP,KYZP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projects the input field dU,dV (or U,V)
% and returns dU,dV pair which is divergence-free.
% The idea is to ensure that the continuity <k, u_k> = 0 
% equation in Fourier transformed form.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fft of a field with divergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FU = AA.*fftn(dU);
FV = AA.*fftn(dV);
FW = AA.*fftn(dW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection so that the wave vector  (kx,ky) is orthogonal to (fu,fv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUU = FU-(KX.*FU + KY.*FV  + KZ.*FW).*KX.*OnePerK;
%FVV = FV-(KX.*FU + KY.*FV  + KZ.*FW).*KY.*OnePerK;
%FWW = FW-(KX.*FU + KY.*FV + KZ.*FW).*KZ.*OnePerK; 

FUU = FU-(KXXP.*FU + KXYP.*FV  + KXZP.*FW);
FVV = FV-(KXYP.*FU + KYYP.*FV  + KYZP.*FW);
FWW = FW-(KXZP.*FU + KYZP.*FV + KZZP.*FW);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse Fourier transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dU = real(ifftn(FUU));
dV = real(ifftn(FVV));
dW = real(ifftn(FWW));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now (dU,dV) is a div free pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
