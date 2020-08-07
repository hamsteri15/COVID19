U=real(ifftn( LP.*fftn(U)));
V=real(ifftn( LP.*fftn(V)));
W=real(ifftn( LP.*fftn(W)));
T=real(ifftn( LP.*fftn(T)));

T = min(1,max(0,T));