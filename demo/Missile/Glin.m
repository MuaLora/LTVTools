function xdot = Glin(t,x,deltaPert,tPert,A,B,Tgrid)

% Interpolate delta
if t<=tPert(1)
    u = deltaPert(1);
elseif t>=tPert(end)
    u = deltaPert(end);
else
    idx = find(t>=tPert(1:end-1) & t<tPert(2:end));
    fac = (t-tPert(idx))/(tPert(idx+1) - tPert(idx));
    u = (1-fac)*deltaPert(idx)+fac*deltaPert(idx+1);
end

% Interpolate state matrices
if t<=Tgrid(1)
    At = A(:,:,1);
    Bt = B(:,:,1);
elseif t>=Tgrid(end)
    At = A(:,:,end);
    Bt = B(:,:,end);
else
    idx = find(t>=Tgrid(1:end-1) & t<Tgrid(2:end));
    fac = (t-Tgrid(idx))/(Tgrid(idx+1) - Tgrid(idx));
    At = (1-fac)*A(:,:,idx)+fac*A(:,:,idx+1);
    Bt = (1-fac)*B(:,:,idx)+fac*B(:,:,idx+1);
end

xdot = At*x+Bt*u;