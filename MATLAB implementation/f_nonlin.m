function dx = f_nonlin(t,x,clf)
u = clf.controller(t,x);
% disp(u);
disp(t);
dx = clf.A*x+clf.B*u;
end