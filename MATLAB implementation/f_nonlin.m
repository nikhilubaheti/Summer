function dx = f_nonlin(t,x,clf)
u = clf.controller(t,x);
clf.u = [clf.u;u];
dx = clf.A*x+clf.B*u;
end