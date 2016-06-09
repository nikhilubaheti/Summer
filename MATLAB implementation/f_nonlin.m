function dx = f_nonlin(t,x,clf)
u = clf.controller(t,x);
if clf.commands.disp_u
    disp(u);
end
if clf.commands.disp_t
    disp(t+clf.t_last);
end
dx = clf.A*x+clf.B*u;
end