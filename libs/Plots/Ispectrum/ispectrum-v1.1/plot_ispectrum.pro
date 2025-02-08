
o=obj_new('data','ispectrum.dat')

openpr,'plot_ispectrum.ps',/std
mu = o->column(1)
int = o->column(2)
ind = where(int ne 0.0)
plot,mu,int,/xlog,/ylog,psym=-sym(1),symsize=0.5,/xstyle
;plot,mu(ind),int(ind),/xlog,/ylog,psym=-sym(1),symsize=0.5,/xstyle
closepr

end
