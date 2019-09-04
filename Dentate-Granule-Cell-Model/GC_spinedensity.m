function [strct] = GC_spinedensity(fac)

if nargin < 1 || isempty(fac)
    fac = 1;
end


strct = struct;


strct.axonh.spines =     struct('scale',1+0*fac);%struct('density',5);
strct.axon.spines =     struct('scale',1+0*fac);%struct('density',5);
strct.soma.spines =     struct('scale',1+0*fac);%struct('density',5);
strct.GCL.spines =      struct('scale',1+0*fac);%struct('density',5);
strct.SGCL.spines =     struct('scale',1+0*fac);%struct('density',5);
strct.adendIML.spines = struct('scale',1+0.5*fac);%struct('density',5);
strct.adendMML.spines = struct('scale',1+1*fac);%2);%struct('density',5);
strct.adendOML.spines = struct('scale',1+1*fac);%2.5);%struct('density',5);



