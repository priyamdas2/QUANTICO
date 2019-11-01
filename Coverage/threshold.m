function a = threshold(a,t)

%if strcmp(method,'soft')
%    x = sign(a).*(abs(a)-t).*(abs(a)>t);
%else
    %x = a.*(abs(a)>t);
    a(abs(a)<=t)=0;
%end
  
