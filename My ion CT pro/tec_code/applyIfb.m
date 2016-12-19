function data=applyIfb(data,ifb);

% corrects the data structure for the IFB
A = 40.3;           % ionosphere constant
f1 = 1.57542e9;     % L1 frequency
f2 = 1.2276e9;      % L2 frequency
c = 0.299792458e9;  % speed of light

F = (f1^2*f2*c) / (A * (f1^2-f2^2));
names=fieldnames(data);
for n=1:length(names)
  name=cell2mat(names(n));
  eval(['data.' name '(:,3)=data.' name '(:,3)+F*f2*ifb;']);
end
