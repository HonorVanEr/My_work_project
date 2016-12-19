function check_splint = check_sample(epochs,splint)

% extract vector with seconds (rounded), shift, and subtract
t1 = [round(epochs(:,6)) ; 0];
t2 = [0 ; round(epochs(:,6))];
t = abs(t1-t2);
t = t(2:length(t-1));

% if t is not within 1 sec of splint or equal to zero, then return 1
[I] = find(t>splint+1 & t<splint-1);
if (find(t(I) ~= 0))
  J = find(t(I) ~= 0);
  check_splint = size(J);
else
  check_splint = 0;
end
