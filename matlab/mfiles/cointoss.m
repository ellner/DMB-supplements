
nt = 100000;
%bins = 20;
p = 0.6;                                      % probability of heads
r = rand(1,nt+1); 
rd = (r < p); % tails = 0, heads = 1
rdiff = rd(1:nt)+2*rd(2:nt+1);                %transitions: tt = 0, ht = 1, th = 2, hh = 3
tht = find(rdiff == 1);                       % heads to tails indices
tth = find(rdiff == 2);                       % tails to heads indices
rn  = min(length(tht),length(tth));           % numbers of transitions
if rd(1) == 0                                 % tails first
   rh = tht(1:rn) - tth(1:rn);                % heads residence times
   rt = [tth(1),tth(2:rn) - tht(1:rn-1)];     % tails residence times
else                                          % heads first
   rh = [tht(1),tht(2:rn) - tth(1:rn-1)];     % heads residence times
   rt = tth(1:rn) - tht(1:rn);                % tails residence times
end;
hhist = histc(rh,[1:max(rh)]);                % heads histogram
thist = histc(rt,[1:max(rt)]);                % tails histogram

hf = figure; bar(hhist);                      %plot heads histogram
tf = figure; bar(thist);                      %plot tails histogram

logh = figure; semilogy([1:length(hhist)],hhist,'b',[1:length(thist)],thist,'r')
