%%% simple plot of rc(w)

w = linspace(0,0.25);
E = 16.5;
A = 0.0198;
n = 2.0455;

rc = zeros(length(w),1);

for i = 1:length(w)
    rc(i) = (E*w(i)^n)/(A^n + w(i)^n);
end

hold on; box on
plot(w,rc, 'Linewidth',2)
xlabel('Oxygen')
ylabel('Climax Growth Rate')
title('Climax Growth Rate')