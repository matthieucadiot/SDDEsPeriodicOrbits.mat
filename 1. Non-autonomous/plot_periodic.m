function [y,mesh] = plot_periodic(a,t)
% Plot of the trigonometric polynomial with Fourier coefficients a = (a_n)_{|n|<=N} on [-pi,pi].

color = [0 0 1];
mesh = (-pi:0.0001:pi);

y = real(sol(a,mesh));

plot(mesh,y,'Color',color,'LineWidth',4)

set(gca,'FontSize',20)

axis tight
xlabel('$$t$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$\bar u(t)$$', 'Interpreter', 'latex', 'FontSize', 30)

end

function s = sol(a,x)
% sum_n a_n e^{inx}
m = (length(a)+1)/2;
s = 0;
for k=-m+1:m-1
    s = s + a(k+m,1)*exp(1i*k*x);
end
end
