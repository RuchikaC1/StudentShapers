dp = 10;
mu = 1;
r = 0.08:1e-3:0.2;
por = 1-pi*(r.^2/0.16);
K = 4*r.^2.*por.^3./(16*5.75*(1-por).^2);
expected = 9*mu*por./(2*K*dp);
actual = [64 272 130 695];
por2 = [0.803650459 0.558213533 0.693203842 0.398679531];


plot(1-por, expected, "k--", LineWidth=1)
hold on
scatter(1-por2, actual, "kx", Linewidth=5)
grid on
xlim([0.15 0.65]);
xlabel("Fibre volume fraction")
ylabel("Infiltration time (s)")


