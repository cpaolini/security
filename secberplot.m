function [sigout]=secberplot(SNR_vec,successRate,command)
semilogy(SNR_vec,successRate)
grid on
% axis([-Inf Inf 10^(-4) 1])
xlabel('normalized SNR(dB) ')
ylabel('security rate')
legend('modified method security rate')
sigout=2;
end
