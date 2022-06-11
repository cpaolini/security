function [sigout]=plotfunc(SNR_vec,successRate,command)
semilogy(SNR_vec,flip(successRate))
grid on
% axis([-Inf Inf 10^(-4) 1])
xlabel('normalized SNR(dB) ')
ylabel('success rate')
legend('Modified method success rate')
sigout=2;
end





