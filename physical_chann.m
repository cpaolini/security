%% Channel
% Send signal over an AWGN channel.
function [data_rx]=channel_d(data_tx,snr)

% snr=10^(-snr/10)
data_rx = awgn(data_tx,snr,'measured');
