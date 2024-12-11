clear;
N = 10^6; % number of bits or symbols
Eb_N0_dB = 0:25; % range of Eb/N0 values in dB
nTx = 2; % number of transmit antennas
nRx = 2; % number of receive antennas

for ii = 1:length(Eb_N0_dB)
    % Transmitter
    ip = rand(1, N) > 0.5; % generating bits 0,1 with equal probability
    s = 2 * ip - 1; % BPSK modulation: 0 -> -1; 1 -> 1
    sMod = kron(s, ones(nRx, 1)); % expand symbol for each Rx antenna
    sMod = reshape(sMod, [nRx, nTx, N/nTx]); % organize as [nRx, nTx, N/nTx] matrix
    
    % Channel and noise
    h = 1/sqrt(2) * (randn(nRx, nTx, N/nTx) + 1i * randn(nRx, nTx, N/nTx)); % Rayleigh channel
    n = 1/sqrt(2) * (randn(nRx, N/nTx) + 1i * randn(nRx, N/nTx)); % AWGN noise
    
    % Transmit signal through the channel with AWGN noise addition
    y = squeeze(sum(h .* sMod, 2)) + 10^(-Eb_N0_dB(ii)/20) * n;
    
    % Receiver with Zero Forcing (ZF) Equalizer
    hCof = zeros(2, 2, N/nTx); 
    hCof(1, 1, :) = sum(h(:, 2, :) .* conj(h(:, 2, :)), 1); % d term
    hCof(2, 2, :) = sum(h(:, 1, :) .* conj(h(:, 1, :)), 1); % a term
    hCof(2, 1, :) = -sum(h(:, 2, :) .* conj(h(:, 1, :)), 1); % c term
    hCof(1, 2, :) = -sum(h(:, 1, :) .* conj(h(:, 2, :)), 1); % b term
    
    hDen = hCof(1, 1, :) .* hCof(2, 2, :) - hCof(1, 2, :) .* hCof(2, 1, :); % ad-bc term
    hDen = reshape(kron(reshape(hDen, 1, N/nTx), ones(2, 2)), 2, 2, N/nTx);
    hInv = hCof ./ hDen; % ZF Equalizer matrix inv(H^H*H)
    
    % Process received symbols
    hMod = reshape(conj(h), nRx, N); % H^H
    yMod = kron(y, ones(1, 2)); % formatting y for equalization
    yMod = sum(hMod .* yMod, 1); % H^H * y
    yMod = kron(reshape(yMod, 2, N/nTx), ones(1, 2));
    yHat = sum(reshape(hInv, 2, N) .* yMod, 1); % inv(H^H*H) * H^H * y
    
    % Hard decision decoding
    ipHat = real(yHat) > 0;
    
    % Error counting
    nErr(ii) = size(find([ip - ipHat]), 2);
end

% Calculate Bit Error Rate (BER)
simBer = nErr / N;

% Theoretical BER calculations
EbN0Lin = 10.^(Eb_N0_dB / 10);
theoryBer_nRx1 = 0.5 * (1 - (1 + 1 ./ EbN0Lin).^(-0.5)); 
p = 1 / 2 - 1 / 2 * (1 + 1 ./ EbN0Lin).^(-0.5);
theoryBerMRC_nRx2 = p.^2 .* (1 + 2 * (1 - p)); 

% Plotting
figure;
semilogy(Eb_N0_dB, theoryBer_nRx1, 'bp-', 'LineWidth', 2);
hold on;
semilogy(Eb_N0_dB, theoryBerMRC_nRx2, 'kd-', 'LineWidth', 2);
semilogy(Eb_N0_dB, simBer, 'mo-', 'LineWidth', 2);
axis([0 25 10^-5 0.5]);
grid on;
legend('Theory (nTx=1, nRx=1)', 'Theory (nTx=1, nRx=2, MRC)', 'Simulated (nTx=2, nRx=2, ZF)');
xlabel('Average Eb/N0 (dB)');
ylabel('Bit Error Rate');
title('BER for BPSK modulation with 2x2 MIMO and ZF equalizer (Rayleigh channel)');
