clc;
clear all;
close all;

N = 10^6;
Eb_No = [0:25];
nRx = 2;
No_Err = zeros(1, length(Eb_No)); % Initialize error array
Sim_Ber = zeros(1, length(Eb_No)); % Initialize BER array

for ii = 1:length(Eb_No)
    RandS = rand(1, N) > 0.5;
    s = 2 * RandS - 1;
    
    STBCode = 1/sqrt(2) * kron(reshape(s, 2, N/2), ones(1, 2));
    Ray_ch = 1/sqrt(2) * (randn(nRx, N) + 1j * randn(nRx, N));
    Wh_gau = 1/sqrt(2) * (randn(nRx, N) + 1j * randn(nRx, N));
    
    y = zeros(nRx, N);
    S_Rx = zeros(nRx * 2, N);
    Nos_Eq = zeros(nRx * 2, N);
    
    for kk = 1:nRx
        No_Ray = kron(reshape(Ray_ch(kk, :), 2, N/2), ones(1, 2));
        temp = No_Ray;
        
        No_Ray(1, [2:2:end]) = conj(temp(2, [2:2:end]));
        No_Ray(2, [2:2:end]) = -conj(temp(1, [2:2:end]));
        
        y(kk, :) = sum(No_Ray .* STBCode, 1) + 10^(-Eb_No(ii)/20) * Wh_gau(kk, :);
        
        S_Rx([2 * kk - 1 : 2 * kk], :) = kron(reshape(y(kk, :), 2, N/2), ones(1, 2));
        Nos_Eq([2 * kk - 1 : 2 * kk], :) = No_Ray;
        
        Nos_Eq(2 * kk - 1, [1:2:end]) = conj(Nos_Eq(2 * kk - 1, [1:2:end]));
        Nos_Eq(2 * kk, [2:2:end]) = conj(Nos_Eq(2 * kk, [2:2:end]));
    end
    
    Nos_EqPower = sum(Nos_Eq .* conj(Nos_Eq), 1);
    EQ_s = sum(Nos_Eq .* S_Rx, 1) ./ Nos_EqPower;
    EQ_s(2:2:end) = conj(EQ_s(2:2:end));
    
    S_EQ_s = real(EQ_s) > 0;
    No_Err(ii) = size(find(RandS - S_EQ_s), 2);
end

Sim_Ber = No_Err / N; % Calculate simulated BER
EbNo_1 = 10.^(Eb_No / 10);
theoryBer_nRx1 = 0.5 .* (1 - 1 * (1 + 1 ./ EbNo_1) .^ (-0.5));
p = 1/2 - 1/2 * (1 + 1 ./ EbNo_1) .^ (-1/2);
theoryBerMRC_nRx2 = p .^ 2 .* (1 + 2 * (1 - p));
pAlamouti = 1/2 - 1/2 * (1 + 2 ./ EbNo_1) .^ (-1/2);
theoryBerAlamouti_nTx2_nRx1 = pAlamouti .^ 2 .* (1 + 2 * (1 - pAlamouti));

figure;
semilogy(Eb_No, theoryBer_nRx1, 'bp-', 'LineWidth', 2);
hold on;
semilogy(Eb_No, theoryBerMRC_nRx2, 'kd-', 'LineWidth', 2);
semilogy(Eb_No, theoryBerAlamouti_nTx2_nRx1, 'c+-', 'LineWidth', 2);
semilogy(Eb_No, Sim_Ber, 'mo-', 'LineWidth', 2); % Plot simulated BER
axis([0 25 10^-5 0.5]);
grid on;

legend('theory (nTx=1, nRx=1)', 'theory (nTx=1, nRx=2, MRC)', 'theory (nTx=2, nRx=1, Alamouti)', 'sim (nTx=2, nRx=2, Alamouti)');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation with 2Tx, 2Rx Alamouti STBC (Rayleigh Channel)');
