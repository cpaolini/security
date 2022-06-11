
clc
clear all
close all

%% MIMO-OFDM Precoding with Phased Arrays


s = rng(61);        % Set RNG state for repeatability

%% System Parameters

% Single-user system with multiple streams
prm.numUsers = 1;          
prm.numSTS = 16;           
prm.numTx = 32;               
prm.numRx = 16;              
prm.bitsPerSubCarrier = 4;   
prm.numDataSymbols = 10;     

prm.fc = 1e9;                   
prm.chanSRate = 100e6;          
prm.ChanType = 'Scattering';    
                                          
prm.NFig = 50;                  

% Array locations and angles
prm.posTx = [0;0;0];            
prm.mobileRange = 300;           

prm.mobileAngle = [33; 0];      
prm.steeringAngle = [30; -20];  
prm.enSteering = true;          

%% 
% Parameters to define the OFDM modulation employed for the system are
% specified below.

prm.FFTLength = 256; 
prm.CyclicPrefixLength = 64; 
prm.numCarriers = 234; 
prm.NumGuardBandCarriers = [7 6];
prm.PilotCarrierIndices = [26 54 90 118 140 168 204 232];
nonDataIdx = [(1:prm.NumGuardBandCarriers(1))'; prm.FFTLength/2+1; ...
              (prm.FFTLength-prm.NumGuardBandCarriers(2)+1:prm.FFTLength)'; ...
              prm.PilotCarrierIndices.';];
prm.CarriersLocations = setdiff((1:prm.FFTLength)',sort(nonDataIdx));

numTx = prm.numTx;
numRx = prm.numRx;
numSTS = prm.numSTS;
prm.numFrmBits = numSTS*prm.numDataSymbols*prm.numCarriers* ...
                 prm.bitsPerSubCarrier*1/3-6; % Account for termination bits

prm.modMode = 2^prm.bitsPerSubCarrier; % Modulation order
% Account for channel filter delay
prm.numPadZeros = 3*(prm.FFTLength+prm.CyclicPrefixLength); 

% Get transmit and receive array information
prm.numSTSVec = numSTS;
[isTxURA,expFactorTx,isRxURA,expFactorRx] = helperArrayInfo(prm,true);

%%

prm.cLight = physconst('LightSpeed');
prm.lambda = prm.cLight/prm.fc;
% Mobile position
[xRx,yRx,zRx] = sph2cart(deg2rad(prm.mobileAngle(1)),...
                         deg2rad(prm.mobileAngle(2)),prm.mobileRange);
prm.posRx = [xRx;yRx;zRx];
[toRxRange,toRxAng] = rangeangle(prm.posTx,prm.posRx);
spLoss = fspl(toRxRange,prm.lambda);
gainFactor = 1;

% noise_pow=1:2:10;
%   for ss=1:1:length(noise_pow)
%         snr_db=noise_pow(ss);
        
        noise_pow=10:10:100;
  for ss=1:1:length(noise_pow)
        prm.NFig=noise_pow(ss);

% Generate the preamble signal
preambleSigSTS = helperGenPreamble(prm);
%   repeat over numTx
preambleSig = zeros(size(preambleSigSTS,1),numTx);
for i = 1:numSTS
    preambleSig(:,(i-1)*expFactorTx+(1:expFactorTx)) = ...
        repmat(preambleSigSTS(:,i),1,expFactorTx);
end

% Transmit preamble over channel
[rxPreSig,chanDelay] = helperApplyChannel(preambleSig,prm,spLoss);

% Front-end amplifier gain and thermal noise
rxPreAmp = phased.ReceiverPreamp( ...
    'Gain',gainFactor*spLoss, ... % account for path loss
    'NoiseFigure',prm.NFig, ...
    'ReferenceTemperature',290, ...
    'SampleRate',prm.chanSRate);
rxPreSigAmp = rxPreAmp(rxPreSig);
rxPreSigAmp = rxPreSigAmp * ...         % scale power
    (sqrt(prm.FFTLength-sum(prm.NumGuardBandCarriers)-1)/(prm.FFTLength));  

% OFDM Demodulation
demodulatorOFDM = comm.OFDMDemodulator( ...
     'FFTLength',prm.FFTLength, ...
     'NumGuardBandCarriers',prm.NumGuardBandCarriers.', ...
     'RemoveDCCarrier',true, ...
     'PilotOutputPort',true, ...
     'PilotCarrierIndices',prm.PilotCarrierIndices.', ...
     'CyclicPrefixLength',prm.CyclicPrefixLength, ...
     'NumSymbols',numSTS, ... % preamble symbols alone
     'NumReceiveAntennas',numRx);

rxOFDM = demodulatorOFDM( ...
    rxPreSigAmp(chanDelay+1:end-(prm.numPadZeros-chanDelay),:));

hD = helperMIMOChannelEstimate(rxOFDM(:,1:numSTS,:),prm); 

% Calculate the feedback weights
v = diagbfweights(hD);

%% 

encoder = comm.ConvolutionalEncoder( ...
    'TrellisStructure',poly2trellis(7,[133 171 165]), ...
    'TerminationMethod','Terminated');

% Generate mapped symbols from bits
txBits = randi([0, 1],prm.numFrmBits,1);
encodedBits = encoder(txBits);

% Bits to QAM symbol mapping
mappedSym = qammod(encodedBits,prm.modMode,'InputType','Bit', ...
    'UnitAveragePower',true);

gridData = reshape(mappedSym,prm.numCarriers,prm.numDataSymbols,numSTS);

preData = complex(zeros(prm.numCarriers,prm.numDataSymbols,numSTS));
for symIdx = 1:prm.numDataSymbols
    for carrIdx = 1:prm.numCarriers
        Q = squeeze(v(carrIdx,:,:));
        normQ = Q * sqrt(numTx)/norm(Q,'fro');      
        preData(carrIdx,symIdx,:) = ...
            squeeze(gridData(carrIdx,symIdx,:)).' * normQ;
    end
end

% OFDM modulation of the data
modulatorOFDM = comm.OFDMModulator( ...
    'FFTLength',prm.FFTLength,...
    'NumGuardBandCarriers',prm.NumGuardBandCarriers.',...
    'InsertDCNull',true, ...
    'PilotInputPort',true,...
    'PilotCarrierIndices',prm.PilotCarrierIndices.',...
    'CyclicPrefixLength',prm.CyclicPrefixLength,...
    'NumSymbols',prm.numDataSymbols,...
    'NumTransmitAntennas',numSTS);

% Multi-antenna 
pilots = helperGenPilots(prm.numDataSymbols,numSTS);

txOFDM = modulatorOFDM(preData,pilots);
txOFDM = txOFDM * (prm.FFTLength/ ...
    sqrt(prm.FFTLength-sum(prm.NumGuardBandCarriers)-1)); % scale power

preambleSigD = helperGenPreamble(prm,v);
txSigSTS = [preambleSigD;txOFDM];

txSig = zeros(size(txSigSTS,1),numTx);
for i = 1:numSTS
    txSig(:,(i-1)*expFactorTx+(1:expFactorTx)) = ...
        repmat(txSigSTS(:,i),1,expFactorTx);
end

%%

% Gain per antenna element 
amplifier = phased.Transmitter('PeakPower',1/numTx,'Gain',0);

% Amplify to achieve peak transmit power for each element
for n = 1:numTx
    txSig(:,n) = amplifier(txSig(:,n));
end

% Transmit antenna array definition 
if isTxURA
    % Uniform Rectangular array
    arrayTx = phased.URA([expFactorTx,numSTS],[0.5 0.5]*prm.lambda, ...
        'Element',phased.IsotropicAntennaElement('BackBaffled',true));
else
    % Uniform Linear array
    arrayTx = phased.ULA(numTx, ...
        'ElementSpacing',0.5*prm.lambda, ...
        'Element',phased.IsotropicAntennaElement('BackBaffled',true));
end

% For evaluating weights for steering  
SteerVecTx = phased.SteeringVector('SensorArray',arrayTx, ...
    'PropagationSpeed',prm.cLight);

% Generate weights for steered direction
wT = SteerVecTx(prm.fc,prm.steeringAngle);

% Radiate along the steered direction, without signal combining
radiatorTx = phased.Radiator('Sensor',arrayTx, ...
    'WeightsInputPort',true, ...
    'PropagationSpeed',prm.cLight, ...
    'OperatingFrequency',prm.fc, ...
    'CombineRadiatedSignals',false);

if prm.enSteering
    txSteerSig = radiatorTx(txSig,repmat(prm.mobileAngle,1,numTx), ...
        conj(wT));
else
    txSteerSig = txSig;
end

%% AWGN channel applying


rxPreAmp = phased.ReceiverPreamp( ...
    'Gain',gainFactor*spLoss, ... % accounts for path loss
    'NoiseFigure',prm.NFig, ...
    'ReferenceTemperature',290, ...
    'SampleRate',prm.chanSRate);

% Front-end amplifier gain and thermal noise
rxSigAmp = rxPreAmp(rxSig);
rxSigAmp = rxSigAmp * ...           % scale power
    (sqrt(prm.FFTLength - sum(prm.NumGuardBandCarriers)-1)/(prm.FFTLength)); 

% Receive array
if isRxURA 
    arrayRx = phased.URA([expFactorRx,numSTS],0.5*prm.lambda, ...
        'Element',phased.IsotropicAntennaElement('BackBaffled',true));
else 
    arrayRx = phased.ULA(numRx, ...
        'ElementSpacing',0.5*prm.lambda, ...
        'Element',phased.IsotropicAntennaElement);
end

SteerVecRx = phased.SteeringVector('SensorArray',arrayRx, ...
    'PropagationSpeed',prm.cLight);

wR = SteerVecRx(prm.fc,toRxAng);

if prm.enSteering
    rxSteerSig = rxSigAmp.*(wR');
else
    rxSteerSig = rxSigAmp;
end

%%

demodulatorOFDM = comm.OFDMDemodulator( ...
     'FFTLength',prm.FFTLength, ...
     'NumGuardBandCarriers',prm.NumGuardBandCarriers.', ...
     'RemoveDCCarrier',true, ...
     'PilotOutputPort',true, ...
     'PilotCarrierIndices',prm.PilotCarrierIndices.', ...
     'CyclicPrefixLength',prm.CyclicPrefixLength, ...
     'NumSymbols',numSTS+prm.numDataSymbols, ... % preamble & data
     'NumReceiveAntennas',numRx);
  
% OFDM Demodulation
rxOFDM = demodulatorOFDM( ...
    rxSteerSig(chanDelay+1:end-(prm.numPadZeros-chanDelay),:));

hD = helperMIMOChannelEstimate(rxOFDM(:,1:numSTS,:),prm);

[rxEq,CSI] = helperMIMOEqualize(rxOFDM(:,numSTS+1:end,:),hD);

scFact = ((prm.FFTLength-sum(prm.NumGuardBandCarriers)-1) ...
         /prm.FFTLength^2)/numTx;
nVar = noisepow(prm.chanSRate,prm.NFig,290)/scFact;
rxSymbs = rxEq(:)/sqrt(numTx);
rxLLRBits = qamdemod(rxSymbs,prm.modMode,'UnitAveragePower',true, ...
    'OutputType','approxllr','NoiseVariance',nVar);

rxLLRtmp = reshape(rxLLRBits,prm.bitsPerSubCarrier,[], ...
                   prm.numDataSymbols,numSTS);
csitmp = reshape(CSI,1,[],1,numSTS);
rxScaledLLR = rxLLRtmp.*csitmp;

decoder = comm.ViterbiDecoder(...
     'InputFormat','Unquantized', ...
     'TrellisStructure',poly2trellis(7, [133 171 165]), ...
     'TerminationMethod','Terminated', ...
     'OutputDataType','double');
rxDecoded = decoder(rxScaledLLR(:));

% Decoded received bits
rxBits = rxDecoded(1:prm.numFrmBits);

%% 

constDiag = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1, ...
    'ShowReferenceConstellation',true, ...
    'ReferenceConstellation', ...
    qammod((0:prm.modMode-1)',prm.modMode,'UnitAveragePower',true), ...
    'ColorFading',false, ...
    'Position',figposition([20 20 35 40]), ...
    'Title','Equalized Symbols', ...
    'EnableMeasurements',true, ...
    'MeasurementInterval',length(rxSymbs));

% Compute and display bit error rate
ber = comm.ErrorRate;
measures = ber(txBits,rxBits);
fprintf('BER = %.5f; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1),measures(3),measures(2));

rng(s); 
[noerr,suc_rate] = biterr(txBits,rxBits);
    
success_rateVEC(1,ss)=suc_rate
  end
  success_rate_v=success_rateVEC/max(success_rateVEC);
  
  command=1;
result=plotfunc(noise_pow/100,success_rate_v,command);


%% plot ===================================================================
result=plotfunc(SNR_vec/length(SNR_vec),success_rate,command);

