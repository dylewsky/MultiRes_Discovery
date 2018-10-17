clear variables; close all; clc
%%
dataDir = '3 Body Example';
inFile = fullfile(dataDir,'mwDMD_sep_recon.mat');
load(inFile);
load(fullfile(dataDir,'mwDMD_params.mat'));
t = tspan;

% Orbital 
% mrComponent = 1;
% % Jupiter
% mrComponent = 2;
% % Saturn
% mrComponent = 3;

nDelays = [16 16 16];
delayStepSizes = [round(wSteps/nDelays(1)), 200, 200];
useLowpass = [1 0 0];

for mrComponent = 1:nComponents
    x = real(xr_sep{mrComponent});
    t = tspan;

    nVars = size(x,1);
    nDelay = nDelays(mrComponent);
    delaySteps = delayStepSizes(mrComponent);


    dt = t(2)-t(1);
  
    
    %% Low Pass Filter
    downSample = 1; %no downsampling
    if useLowpass(mrComponent) == 1
        filtTrunc = ceil([50000 0]/downSample); %additional steps to cut off from either end after filtering
        fcutlow=0.0002;   %low cut frequency in Hz
        fs = 1/(t(2)-t(1));
        [b,a]=butter(3,(fcutlow/fs)*2,'low');

        xBand = zeros(size(x));
        for j = 1:size(x,1)
            xBand(j,:) = filter(b, a, x(j,:));
        end

        filtShift = ceil(17943/downSample);
        % % Determine optimal filter shift
        % shiftIterations = 5;
        % testShifts = 1:5;
        % for shiftIter = 1:shiftIterations
        %     bestNorm = norm(yBand-y);
        %     optJ = 0;
        %     for j = 1:length(testShifts)
        %         thisShift = round(10^testShifts(j));
        %         testNorm = norm(y(:,1:end-thisShift) - yBand(:,thisShift+1:end));
        % %         figure
        % %         plot(y(:,1:end-thisShift).','k')
        % %         hold on
        % %         plot(yBand(:,thisShift+1:end).','r')
        % %         hold off
        % %         disp([shiftIter, j])
        %         if testNorm < bestNorm
        %             optJ = j;
        %             bestNorm = testNorm;
        %         end
        %     end
        %     if shiftIter < shiftIterations
        %         testShiftStep = testShifts(2)-testShifts(1);
        %         testShifts = linspace(testShifts(optJ)-testShiftStep/2, testShifts(optJ)+testShiftStep/2, length(testShifts));
        %     end
        % end
        % filtShift = round(10^(testShifts(optJ)));
        % clear('testShifts','testNorm','bestNorm','shiftIter','shiftIterations','testShiftStep');

        xBand = xBand(:,filtShift+1:end);
        t = t(1:end-filtShift);
        % figure
        % plot(tspan,y(:,1:end-filtShift).','k')
        % hold on
        % plot(tspan,yBand.','r')
        % hold off
        x = xBand; clear('yBand');

        x = x(:,filtTrunc(1)+1:end-filtTrunc(2));
        t = t(filtTrunc(1)+1:end-filtTrunc(2));
    end

    %% EIGEN-TIME DELAY COORDINATES
    stackmax = nVars*nDelay; % Number of shift-stacked rows
    r=16; % Rank of HAVOK Model
    H = zeros(stackmax,size(x,2)-(nDelay-1)*delaySteps);
    for k=1:nDelay
        delayInds = ((k-1)*nVars + 1) : k*nVars;
        H(delayInds,:) = x(:,(k-1)*delaySteps+1:(size(x,2)-(nDelay)*delaySteps + k*delaySteps));
    end

    [U,S,V] = svd(H,'econ'); % Eigen delay coordinates

    %% COMPUTE DERIVATIVES (4TH ORDER CENTRAL DIFFERENCE)
    dV = zeros(length(V)-5,r);
    for i=3:length(V)-3
        for k=1:r
            dV(i-2,k) = (1/(12 * dt)) * (-V(i+2,k)+8 * V(i+1,k)-8 * V(i-1,k)+V(i-2,k));
        end
    end
    % trim first and last two that are lost in derivative
    V = V(3:end-3,1:r);

    %% BUILD HAVOK REGRESSION MODEL ON TIME DELAY COORDINATES
    Xi = V\dV;
    A = Xi(1:r-1,1:r-1)';
    B = Xi(end,1:r-1)';

    %% FIGURES
    figure('Position',[100 75 600 900])
    subplot(3,4,1:3)
    imagesc(A)
    caxis([min(min(Xi)) max(max(Xi))]);
    title('A Matrix')

    subplot(3,4,4)
    imagesc(B)
    colorbar
    caxis([min(min(Xi)) max(max(Xi))]);
    title('B Vector')

    subplot(3,4,5:8)
    plot(t,x)
    title('Input Data')
    xlabel('t')
    ylabel('x')

    subplot(3,4,9:12)
    plot(t(1:size(V,1)),V(:,end))
    title('Control Input (v_r)')
    xlabel('t')
    ylabel('v_r')

    %% Save Model
    outFile = ['MR_HAVOK_model_' num2str(mrComponent) '.mat'];
    lowpass = useLowpass(mrComponent);
    U = U(:,1:r);
    S = S(1:r,1:r);
    save(fullfile(dataDir,outFile),'A','B','nDelay','delaySteps','dt','r','U','S','lowpass');
end