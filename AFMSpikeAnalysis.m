clear
clc
close all
set ( 0 , 'DefaultFigureWindowStyle' , 'Docked' ) ;

%% More quickly process pre-sectioned segments of CAFM 0D (time-domain) data.
% Files should be sections of data cut from full measurements,
% corresponding to regions of constant stress (i.e. constant current or
% voltage, or other data would also work). This is for spike analysis, i.e.
% during stress, spiking occurs, and this script allows the user to look
% at the distributions of spiking events and perform FFT. Typical input
% data is absolute, sensitivity-corrected current/voltage/time floating
% points, but the script should readily be adjusted for other data types.
% As written, input files should be 3 ascii columns without a header;
% 1) time; 2) current; 3) voltage. Can process multiple files at once rather
% than altogether. If files names contain a metric denoting the test
% condition (using the same unit magnitude for each title) then they
% should get sorted accordingly.

%% ADDITIONAL MATLAB FUNCTIONS REQUIRED.
% distinguishable_colors.m, unless this colour set isn't used. It's not
% essential.

%% Get files and locations.
addpath ( cd ) ;
[ FileGroup , DataPath ] = uigetfile ( '*.txt' , 'DialogTitle' ,...
    'Select files:' , 'MultiSelect' , 'on' ) ; % Gets file names and location.
% Used to account for the case where only a single file is selected.
SingleFile = double ( ischar ( FileGroup ) ) ;
if  SingleFile > 0
    NFC = 1 ;
    FileSet = strcat ( DataPath , FileGroup ) ;
    % Chosen file to determine number of columns, set to 1 to use first specified file.
    FileChoice = fopen ( FileSet , 'r' ) ;
else
    NFC = length ( FileGroup ) ;  % Number of files to import.
    FileSet = repmat ( { '' } , 1 , NFC ) ; % Generates cell array to place filenames into.
% Generates array of files for analysis.
    for i = 1 : NFC
        % Concatenates path and file strings and adds to output array.
        FileSet ( i ) = strcat ( DataPath , FileGroup ( i ) ) ;
    end
    % Chosen file to determine number of columns, set to 1 to use first specified file.
    FileChoice = fopen ( char ( FileSet ( : , 1 ) ) , 'r' ) ;
end
cd ( DataPath ) ;

%% Import data.
DataArray = repmat ( { ' ' } , NFC , 1 ) ; % Recipient array for all data.
SetPoints = zeros ( NFC , 1 ) ; % Setpoint values.
SampleFrequency = zeros ( NFC , 1 ) ; % Matrix for sample lengths per file.
AnalysisSection = repmat ( { ' ' } , NFC , 1 ) ; % For cropped analysis sections.
VIterationRange = [ 0.05 0.1 0.2 0.5 ]; % Range of spike voltage prominence thresholds.
NumProms = numel ( VIterationRange ) ; % Number of prominences used.
for i = 1 : NFC
    if NFC == 1 
        FileName = char ( FileSet ) ; % Choose single file.
        SetPoint = textscan ( FileGroup { i } , '%s%f%s' ) ;
        SetPoints ( i ) = SetPoint { 2 } ; % Get setpoint.
    else
        FileName = char ( FileSet ( : , i ) ) ; % Choose file from set.
        SetPoint = textscan ( FileGroup { i } , '%s%f%s' ) ;
        SetPoints ( i ) = SetPoint { 2 } ; % Get setpoint.
    end
    % Read in data from files.
    FileChoice = fopen ( FileName , 'r' ) ;
	DataArray { i } = textscan ( FileChoice , '%f%f%f' , 'Delimiter' , '' ) ;
    % Offset time data to start at 0 s.
    DataArray { i } { 1 } = DataArray { i } { 1 } - DataArray { i } { 1 } ( 1 ) ;
    fclose ( FileChoice ) ;
    % Produce smoothed data (background) and filtered, rectified data.
    SampleFrequency ( i ) = numel ( DataArray { i } { 1 } ) / DataArray { i } { 1 } ( end ) ;
    FilterWidth = round ( SampleFrequency ( i ) ) ;
    DataArray { i } { 4 } = smooth ( DataArray { i } { 3 } , FilterWidth , 'sgolay' , 1 ) ;
    DataArray { i } { 5 } = abs ( highpass ( DataArray { i } { 3 } , 2 , SampleFrequency ( i ) ) ) ;
    % Produce section for analysis, cropping first and last 10 seconds.
    AnalysisSection { i } = DataArray { i } { 3 } ( 10 * FilterWidth : end -...
        10 * FilterWidth ) ;
end 
% Higher resolution data range to use for fitting, if necessary.
VFittingRange = SetPoints ( 1 ) : 0.01 : SetPoints ( end ) ;
%% Process data to get ditributions of peak parameters.
OutputData = zeros (NumProms , 17 , NFC ) ;
% Data columns in output array:
% 1) spikes per second (SpS), 2) modal spike separation (max in histogram)
% 3) frequency of SpS (from histogram), 4) SpS distribution µ, 5) SpS distribution ?.
% 6) modal spike width (max in histogram), 7) frequency of max width (from
% histogram), 8) width distribution µ, 9) width distribution ?, 10) modal
% spike prominence (max in histogram), 11) frequency of max prominence
% (from histogram), 12) prominence µ, 13) prominence ?, 14) modal spike
% height (max in histogram), 15) frequency of max height (from histogram),
% 16) height µ, 17) height ?.

MinPeakHeight = 0 ; % Minimum absolute peak height in V.
MinPeakWidth = 0 ; % Minimum peak width threshold in s.
MaxPeakWidth = 5 ; % Maximum peak width threshold in s.
DistType = 'LogNormal' ;  % Type of distribution to fit.
% For lognormal, fitting parameters µ (median) and ? (scatter) are used.
for i = 1 : NFC
    % Find spikes (peaks) in the data. This section looks through each file
    % and iterates findpeaks.m through a range of spike prominences, defined
    % above as VIterationRange. At each prominence, a histogram of
    % the resulting spike spacings, widths, prominences and heights are
    % generated. From this, the modal value (i.e. histogram max) is taken,
    % along with the frequency at this value and the fitting parameters.
    % These are then output and plotted to look at the change in
    % distributions as a function of the applied bias.
    for j = 1 : NumProms
        [ VoltagePeaks , VoltageLocs , VoltageWidths , VoltageProms ] = ...
            findpeaks ( DataArray { i } { 3 } , SampleFrequency ( i ) , ...
            'MinPeakProminence' , VIterationRange ( j ) , 'MinPeakHeight' ,...
            MinPeakHeight , 'WidthReference' , 'HalfProm' , 'MinPeakWidth' ,...
            MinPeakWidth , 'MaxPeakWidth' , MaxPeakWidth ) ;
        % Get distributions of values.
        % Number of peaks per second.
        OutputData ( j , 1 , i ) = numel ( VoltagePeaks ) / DataArray { i } { 1 } ( end ) ;
        %% Spike separation.
        VPeakSeps = diff ( VoltageLocs ) ;
        if numel ( VPeakSeps ) < 2
            % Use length of sampling region as time between spikes in the
            % case that less than two spikes are detected.
            OutputData ( j , 2 , i ) = DataArray { i } { 1 } ( end ) ;
            OutputData ( j , 3 , i ) = 0 ;
            OutputData ( j , 4 , i ) = 0 ;
            OutputData ( j , 5 , i ) = 0 ;
        else
            VPeakSeps = diff ( VoltageLocs ) ;
            SepBins = round ( max ( VPeakSeps ) / ( min ( VPeakSeps ) / 2 ) ) ;
            % Ensure there aren't so many bins that Matlab can't handle it.
            if SepBins > 512
                SepBins = 512 ;
            end
            SepHist = histfit ( VPeakSeps , SepBins , DistType ) ;
            [ SMax , SIndex ] = max ( SepHist ( 1 ) . YData ) ;
            OutputData ( j , 2 , i ) = SepHist ( 1 ) . XData ( SIndex ) ;
            OutputData ( j , 3 , i ) = SMax ;
            SepDist = fitdist ( VPeakSeps , DistType ) ;
            OutputData ( j , 4 , i ) = SepDist . mu ;
            OutputData ( j , 5 , i ) = SepDist . sigma ;
            xlabel ( 'Peak separation/s' ) ;
            ylabel ( 'Counts' ) ;
            title ( strcat ( FileGroup { i } ( 1 : end - 4 ) , { ' at ' } ,...
                num2str ( VIterationRange ( j ) ) , ' V prominence' ) ) ;
            set ( gca , 'FontSize' , 16 ) ;
            set ( gcf , 'Color' , 'w' ) ;
            % Histograms can be plotted and/or saved.
%             saveas ( gcf , char ( strcat ( 'SpikeSepHist,' , { ' ' } ,...
%                 num2str ( VPromThresh ) , ' V prominence.tif' ) ) ) ;
        end
        %% Spike widths.
        if numel ( VoltageWidths ) < 2
            OutputData ( j , 6 , i ) = 0 ;
            OutputData ( j , 7 , i ) = 0 ;
            OutputData ( j , 8 , i ) = 0 ;
            OutputData ( j , 9 , i ) = 0 ;
        else
            WidthBins = round ( max ( VoltageWidths ) / ( min ( VoltageWidths ) / 2 ) ) ;
            if WidthBins > 512
                WidthBins = 512 ;
            end
            WidthHist = histfit ( VoltageWidths , WidthBins , DistType ) ;
            [ WMax , WIndex ] = max ( WidthHist ( 1 ) . YData ) ;
            OutputData ( j , 6 , i ) = WidthHist ( 1 ) . XData ( WIndex ) ;
            OutputData ( j , 7 , i ) = WMax ;
            WidthDist =  fitdist ( VoltageWidths , DistType )  ;
            OutputData ( j , 8 , i ) = WidthDist . mu ;
            OutputData ( j , 9 , i ) = WidthDist . sigma ;
            xlabel ( 'Peak width/s' ) ;
            ylabel ( 'Counts' ) ;
            title ( strcat ( FileGroup { i } ( 1 : end - 4 ) , { ' at ' } ,...
                num2str ( VIterationRange ( j ) ) , ' V prominence' ) ) ;
            set ( gca , 'FontSize' , 16 ) ;
            set ( gcf , 'Color' , 'w' ) ;            
%             saveas ( gcf , char ( strcat ( 'SpikeWidthHist,' , { ' ' } ,...
%                 num2str ( VPromThresh ) , ' V prominence.tif' ) ) ) ;
        end
        %% Spike prominences.
        if numel ( VoltageProms ) < 2
            OutputData ( j , 10 , i ) = 0 ;
            OutputData ( j , 11 , i ) = 0 ;
            OutputData ( j , 12 , i ) = 0 ;
            OutputData ( j , 13 , i ) = 0 ;
        else
            PromsBins = round ( max ( VoltageProms ) / ( min ( VoltageProms ) / 2 ) ) ;
            if PromsBins > 512
                PromsBins = 512 ;
            end
            PromsHist = histfit ( VoltageProms , PromsBins , DistType ) ;
            [ PMax , PIndex ] = max ( PromsHist ( 1 ) . YData ) ;
            OutputData ( j , 10 , i ) = PromsHist ( 1 ) . XData ( PIndex ) ;
            OutputData ( j , 11 , i ) = PMax ;
            PromsDist =  fitdist ( VoltageProms , DistType )  ;
            OutputData ( j , 12 , i ) = PromsDist . mu ;
            OutputData ( j , 13 , i ) = PromsDist . sigma ;
        end
        %% Spike heights.
        if numel ( VoltagePeaks ) < 2
            OutputData ( j , 14 , i ) = 0 ;
            OutputData ( j , 15 , i ) = 0 ;
            OutputData ( j , 16 , i ) = 0 ;
            OutputData ( j , 17 , i ) = 0 ;
        else
            HeightsBins = round ( max ( VoltagePeaks ) / ( min ( VoltagePeaks ) / 2 ) ) ;
            if HeightsBins > 512
                HeightsBins = 512 ;
            end
            HeightsHist = histfit ( VoltagePeaks , HeightsBins , DistType) ;
            [ HMax , HIndex ] = max ( HeightsHist ( 1 ) . YData ) ;
            OutputData ( j , 14 , i ) = HeightsHist ( 1 ) . XData ( HIndex ) ;
            OutputData ( j , 15 , i ) = HMax ;
            HeightsDist =  fitdist ( VoltagePeaks , DistType )  ;
            OutputData ( j , 16 , i ) = HeightsDist . mu ;
            OutputData ( j , 17 , i ) = HeightsDist . sigma ;
        end
    end
end
close ;

%% Plot and fit data.
% Option of rising edge threshold fitting algorithms.
iFit = 1 ;
if iFit == 1
    % Hyperbolic tangent, iFit = 1 (Mehonic, Front. Neurosci., vol. 10, no.
    % 57, pp. 1?10, 2016.)
    RisingThresholdFit = fittype ( 'a*(1+(tanh((b*x)-c)))' , 'Independent' ,...
        'x' , 'Dependent' , 'y' ) ;
    SpSUpper = [ 10000 10000 10000 ] ;
    SpSLower = [ -10000 -10000 -10000 ] ;
    SpSTitle =  'SpS tanh fit = p_1 [ 1 + tanh ( p_2i_{set} - p_0 )]'  ;
elseif iFit == 2
    % Sigmoid, iFit = 2, probabilistic early deep neural network
    % "squishification" (3Blue1Brown, YouTube).
    RisingThresholdFit =  fittype ( '(a/(1+exp(b*-x)))+c' , 'Independent' ,...
        'x' , 'Dependent' , 'y' ) ;
    SpSUpper = [ 10000 10000 10000 ] ;
    SpSLower = [ -10000 -10000 -10000 ] ;
    SpSTitle =  'SpS sigmoid fit = [p_1 / (1 + e^{-p_2x})] + p_0'  ;
elseif iFit == 3
    % Softplus ReLU, iFit = 3 (rectified linear unit), more current
    % "squishifier" (3Blue1Brown, Youtube).
    RisingThresholdFit = fittype ( '(a*log(1+exp(b*x)))+c' , 'Independent' ,...
        'x' , 'Dependent' , 'y' ) ;
    SpSUpper = [ 10000 10000 10000 ] ;
    SpSLower = [ -10000 -10000 -10000 ] ;
    SpSTitle =  'SpS softplus ReLU fit = p_1ln(1+e^{p_2x}) + p_0' ;
end
FallingThresholdFit = fittype ( 'a*(1-(tanh((b*x)-c)))' , 'independent' ,...
    'x' , 'dependent' , 'y' ) ;
%% Number of spikes per second (from entire range).
figure ;
iColours = 3 ; % Choose colour set to plot with.
if iColours == 1 % Creates rainbow colour set.
    ChosenColours = jet (NumProms ) ;
elseif iColours == 2 % Creates the most distinct colour set.
    ChosenColours = distinguishable_colors ( NumProms ) ;
elseif iColours == 3 % Creates another colour set.
    ChosenColours = lines ( NumProms ) ;
end
SpS = zeros ( NumProms , 4 ) ;
SpSFitOptions = fitoptions ( 'Method' , 'NonlinearLeastSquares' ) ;
SpSFitOptions . Display = 'Off' ;
SpSFitOptions . Robust = 'LAR' ;
SpSFitOptions . StartPoint = [ 0 0 0 ] ;
SpSFitOptions . Upper = SpSUpper ;
SpSFitOptions . Lower = SpSLower ;
for i = 1 : NumProms
    SpikesPerSecond = reshape ( OutputData ( i , 1 , : ) , NFC , 1 ) ;
    semilogx ( SetPoints , SpikesPerSecond , 'o' , 'Color' ,...
        ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spikes per second' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    [ SpSFit , SpSGoF ] = fit ( SetPoints , SpikesPerSecond , RisingThresholdFit ,...
        SpSFitOptions ) ;
    %% For all fits.
    SpS ( i , 1 ) = SpSFit . a ;
    SpS ( i , 2 ) = SpSFit . b ;
    SpS ( i , 3 ) = SpSFit . c ;
    SpS ( i , end ) = SpSGoF . rsquare ;
    % For tanh fits - this is a really good fit.
    if iFit == 1
        SpSModel = SpS ( i , 1 ) .* ( 1 + tanh ( ( SpS ( i , 2 ) .* VFittingRange ) -...
            SpS ( i , 3 ) ) ) ;
    % For sigmoid fits - this is also a really good fit.
    elseif iFit == 2
        SpSModel = ( SpS ( i , 1 ) ./ ( 1 + exp ( SpS ( i , 2 ) .* ( -1 ) .*...
            VFittingRange ) ) ) + SpS ( i , 3 ) ;
    % For softplus ReLU fits - this is a poor fit, just increases rapidly with
    % no plateau.
    elseif iFit == 3
        SpSModel = ( SpS ( i , 1 ) .* log ( 1 + exp ( ( SpS ( i , 2 ) .*...
        VFittingRange ) ) ) ) + SpS ( i , 3 ) ;
    end
    semilogx ( VFittingRange , SpSModel  , 'Color' , ChosenColours ( i , : ) ) ;
    %% Also plot change in fitting parameters at end of loop.
    if i == NumProms
        figure ;
        plot ( VIterationRange , SpS ( : , 1 ) / SpS ( 1 , 1 )  ,...
            VIterationRange , SpS ( : , 2 ) / SpS ( 1 , 2 ) ,...
            VIterationRange , SpS ( : , 3 ) / SpS ( 1 , 3 ) ,...
            VIterationRange , SpS ( : , 4 ) , 'LineWidth' , 2 ) ;
        xlabel ( 'Prominence threshold/V' ) ;
        ylabel ( 'Normalized parameter value' ) ;
        legend ( 'p_1' , 'p_2' , 'p_0' , 'Goodness of fit' , 'Location' , 'NorthWest' ) ;
        legend ( 'BoxOff' ) ;
        set ( gca , 'FontSize' , 14 ) ;
        set ( gcf , 'Color' , 'w' ) ;
        title ( SpSTitle ) ;
        % p1 decreases power/logarithmically with threshold, just
        % indicating fewer spikes per second at higher thresholds. p2
        % increases linearly, which follows because the current se. p0
        % doesn't change much and likely just represents the poor fitting
        % at very low current bias.
    end
end
%% Modal time between spikes.
figure ;
TbS = zeros ( NumProms , 4 ) ;
TbSFitOptions = fitoptions ( 'Method' , 'NonlinearLeastSquares' ) ;
TbSFitOptions . Display = 'Off' ;
TbSFitOptions . Robust = 'LAR' ;
for i = 1 : NumProms
    TimeBetweenSpikes = reshape ( OutputData ( i , 2 , : ) , NFC , 1 ) ;
    loglog ( SetPoints , TimeBetweenSpikes , 'o' , 'Color' , ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal time between spikes/s' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    [ TbSFit , TbSGoF ] = fit ( SetPoints , TimeBetweenSpikes , 'Power2' ,...
         TbSFitOptions ) ;
    TbS ( i , 1 ) = TbSFit . a ;
    TbS ( i , 2 ) = TbSFit . b ;
    TbS ( i , 3 ) = TbSFit . c ;
    TbS ( i , 4 ) = TbSGoF . rsquare ;
    loglog ( VFittingRange , ( TbS ( i , 1 ) .* ( VFittingRange .^ TbS ( i , 2 ) ) ) +...
        TbS ( i , 3 ) , 'Color' , ChosenColours ( i , : ) ) ;
    % Also plot change in fitting parameters at end of loop.
    if i == NumProms
        figure ;
        plot ( VIterationRange , TbS ( : , 1 ) , VIterationRange , TbS ( : , 2 ) ,...
            'LineWidth' , 2 ) ;
        xlabel ( 'Prominence threshold/V' ) ;
        ylabel ( 'Parameter value' ) ;
        set ( gca , 'FontSize' , 14 ) ;
        set ( gcf , 'Color' , 'w' ) ;
        title ( 'TbS threshold fitting' ) ;
        legend ( 'a' , 'b' , 'Location' , 'NorthWest' ) ;
    end
end
%% Modal spike width.
figure ;
SW = zeros ( NumProms , 4 ) ;
SWFitOptions = fitoptions ( 'Method' , 'NonlinearLeastSquares' ) ;
SWFitOptions . Display = 'Off' ;
SWFitOptions . Robust = 'LAR' ;
SWFitOptions . StartPoint = [ 0 0 0 ] ;
for i = 1 : NumProms
    ModalSpikeWidths = reshape ( OutputData ( i , 6 , : ) , NFC , 1 ) ;
    loglog ( SetPoints , ModalSpikeWidths , 'o' , 'Color' , ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal spike width/s' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
	[ SWFit , SWFitGoF ] = fit ( SetPoints , ModalSpikeWidths , '(a*exp(-b*x))+c' ,...
         SWFitOptions ) ; % Need to choose the best fitting option here.
     % (a*log(1+exp(b*x)))+c, the softplus ReLU, works well.
     % (a*exp(-b*x))+c, a single exponential, works well.
     % (a/(1+exp(b*-x)))+c, a sigmoid, also works well.
    SW ( i , 1 ) = SWFit . a ;
    SW ( i , 2 ) = SWFit . b ;
    SW ( i , 3 ) = SWFit . c ;
    SW ( i , 4 ) = SWFitGoF . rsquare ;
    loglog ( VFittingRange , ( SW ( i , 1 ) .* exp ( ( -1 ) .* SW ( i , 2 ) .*...
        VFittingRange ) ) + SW ( i , 3 ) , 'Color' , ChosenColours ( i , : ) ) ;
end
%% Modal energy per spike, as spike width x current setpoint x prominence.
% Using prominence is an idealised case; this assumes that spikes aren't
% typically much greater than the desired promience, in which case the
% energy per SynOps would increase.
figure ;
EpS = zeros ( NumProms , 4 ) ;
EpSFitOptions = fitoptions ( 'Method' , 'NonLinearLeastSquares' ) ;
EpSFitOptions . Robust = 'LAR' ;
% EpSFitOptions . Lower = [ -100 -100 ] ;
% EpSFitOptions . Upper = [ 100 100 ] ;
MeanVoltage = zeros ( NFC , 1 ) ;
for i = 1 : NFC
    MeanVoltage ( i ) = mean ( DataArray { i } { 4 } ( : ) ) ;
end
for i = 1 : NumProms
    EnergyPerSpike = reshape ( OutputData ( i , 6 , : ) , NFC , 1 ) .*...
        SetPoints .* VIterationRange ( i ) ;
    loglog ( SetPoints , EnergyPerSpike , 'o' , 'Color' , ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Energy per SynOps/nJ' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
	[ EpSFit , EpSFitGoF ] = fit ( SetPoints , EnergyPerSpike , 'Power2' ,...
        EpSFitOptions ) ; % Linear fitting works best here.
    EpS ( i , 1 ) = EpSFit . a ;
    EpS ( i , 2 ) = EpSFit . b ;
    SpS ( i , 3 ) = EpSFit . c ;
    EpS ( i , 4 ) = EpSFitGoF . rsquare ;
    loglog ( VFittingRange , ( EpS ( i , 1 ) .* ( VFittingRange  .^ EpS ( i , 2 ) ) ) +...
        EpS ( i , 3 ) , 'Color' , ChosenColours ( i , : ) ) ;
end
%% Spike separation distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 4 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike separation µ/-' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 5 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , ChosenColours ( i , : ) ) ;
        ylabel ( 'Spike separation {\sigma}/--' ) ;
end
%% Spike width distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 8 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike width µ/-' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 9 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , ChosenColours ( i , : ) ) ;
        ylabel ( 'Spike width {\sigma}/--' ) ;
end
%% Spike prominence distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 12 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike prominence µ/-' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 13 , : ) , NFC , 1 ) ) , '--' ,....
        'Color' , ChosenColours ( i , : ) ) ;
        ylabel ( 'Spike prominence {\sigma}/--' ) ;
end
%% Modal spike heights.
figure ;
for i = 1 : NumProms
    ModalSpikeHeights = reshape ( OutputData ( i , 14 , : ) , NFC , 1 ) ;
    semilogx ( SetPoints , ModalSpikeHeights  , 'Color' ,...
        ChosenColours ( i , : ) ) ;
	xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal spike height/V' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
end
%% Spike height distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 16 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike height µ/-' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 17 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , ChosenColours ( i , : ) ) ;
        ylabel ( 'Spike height {\sigma}/--' ) ;
end
