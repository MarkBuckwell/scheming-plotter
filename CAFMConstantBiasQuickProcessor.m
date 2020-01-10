clear
clc
close all
set ( 0 , 'DefaultFigureWindowStyle' , 'Docked' ) ;

%% More quickly process pre-sections constant bias CAFM 0D data. Files
% previously exported from matlab, as sections cut from full measurements,
% ready for analysis via FFT and spike distributions. Data is already
% absolute and has been sensitivity-corrected. Input files are 3 ascii
% columns without a header; 1) time; 2) current; 3) voltage. Process
% multiple files at once rather than altogether. If files are named such
% that the current bias of each measurement is in the title (using the
% same unit magnitude for each title) then they should get sorted accordingly.

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
VIterationRange = 0.05 : 0.05 : 1 ; % Range of spike voltage thresholds.
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

%% Process data to get ditributions of peak parameters.
OutputData = zeros ( numel ( VIterationRange ) , 17 , NFC ) ;
MinPeakHeight = 0 ; % Minimum absolute peak height in V.
MinPeakWidth = 0 ; % Minimum peak width threshold in s.
MaxPeakWidth = 5 ; % Maximum peak width threshold in s.
DistType = 'LogNormal' ;  % Type of distribution to fit.
for i = 1 : NFC
    % Find peaks in the data.
    for j = 1 : numel ( VIterationRange )
        [ VoltagePeaks , VoltageLocs , VoltageWidths , VoltageProms ] = ...
            findpeaks ( DataArray { i } { 3 } , SampleFrequency ( i ) , ...
            'MinPeakProminence' , VIterationRange ( j ) , 'MinPeakHeight' , MinPeakHeight ,...
            'WidthReference' , 'HalfProm' , 'MinPeakWidth' , MinPeakWidth ,...
            'MaxPeakWidth' , MaxPeakWidth ) ;
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
            OutputData ( i , 14 ) = 0 ;
            OutputData ( i , 15 ) = 0 ;
            OutputData ( i , 16 ) = 0 ;
            OutputData ( i , 17 ) = 0 ;
        else
            HeightsBins = round ( max ( VoltagePeaks ) / ( min ( VoltagePeaks ) / 2 ) ) ;
            if HeightsBins > 512
                HeightsBins = 512 ;
            end
            HeightsHist = histfit ( VoltagePeaks , HeightsBins , DistType) ;
            [ HMax , HIndex ] = max ( HeightsHist ( 1 ) . YData ) ;
            OutputData ( i , 14 ) = HeightsHist ( 1 ) . XData ( HIndex ) ;
            OutputData ( i , 15 ) = HMax ;
            HeightsDist =  fitdist ( VoltagePeaks , DistType )  ;
            OutputData ( i , 16 ) = HeightsDist . mu ;
            OutputData ( i , 17 ) = HeightsDist . sigma ;
        end
    end
end
close ;

%% Plot and fit data.
JetSet = jet ( numel ( VIterationRange ) ) ; % Creates rainbow colour set.

%% Option of rising edge threshold fitting algorithms.
iFit = 1 ;
% Hyperbolic tangent, iFit = 1 (Mehonic, Kenyon, Frontiers).
RisingTanH = fittype ( 'a*(1+(tanh((b*x)-c)))' , 'Independent' , 'x' , 'Dependent' , 'y' ) ;
% Sigmoid, iFit = 2, probabilistic early deep neural network "squishification" (3Blue1Brown, YouTube).
RisingSigmoid = fittype ( '(a/(1+exp(b*-x)))+c' , 'Independent' , 'x' , 'Dependent' , 'y' );
% Softplus ReLU, iFit = 3 (rectified linear unit), more current "squishifier" (3Blue1Brown, Youtube).
RisingSoftReLU = fittype ( '(a*log(1+exp(b*x)))+c' , 'Independent' , 'x' , 'Dependent' , 'y' ) ;

FallingThresholdFit = fittype ( 'a*(1-(tanh((b*x)-c)))' , 'independent' , 'x' , 'dependent' , 'y' ) ;
%% Number of spikes per second (from entire range).
figure ;
SpS = zeros ( numel ( VIterationRange ) , iFit + 1 ) ;
% Threhsold fitting as in ?Fig 7, A. Mehonic and A. J. Kenyon, ?Emulating the
% electrical activity of the neuron using a silicon oxide RRAM cell,? Front.
% Neurosci., vol. 10, no. 57, pp. 1?10, 2016.
for i = 1 : numel ( VIterationRange )
    SpikesPerSecond = reshape ( OutputData ( i , 1 , : ) , NFC , 1 ) ;
    semilogx ( SetPoints , SpikesPerSecond , 'o' , 'Color' , JetSet ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spikes per second' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    [ SpSFit , SpSGoF ] = fit ( SetPoints , SpikesPerSecond , RisingThresholdFit ,...
        'StartPoint' , [ 0 0 0 ] ) ;
    
    %% For all fits
    SpS ( i , 1 ) = SpSFit . a ;
    SpS ( i , 2 ) = SpSFit . b ;
    SpS ( i , 3 ) = SpSFit . c ;
    SpS ( i , end ) = SpSGoF . rsquare ;
    % For tanh fits - this is a really good fit.
    if iFit == 1
        SpSModel = SpS ( i , 1 ) .* ( 1 + tanh ( ( SpS ( i , 2 ) .* SetPoints ) -...
            SpS ( i , 3 ) ) ) ;
    % For sigmoid fits - this is a poor fit, doesn't capture the range
    % properly.
    elseif iFit == 2
        SpSModel = ( SpS ( i , 1 ) ./ ( 1 + exp ( SpS ( i , 2 ) .* ( -1 ) .*...
            SetPoints ) ) ) + SpS ( i , 3 ) ;
    % For softplus ReLU fits - this is a poor fit, just increases rapidly with
    % no plateau.
    elseif iFit == 3
        SpSModel = ( SpS ( i , 1 ) .* log ( 1 + exp ( ( SpS ( i , 2 ) .*...
        SetPoints ) ) ) ) + SpS ( i , 3 ) ;
    end
        
    semilogx ( SetPoints , SpSModel  , 'Color' , JetSet ( i , : ) ) ;
    %% Also plot change in fitting parameters at end of loop.
    if i == numel ( VIterationRange )
        figure ;
        plot ( VIterationRange , SpS ( : , 1 ) / SpS ( 1 , 1 )  ,...
            VIterationRange , SpS ( : , 2 ) / SpS ( 1 , 2 ) ,...
            VIterationRange , SpS ( : , 3 ) / SpS ( 1 , 3 ) , 'LineWidth' , 2 ) ;
        xlabel ( 'Prominence threshold/V' ) ;
        ylabel ( 'Normalized parameter value' ) ;
        set ( gca , 'FontSize' , 14 ) ;
        set ( gcf , 'Color' , 'w' ) ;
        title ( 'SpS threshold fitting, SpS = p_1 [ 1 + tanh ( p_2i_{set} - p_0 )]' ) ;
        % p1 decreases power/logarithmically with threshold, just
        % indicating fewer spikes per second at higher thresholds. p2
        % increases linearly, which follows because the current se. p0
        % doesn't change much and likely just represents the poor fitting
        % at very low current bias.
        legend ( 'p_1' , 'p_2' , 'p_0' , 'Location' , 'NorthWest' ) ;
    end
end
%% Modal time between spikes.
figure ;
TbS = zeros ( numel ( VIterationRange ) , 4 ) ;
TbSFitOptions = fitoptions ( 'Method' , 'NonlinearLeastSquares' ) ;
TbSFitOptions . Display = 'Off' ;
TbSFitOptions . Robust = 'Bisquare' ;
for i = 1 : numel ( VIterationRange )
    TimeBetweenSpikes = reshape ( OutputData ( i , 2 , : ) , NFC , 1 ) ;
    loglog ( SetPoints , TimeBetweenSpikes , 'o' , 'Color' , JetSet ( i , : ) ) ;
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
    semilogx ( SetPoints , ( TbS ( i , 1 ) .* ( SetPoints .^ TbS ( i , 2 ) ) ) +...
        TbS ( i , 3 ) , 'Color' , JetSet ( i , : ) ) ;
    % Also plot change in fitting parameters at end of loop.
    if i == numel ( VIterationRange )
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
for i = 1 : numel ( VIterationRange )
    loglog ( SetPoints , ( reshape ( OutputData ( i , 6 , : ) , NFC , 1 ) ) ,...
        'Color' , JetSet ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal spike width/s' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
end
%% Modal energy per spike, as spike width x current setpoint x mean voltage.
figure ;
MeanVoltage = zeros ( NFC , 1 ) ;
for i = 1 : NFC
    MeanVoltage ( i ) = mean ( DataArray { i } { 4 } ( : ) ) ;
end
for i = 1 : numel ( VIterationRange )
    loglog ( SetPoints , ( reshape ( OutputData ( i , 6 , : ) , NFC , 1 ) ) .*...
        SetPoints .* VIterationRange ( i ) , 'Color' , JetSet ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Energy per SynOps/nJ' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
end
%% Spike separation distribution.
figure ;
for i = 1 : numel ( VIterationRange )
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 4 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , JetSet ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike separation µ/solid' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 5 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , JetSet ( i , : ) ) ;
        ylabel ( 'Spike separation sigma/dashed' ) ;
end
%% Spike width distribution.
figure ;
for i = 1 : numel ( VIterationRange )
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 8 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , JetSet ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike width µ/solid' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 9 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , JetSet ( i , : ) ) ;
        ylabel ( 'Spike width sigma/dashed' ) ;
end
%% Spike prominence distribution.
figure ;
for i = 1 : numel ( VIterationRange )
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 12 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , JetSet ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike prominence µ/solid' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 13 , : ) , NFC , 1 ) ) , '--' ,....
        'Color' , JetSet ( i , : ) ) ;
        ylabel ( 'Spike prominence sigma/--' ) ;
end
%% Spike height distribution.
figure ;
for i = 1 : numel ( VIterationRange )
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 16 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , JetSet ( i , : ) ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike height µ/solid' ) ;
    set ( gca , 'FontSize' , 14 ) ;
    set ( gcf , 'Color' , 'w' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 17 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , JetSet ( i , : ) ) ;
        ylabel ( 'Spike height sigma/dashed' ) ;
end
