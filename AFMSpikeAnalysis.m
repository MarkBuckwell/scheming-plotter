clear
clc
close all
set ( 0 , 'DefaultFigureWindowStyle' , 'Docked' ) ;

%% More quickly process pre-sectioned segments of CAFM 0D (time-domain) data.
% Files should be sections of data cut from full measurements,
% corresponding to regions of constant stress (i.e. constant current or
% voltage, or other data would also work). This is for spike analysis, i.e.
% during stress, spiking occurs, and this script allows the user to look
% at the distributions of spiking events. Typical input
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
% VIterationRange = 0 : 0.1 : 5 ;
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
    DataArray { i } { 5 } = highpass ( DataArray { i } { 3 } , 2 , SampleFrequency ( i ) ) ;
    DataArray { i } { 6 } = abs ( highpass ( DataArray { i } { 3 } , 2 , SampleFrequency ( i ) ) ) ;
    % Produce section for analysis, cropping first and last 10 seconds.
    AnalysisSection { i } = DataArray { i } { 3 } ( 10 * FilterWidth : end -...
        10 * FilterWidth ) ;
end 
% Higher resolution data range to use for fitting, if necessary.
VFittingRange = SetPoints ( 1 ) : 0.01 : SetPoints ( end ) ;

%% Process data to get ditributions of peak parameters.
OutputData = zeros ( NumProms , 19 , NFC ) ;
% Data columns in output array:
% 1) spikes per second (SpS), 2) modal spike separation (max in histogram)
% 3) frequency of SpS (from histogram), 4) SpS distribution µ, 5) SpS distribution sigma.
% 6) modal spike width (max in histogram), 7) mean spike width
% 8) frequency of max width (from histogram), 9) width distribution µ,
% 10) width distribution sigma, 11) modal spike prominence (max in histogram),
% 12) frequency of max prominence (from histogram), 13) prominence µ,
% 14) prominence sigma, 15) modal spike height (max in histogram),
% 16) frequency of max height (from histogram), 17) height µ,
% 18) height sigma, 1() mean spike height.
FFTOutput = repmat ( { ' ' } , NFC , NumProms ) ;
% Output array for converting spikes to binary, for cleaner FFT analysis.

FindIn = 3 ; % Data column to find peaks in.
MinPeakHeight = 0 ; % Minimum absolute peak height in V.
MinPeakWidth = 0 ; % Minimum peak width threshold in s.
MaxPeakWidth = 0.5 ; % Maximum peak width threshold in s.
DistType = 'LogNormal' ;  % Type of distribution to fit.
SaveHists = 0 ; % Histograms can be saved if SaveHists == 1.
if SaveHists == 1 % New folder to put histograms in, if not already present.
    if ~exist ( 'Spike Histograms' , 'File' ) == 0
        mkdir ( 'Spike Histograms' ) ;
    end
    cd ./'Spike Histograms'
end
% For lognormal, fitting parameters µ (median) and sigma (scatter) are used.
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
            findpeaks ( DataArray { i } { FindIn } , SampleFrequency ( i ) , ...
            'MinPeakProminence' , VIterationRange ( j ) , 'MinPeakHeight' ,...
            MinPeakHeight , 'WidthReference' , 'HalfProm' , 'MinPeakWidth' ,...
            MinPeakWidth , 'MaxPeakWidth' , MaxPeakWidth ) ;
        if SaveHists == 1
            findpeaks ( DataArray { i } { FindIn } , SampleFrequency ( i ) , ...
            'MinPeakProminence' , VIterationRange ( j ) , 'MinPeakHeight' ,...
            MinPeakHeight , 'WidthReference' , 'HalfProm' , 'MinPeakWidth' ,...
            MinPeakWidth , 'MaxPeakWidth' , MaxPeakWidth ) ; 
            title ( strcat ( FileGroup { i } ( 1 : end - 4 ) , { ' at ' } ,...
                num2str ( VIterationRange ( j ) ) , ' V_p' ) ) ;
            xlabel ( 'Time/s' ) ;
            ylabel ( 'Voltage/V' ) ;
            set ( gca , 'FontSize' , FontSize ) ;
            set ( gcf , 'Color' , 'w' ) ;
            saveas ( gcf , char ( strcat ( 'Spikesection,' ,  { ' ' } ,...
                FileGroup { i } ( 1 : end - 4 ) ,...
                num2str ( VIterationRange ( j ) ) , ' V prominence.tif' ) ) ) ;
        end
        % Get distributions of values. Binning of separations and widths is
        % applied such that the resolution of the bins is around half the
        % smallest difference between values. Binning of prominences and
        % heights is applied such that the bins are around 0.1 V.
        % Number of peaks per second.
        OutputData ( j , 1 , i ) = numel ( VoltagePeaks ) /...
            DataArray { i } { 1 } ( end ) ;
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
            if SaveHists == 1
                title ( strcat ( FileGroup { i } ( 1 : end - 4 ) , { ' at ' } ,...
                    num2str ( VIterationRange ( j ) ) , ' V_p' ) ) ;
                xlabel ( 'Spike separation/s' ) ;
                ylabel ( 'Counts' ) ;
                set ( gca , 'FontSize' , FontSize ) ;
                set ( gcf , 'Color' , 'w' ) ;
                saveas ( gcf , char ( strcat ( 'SepsHist,' ,  { ' ' } ,...
                    FileGroup { i } ( 1 : end - 4 ) ,...
                    num2str ( VIterationRange ( j ) ) , ' V prominence.tif' ) ) ) ;
            end
        end
        %% Spike widths.
        if numel ( VoltageWidths ) < 2
            OutputData ( j , 6 , i ) = 0 ;
            OutputData ( j , 7 , i ) = 0 ;
            OutputData ( j , 8 , i ) = 0 ;
            OutputData ( j , 9 , i ) = 0 ;
        else
            WidthBins = round ( max ( VoltageWidths ) /...
                ( min ( VoltageWidths ) / 2 ) ) ;
            if WidthBins > 512
                WidthBins = 512 ;
            end
            WidthHist = histfit ( VoltageWidths , WidthBins , DistType ) ;
            [ WMax , WIndex ] = max ( WidthHist ( 1 ) . YData ) ;
            OutputData ( j , 6 , i ) = WidthHist ( 1 ) . XData ( WIndex ) ;
            OutputData ( j , 7 , i ) = mean ( VoltageWidths ) ;
            OutputData ( j , 8 , i ) = WMax ;
            WidthDist =  fitdist ( VoltageWidths , DistType )  ;
            OutputData ( j , 9 , i ) = WidthDist . mu ;
            OutputData ( j , 10 , i ) = WidthDist . sigma ;
            if SaveHists == 1
            title ( strcat ( FileGroup { i } ( 1 : end - 4 ) , { ' at ' } ,...
                num2str ( VIterationRange ( j ) ) , ' V_p' ) ) ;
            xlabel ( 'Spike width/s' ) ;
            ylabel ( 'Counts' ) ;
            set ( gca , 'FontSize' , FontSize ) ;
            set ( gcf , 'Color' , 'w' ) ;
                saveas ( gcf , char ( strcat ( 'WidthHist,' , { ' ' } ,...
                    FileGroup { i } ( 1 : end - 4 ) ,...
                    num2str ( VIterationRange ( j ) ) , ' V prominence.tif' ) ) ) ;
            end
        end
        %% Spike prominences.
        if numel ( VoltageProms ) < 2
            OutputData ( j , 11 , i ) = 0 ;
            OutputData ( j , 12 , i ) = 0 ;
            OutputData ( j , 13 , i ) = 0 ;
            OutputData ( j , 14 , i ) = 0 ;
        else
            PromsBins = round ( max ( VoltageProms ) /...
                ( min ( VoltageProms ) ) / 0.1 ) ;
            if PromsBins > 512
                PromsBins = 512 ;
            end
            PromsHist = histfit ( VoltageProms , PromsBins , DistType ) ;
            [ PMax , PIndex ] = max ( PromsHist ( 1 ) . YData ) ;
            OutputData ( j , 11 , i ) = PromsHist ( 1 ) . XData ( PIndex ) ;
            OutputData ( j , 12 , i ) = PMax ;
            PromsDist =  fitdist ( VoltageProms , DistType )  ;
            OutputData ( j , 13 , i ) = PromsDist . mu ;
            OutputData ( j , 14 , i ) = PromsDist . sigma ;
            if SaveHists == 1
            title ( strcat ( FileGroup { i } ( 1 : end - 4 ) , { ' at ' } ,...
                num2str ( VIterationRange ( j ) ) , ' V_p' ) ) ;
            xlabel ( 'Spike prominences/V' );
            ylabel ( 'Counts' ) ;
            set ( gca , 'FontSize' , FontSize ) ;
            set ( gcf , 'Color' , 'w' ) ;
            saveas ( gcf , char ( strcat ( 'PromsHist,' , { ' ' } ,...
                    FileGroup { i } ( 1 : end - 4 ) ,...
                    num2str ( VIterationRange ( j ) ) , ' V prominence.tif' ) ) ) ;
            end
        end
        %% Spike heights.
        if numel ( VoltagePeaks ) < 2
            OutputData ( j , 15 , i ) = 0 ;
            OutputData ( j , 16 , i ) = 0 ;
            OutputData ( j , 17 , i ) = 0 ;
            OutputData ( j , 18 , i ) = 0 ;
            OutputData ( j , 19 , i ) = 0 ;
        else
            HeightsBins = round ( max ( VoltagePeaks ) /...
                ( min ( VoltagePeaks ) ) / 0.1) ;
            if HeightsBins > 512
                HeightsBins = 512 ;
            end
            % Normal distribution might be a more appropriate fit for the
            % spike heights. Also, double peaks appear if there has been a
            % jump in the feedback during the measurement.
            HeightsHist = histfit ( VoltagePeaks , HeightsBins , 'Normal' ) ;
            [ HMax , HIndex ] = max ( HeightsHist ( 1 ) . YData ) ;
            OutputData ( j , 15 , i ) = HeightsHist ( 1 ) . XData ( HIndex ) ;
            OutputData ( j , 16 , i ) = HMax ;
            HeightsDist =  fitdist ( VoltagePeaks , DistType )  ;
            OutputData ( j , 17 , i ) = HeightsDist . mu ;
            OutputData ( j , 18 , i ) = HeightsDist . sigma ;
            OutputData ( j , 19 , i ) = mean ( VoltagePeaks ) ; % Mean peak height.
            if SaveHists == 1
                title ( strcat ( FileGroup { i } ( 1 : end - 4 ) , { ' at ' } ,...
                    num2str ( VIterationRange ( j ) ) , ' V_p' ) ) ;
                xlabel ( 'Spike heights/V' );
                ylabel ( 'Counts' ) ;
                set ( gca , 'FontSize' , FontSize ) ;
                set ( gcf , 'Color' , 'w' ) ;    
                saveas ( gcf , char ( strcat ( 'HeightsHist,' , { ' ' } ,...
                        FileGroup { i } ( 1 : end - 4 ) ,...
                        num2str ( VIterationRange ( j ) ) , ' V prominence.tif' ) ) ) ;
            end
        end
        % Get peak locations as binary values.
        [ VoltagePeaks , BinaryLocs ] = findpeaks ( DataArray { i } { FindIn } ) ;
        BinaryPeaks = zeros ( size ( DataArray { i } { FindIn } , 1 ) , 1 ) ;
        BinaryPeaks ( BinaryLocs ) = 1 ;
        FFTOutput { i , j } = BinaryPeaks ;
    end
end
close ; % Closes histogram figure.
if SaveHists == 1 % Go back to main folder if histograms have been saved.
    cd ../ ;
end

%% Plotting options and fit data.
% Hyperbolic tangent, iFit = 1 (Mehonic, Front. Neurosci., vol. 10, no.
% 57, pp. 1?10, 2016.)
tanhFit = 'a*(1+(tanh((b*(x-c)))))' ;
tanhTitle =  ' tanh fit = p_1 [ 1 + tanh ( p_2(I_{set} - I_0 ))]'  ;
% Sigmoid, iFit = 2, probabilistic early deep neural network
% "squishification" (3Blue1Brown, YouTube).
SigmoidFit =  '(a/(1+exp(b*-x)))+c' ;
SigmoidTitle =  ' sigmoid fit = [p_1 / (1 + e^{-p_2.I_{set}})] + I_0'  ;
% Softplus ReLU, iFit = 3 (rectified linear unit), more current
% "squishifier" (3Blue1Brown, Youtube).
SoftplusReLUFit = '(a*log(1+exp(b*x)))+c' ;
SoftplusReLUETitle =  ' softplus ReLU fit = p_1ln(1+e^{p_2.I_{set}}) + I_0' ;
% Choose fit type for spikes per second data here.
iFit = 1 ; % Fit type index.
if iFit == 1
    RisingThresholdFit = fittype ( tanhFit ,  'Independent' , 'x' ,...
        'Dependent' , 'y' ) ;
    FitTitle = tanhTitle ;
elseif iFit == 2
    RisingThresholdFit = fittype ( SigmoidFit ,  'Independent' , 'x' ,...
        'Dependent' , 'y' ) ;
    FitTitle = SigmoidTitle ;
elseif iFit == 3
    RisingThresholdFit = fittype ( SoftplusReLUFit,  'Independent' , 'x' ,...
        'Dependent' , 'y' ) ;
    FitTitle = SoftplusReLUTitle ;
end
FallingThresholdFit = fittype ( 'a*(1-(tanh((b*x)-c)))' , 'independent' ,...
    'x' , 'dependent' , 'y' ) ;
iColours = 3 ; % Choose colour set to plot with.
if iColours == 1 % Creates rainbow colour set.
    ChosenColours = jet (NumProms ) ;
elseif iColours == 2 % Creates the most distinct colour set.
    ChosenColours = distinguishable_colors ( NumProms ) ;
elseif iColours == 3 % Creates another colour set.
    ChosenColours = lines ( NumProms ) ;
end
LineWidth = 1.2 ; % Choose line width for all plots.
MarkerSize = 5 ; % Choose marker size for all plots.
FontSize = 14 ; % Choose font size for all plots.
FigLineWidth = 1.2 ; % Choose figure box line width for all plots.
FontName = 'Helvetica' ; % Choose font for all plots.

%% Number of spikes per second (from entire range).
% i.e. taking the total number of detected spikes and dividing by the total
% sampling time.
figure ;
SpS = zeros ( NumProms , 4 ) ;
SpSFitOptions = fitoptions ( 'Method' , 'NonlinearLeastSquares' ) ;
SpSFitOptions . Display = 'Off' ;
SpSFitOptions . Robust = 'LAR' ;
SpSFitOptions . StartPoint = [ 15 2 0.1 ] ;
SpSFitOptions . Upper = [ 100 100 0.5 ] ;
SpSFitOptions . Lower = [ -100 -100 -0.5 ] ;
for i = 1 : NumProms
    SpikesPerSecond = reshape ( OutputData ( i , 1 , : ) , NFC , 1 ) ;
    semilogx ( SetPoints , SpikesPerSecond , 'o' , 'Color' ,...
        ChosenColours ( i , : ) , 'MarkerSize' , MarkerSize , 'LineWidth' , LineWidth ) ;
    xlabel ( '{\itI_{set}} [nA]' ) ;
    ylabel ( '\itS_{pS}' ) ;
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
        SpSModel = SpS ( i , 1 ) .* ( 1 + tanh ( ( SpS ( i , 2 ) .* ( VFittingRange -...
            SpS ( i , 3 ) ) ) ) ) ;
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
    semilogx ( VFittingRange , SpSModel  , 'Color' , ChosenColours ( i , : ),...
        'LineWidth' , LineWidth ) ;
    %% Also plot change in fitting parameters at end of loop.
    if i == NumProms
        figure ;
        plot ( VIterationRange , SpS ( : , 1 ) / SpS ( 1 , 1 )  ,...
            VIterationRange , SpS ( : , 2 ) / SpS ( 1 , 2 ) ,...
            VIterationRange , SpS ( : , 3 ) / SpS ( 1 , 3 ) ,...
            VIterationRange , SpS ( : , 4 ) , 'LineWidth' , LineWidth ) ;
        xlabel ( 'Prominence threshold/V' ) ;
        ylabel ( 'Normalized parameter value' ) ;
        legend ( 'p_1' , 'p_2' , 'I_0' , 'Goodness of fit' , 'Location' , 'NorthWest' ) ;
        legend ( 'BoxOff' ) ;
        title ( strcat ( 'SpS' , FitTitle ) ) ;
    end
end

%% Additional plotting for SpS fitting parameters.
iSpS = 0 ;
if iSpS == 1
    figure ;
    semilogy ( SpS ( : , 3 ) , SpS ( : , 1 ) , 'LineWidth' , 3 ) ;
    title ('{\it S_{pS} = S_{pSmax} [ 1 + tanh ( p(I_{set} - I_0 ))]}/2');
    ylabel ( '{\it S_{pSmax}}');
    xlabel ( '{\it I_0} [nA]');

    figure ;
    plot ( SpS ( : , 3 ) , SpS ( : , 2 ) , 'LineWidth' , 3 ) ;
    title ('{\it S_{pS} = S_{pSmax} [ 1 + tanh ( p(I_{set} - I_0 ))]}/2');
    ylabel ( '{\it p}');
    xlabel ( '{\it I_0} [nA]');

    figure ;
    plot ( VIterationRange , SpS ( : , 2 ) , 'LineWidth' , LineWidth ) ;
    title ('{\it S_{pS} = S_{pSmax} [ 1 + tanh ( p(I_{set} - I_0 ))]}/2');
    ylabel ( '{\it p}');
    xlabel ( '{\it V_p} [V]');

    figure ;
    semilogy ( VIterationRange , SpS ( : , 1 ) , 'LineWidth' , 3 ) ;
    title ('{\it S_{pS} = S_{pSmax} [ 1 + tanh ( p(I_{set} - I_0 ))]}/2');
    ylabel ( '{\it S_{pSmax}}');
    xlabel ( '{\it V_p} [V]');

    figure ;
    plot ( VIterationRange , SpS ( : , 3 ) , 'LineWidth', LineWidth );
    title ('{\it S_{pS} = S_{pSmax} [ 1 + tanh ( p(I_{set} - I_0 ))]}/2');
    ylabel ( '{\it I_0} [nA]');
    xlabel ( '{\it V_p} [V]');
end

%% Modal time between spikes.
% i.e. how close together in time can we expect spikes to be.
figure ;
TbS = zeros ( NumProms , 4 ) ;
TbSFitOptions = fitoptions ( 'Method' , 'NonlinearLeastSquares' ) ;
TbSFitOptions . Display = 'Off' ;
TbSFitOptions . Robust = 'LAR' ;
for i = 1 : NumProms
    TimeBetweenSpikes = reshape ( OutputData ( i , 2 , : ) , NFC , 1 ) ;
    loglog ( SetPoints , TimeBetweenSpikes , 'o' , 'Color' ,...
        ChosenColours ( i , : ) , 'MarkerSize' , MarkerSize , 'LineWidth' , LineWidth ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal time between spikes/s' ) ;
    hold on
    [ TbSFit , TbSGoF ] = fit ( SetPoints , TimeBetweenSpikes , 'Power2'  ,...
         TbSFitOptions ) ;
    TbS ( i , 1 ) = TbSFit . a ;
    TbS ( i , 2 ) = TbSFit . b ;
    TbS ( i , 3 ) = TbSFit . c ;
    TbS ( i , 4 ) = TbSGoF . rsquare ;
    loglog ( VFittingRange , ( TbS ( i , 1 ) .* ( VFittingRange .^ TbS ( i , 2 ) ) ) +...
        TbS ( i , 3 ) , 'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
    % Also plot change in fitting parameters at end of loop.
    if i == NumProms
        figure ;
        plot ( VIterationRange , TbS ( : , 1 ) , VIterationRange , TbS ( : , 2 ) ,...
            VIterationRange , TbS ( : , 3 ) , 'LineWidth' , LineWidth ) ;
        xlabel ( 'Prominence threshold/V' ) ;
        ylabel ( 'Parameter value' ) ;
        title ( 'TbS threshold fitting' ) ;
        legend ( 'a' , 'b' , 'c' , 'Location' , 'NorthWest' ) ;
    end
end

%% Effective spike frequency.
% i.e. at a given time, what frequency is the spiking pattern equivalent to.
figure ;
ESF = zeros ( NumProms , 4 ) ;
ESFFitOptions = fitoptions ( 'Method' , 'NonlinearLeastSquares' ) ;
ESFFitOptions . Display = 'Off' ;
ESFFitOptions . Robust = 'LAR' ;
ESFFitOptions . StartPoint = [ 0 0 0 ] ;
for i = 1 : NumProms
    EffectiveSpikeFreq = 1 ./ reshape ( OutputData ( i , 2 , : ) , NFC , 1 ) ;
    semilogx ( SetPoints , EffectiveSpikeFreq , 'o' , 'Color' ,...
        ChosenColours ( i , : ) , 'MarkerSize' , MarkerSize , 'LineWidth' , LineWidth ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Effective spike frequency/Hz' ) ;
    hold on
    [ ESFFit , ESFGoF ] = fit ( SetPoints , EffectiveSpikeFreq , tanhFit  ,...
         ESFFitOptions ) ;
    ESF ( i , 1 ) = ESFFit . a ;
    ESF ( i , 2 ) = ESFFit . b ;
    ESF ( i , 3 ) = ESFFit . c ;
    ESF ( i , 4 ) = ESFGoF . rsquare ;
    % tanh fitting seems the most appropriate, here.
    ESFModel = ESF ( i , 1 ) .* ( 1 + tanh ( ( ESF ( i , 2 ) .* VFittingRange ) -...
     	ESF ( i , 3 ) ) ) ;
%     ESFModel = ( ESF ( i , 1 ) ./ ( 1 + exp ( ESF ( i , 2 ) .* ( -1 ) .*...
%     	VFittingRange ) ) ) + ESF ( i , 3 ) ;
%     ESFModel = ( ESF ( i , 1 ) .* log ( 1 + exp ( ( ESF ( i , 2 ) .*...
%         VFittingRange ) ) ) ) + ESF ( i , 3 ) ;
    semilogx ( VFittingRange , ESFModel , 'Color' , ChosenColours ( i , : ) ,...
        'LineWidth' , LineWidth ) ;
    % Also plot change in fitting parameters at end of loop.
    if i == NumProms
        figure ;
        plot ( VIterationRange , ESF ( : , 1 ) , VIterationRange , ESF ( : , 2 ) ,...
            VIterationRange , ESF ( : , 3 ) , 'LineWidth' , LineWidth ) ;
        xlabel ( 'Prominence threshold/V' ) ;
        ylabel ( 'Parameter value' ) ;
        title ( 'ESF threshold fitting' ) ;
        legend ( 'a' , 'b' , 'c' , 'Location' , 'NorthWest' ) ;
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
    loglog ( SetPoints , ModalSpikeWidths , 'o' , 'Color' ,...
        ChosenColours ( i , : ) , 'MarkerSize' , MarkerSize , 'LineWidth' , LineWidth ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal spike width/s' ) ;
    hold on
	[ SWFit , SWFitGoF ] = fit ( SetPoints , ModalSpikeWidths , '(a*exp(-b*x))+c' ,...
         SWFitOptions ) ; % Need to choose the best fitting option here.
	 % (a*exp(-b*x))+c, a single exponential, works well.
     % The softplus ReLU also works well.
     % A sigmoid, also works well.
    SW ( i , 1 ) = SWFit . a ;
    SW ( i , 2 ) = SWFit . b ;
    SW ( i , 3 ) = SWFit . c ;
    SW ( i , 4 ) = SWFitGoF . rsquare ;
    loglog ( VFittingRange , ( SW ( i , 1 ) .* exp ( ( -1 ) .* SW ( i , 2 ) .*...
        VFittingRange ) ) + SW ( i , 3 ) , 'Color' , ChosenColours ( i , : ),...
        'LineWidth' , LineWidth ) ;
end

%% Modal energy per spike.
% Energy per spike = modal spike FWHM x current setpoint x mean spike height.
% Using prominence is an idealised case; this assumes that spikes aren't
% typically much greater than the desired promience, in which case the
% energy per SynOps would increase.
figure ;
EpS = zeros ( NumProms , 4 ) ;
EpSFitOptions = fitoptions ( 'Method' , 'LinearLeastSquares' ) ;
EpSFitOptions . Robust = 'LAR' ;
EpSFitOptions . Lower = [ -0.05 -0.05 ] ;
EpSFitOptions . Upper = [ 0.05 0.05 ] ;
for i = 1 : NumProms
    EnergyPerSpike = reshape ( OutputData ( i , 6 , : ) , NFC , 1 ) .*...
        SetPoints .* reshape ( OutputData ( i , 19 , : ) , NFC , 1 ) ;
    loglog ( SetPoints , EnergyPerSpike , 'o' , 'Color' , ChosenColours ( i , : ),...
        'MarkerSize' , MarkerSize , 'LineWidth' , LineWidth ) ;
    xlabel ( '{\it I_{set}} [nA]' ) ;
    ylabel ( '{\itE_{pS}} [nJ]' ) ;
    hold on
	[ EpSFit , EpSFitGoF ] = fit ( SetPoints , EnergyPerSpike , 'Poly1' ,...
        EpSFitOptions ) ; % Linear fitting works best here.
    EpS ( i , 1 ) = EpSFit . p1 ;
    EpS ( i , 2 ) = EpSFit . p2 ;
    EpS ( i , 4 ) = EpSFitGoF . rsquare ;
    loglog ( VFittingRange , ( EpS ( i , 1 ) .* VFittingRange ) + EpS ( i , 2 ) ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
end

%% Spike separation distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 4 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike separation µ/-' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 5 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
        ylabel ( 'Spike separation {\sigma}/--' ) ;
end

%% Spike width distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 9 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike width µ/-' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 10 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
        ylabel ( 'Spike width {\sigma}/--' ) ;
end

%% Spike prominence distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 13 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike prominence µ/-' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 14 , : ) , NFC , 1 ) ) , '--' ,....
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
        ylabel ( 'Spike prominence {\sigma}/--' ) ;
end

%% Modal spike heights.
figure ;
for i = 1 : NumProms
    ModalSpikeHeights = reshape ( OutputData ( i , 15 , : ) , NFC , 1 ) ;
    semilogx ( SetPoints , ModalSpikeHeights  , 'Color' ,...
        ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
	xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Modal spike height/V' ) ;
    hold on
end

%% Mean spike heights.
figure ;
for i = 1 : NumProms
    MeanSpikeHeights = reshape ( OutputData ( i , 19 , : ) , NFC , 1 ) ;
    semilogx ( SetPoints , MeanSpikeHeights  , 'Color' ,...
        ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
	xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Mean spike height/V' ) ;
    hold on
end

%% Spike height distribution.
figure ;
for i = 1 : NumProms
    yyaxis left
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 17 , : ) , NFC , 1 ) ) , '-' ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
    xlabel ( 'Current bias/nA' ) ;
    ylabel ( 'Spike height µ/-' ) ;
    hold on
    yyaxis right
    semilogx ( SetPoints , ( reshape ( OutputData ( i , 18 , : ) , NFC , 1 ) ) , '--' ,...
        'Color' , ChosenColours ( i , : ) , 'LineWidth' , LineWidth ) ;
        ylabel ( 'Spike height {\sigma}/--' ) ;
end

%% Autocorrelation.
iAuto = 1 ;
if iAuto == 1
    for i = 1 : size ( DataArray , 1 )
        figure ;
        AutoSpikes = autocorr ( DataArray { i } { 6 } , 'NumLags' ,...
            size ( DataArray { i } { 3 } , 1 ) - 1 ) ;
        AutoTime = ( 0 : 1 : size ( AutoSpikes ) - 1 ) .* ( 1 / SampleFrequency ( i , 1 ) ) ;
        plot ( AutoTime , AutoSpikes ) ;
        xlabel ( 'Lag time [s]' ) ;
        ylabel ( 'ACF [a.u.]' ) ;
    end
end

%% FFT from binary spike data.
iFFT = 0 ;
if iFFT == 1
    figure ;
    for i = 1 : NFC
        for j = 1 : NumProms
        SpikeFFT = fft ( FFTOutput { i , j } ) ; % FFT of binary spikes.
        FFTFrequency = ( 1 ./ SampleFrequency ( i ) ) *...
            ( 0 : ( length ( SpikeFFT ) / 2 ) ) / length ( SpikeFFT ) ;
          SpikeSpectrum2 = abs ( SpikeFFT / length ( SpikeFFT ) ) ; % Two-sided spectrum.
          SpikeSpectrum1 = SpikeSpectrum2 ( 1 : ( length ( SpikeFFT ) / 2 + 1 ) ) ; % Single-sided spectrum.
        SpikeSpectrum1 ( 2 : end - 1 ) = 2 * SpikeSpectrum1 ( 2 : end - 1 ) ;
        semilogx ( FFTFrequency , SpikeSpectrum1 , 'LineWidth' , LineWidth ) ;
        end
    end
end
%% Format all figures.
iF = 1 ;
if iF == 1
    NumPlots = findobj ( 'Type' , 'Figure' ) ;
    for i = 1 : length ( NumPlots )
        figure ( NumPlots ( i ) . Number ) ;
        set ( gca , 'FontSize' , FontSize , 'FontName' , FontName ,...
            'LineWidth' , FigLineWidth , 'Tickdir' , 'Out' , 'Box' , 'Off' ) ;
        set ( gcf , 'Color' , 'w' ) ;
    end
end
