function [ ThermalArray , ThermLocation ] = CantileverThermalProcessorFunction ( DeflectionSensitivity , SensitivityError )

%% Imports and processes sets of .txt peak data files cotaining two columns.
% The first column should be a scale and the second an intensity. A single peak per
% file is assumed, though this may be a convolution of smaller peaks. The data is
% fitted using peak.m. The calibration section allows a set of peaks, ideally
% over a range of change of a single sampling parameter, to be fitted in
% order to assess the effect of that parameter on the peak height and area.
% The subsequent fitting section should be used to fit a set up single
% peaks, ideally with the same sampling parameters. These will be summed
% and fitted in order to find the area and height. These values will then be
% compared to the calibration value of the chosen sampling parameter (whether
% using a set of calibration curves or writing in the fitting parameters) in
% order to scale the data to the real value of this parameter, assuming it is not
% well-determined at the time of sampling. The the mean, standard deviation and
% standard error on the peak positions are not calculated as, assuming multiple
% peak files are input, these data are summed in order to produce a single peak data
% set with improved signal to noise ratio. Designed for use in extracting precise
% peak areas from AFM thermal tune data from for example, QNM mode.

% Takes deflection sensitivity and error as inputs from characerisation suite script.

% Gets data from selected files.
addpath ( cd ) ;

% Calibration sensitivity values.
SensInputArray = 0 : 5 : 125 ;
SensInputArray ( 1 ) = 1 ;
SensInputArray = horzcat ( 0 , SensInputArray ) ;
SensInputArrayCrop = SensInputArray' ;
HighResInputArray = ( SensInputArray ( 1 ) : 0.1 : SensInputArray ( end ) )' ;

Background = 1 ; % 1 = linear background, 2 = quadratic background.

% Inital fitting parameters.
FileName = strcat ( 'parameters' , '.txt' ) ;
if exist ( FileName , 'file' ) > 0 % Checks to see if peakfit.m parameters file exists.

    ParametersFile = fopen ( FileName ) ; % If present, parameters are opened.
    FileParameters = textscan ( ParametersFile , '%f' ) ;
    Parameters = FileParameters { 1 } ;

    NumPeaks = Parameters ( 1 ) ; % Number of peaks to fit to data.
    PeakShape = Parameters ( 2 ) ; % 1 = Gaussian, 2 = Lorentzian, 5 = exponentially-broadened Gaussian, 13 = Gaussian/Lorentzian blend,
                                                           % 18 = exponentially-broadened Lorentzian, 20 = Voigt.
    NumTrials = Parameters ( 3 ) ; % Number of trial fits.

    % Exponential/Voigt broadening factor. For Voigt profiel, cannot be negative.
    if PeakShape == 20

        Extra = abs ( Parameters ( 4 ) ) ;

    else

        Extra = Parameters ( 4 ) ;

    end

    Bipolar = Parameters ( 5 ) ;
    Variance = Parameters ( 6 ) ;
    CropSize = Parameters ( 7 ) ;

else % Otherwise, use basic/standard set of parameters.

    NumPeaks = 1 ;  % Number of peaks to fit to data.
    PeakShape = 2 ; % 1 = Gaussian, 2 = Lorentzian, 5 = exponentially-broadened Gaussian, 13 = Gaussian/Lorentzian blend,
                                   % 18 = exponentially-broadened Lorentzian, 20 = Voigt.
    NumTrials = 20 ; % Number of trial fits.
    Extra = 1 ; % Exponential/Voigt broadening factor.
    Bipolar = 0 ; % 0 = all peaks positive, 1 = can use negative peaks.
    Variance = 1 ; % Variation between fit attempts when NumPeaks > 1.
    CropSize = 0 ; % 0 = no data cropping, > 0 = crop data to this many points either side of maximum y value.

end

iSet = 1 ; % First run, set at 0 to input calibration files or 1 to use in-built calibration and proceed to data fitting.
if iSet == 0

    [ NFC , DataArray ] = CantileverThermalFileSelector ( iSet ) ;

    %% For calibrating data scale with input deflection sensitivity.

    IntensityArray = DataArray ( : , 2 , : ) ; % Array of just intensities.
    IntensityArray = reshape ( IntensityArray , length ( IntensityArray ) , NFC ) ;
    MaxArray = max ( IntensityArray ) ;

    SensFit = polyfit ( SensInputArray , MaxArray , 2 ) ; % Quadratic fit to the sensitivity vs. max signal.
    FitArray =  ( SensFit ( 1 ) .* ( SensInputArray .^ 2 ) + ( SensFit ( 2 ) .* SensInputArray )...
        + SensFit ( 3 ) ) ; % Fitting data.

    figure ( 2 ) ;
    plot ( SensInputArray , MaxArray , SensInputArray , FitArray ) ;
    xlabel ( 'Deflection sensitivity input (mV)' ) ;
    ylabel ( 'Maximum singnal intensity (fm^2/Hz)' ) ;

    [ fMax , fIndex ] = max ( DataArray ( : , 2 , 2 ) ) ; % Find the peak and sample around it.
    fCropSize = 750 ;
    Xnmin = fIndex - fCropSize ;
    Xnmax = fIndex + fCropSize ;
    QuantArray = zeros ( NFC , 3 ) ; % Fitted peak size/height quantification array.

    for i = 2 : NFC

        Xn = DataArray ( Xnmin : Xnmax , 1 , i ) ;
        Yn = DataArray ( Xnmin : Xnmax , 2 , i ) ;
        figure ( 4 ) ;
        [ MultiResults , MultiGoodness , MultiBase , MultiCoeff , MultiResidual , Xj , Yj ] = peakfit ( [ Xn Yn ]...
            , 0 , 0 , NumPeaks , PeakShape , Extra , NumTrials , 0 , Background , 0 , 1 , Bipolar , 0 , Variance ) ;
    %     close ( gcf ) ;
        MultiB = ( MultiBase ( 1 ) * Xj' ) + MultiBase ( 2 ) ;                                                                                                                                      
        YjC = Yj' +  MultiB ;
        Yj = Yj - min ( Yj ) ;
        figure ( 3 ) ;
        hold on
        plot ( Xj , Yj ) ;
        QuantArray ( i , 1 ) = max ( Yj ) ; % Peak height.
        QuantArray ( i , 2 ) = sum ( Yj ) ; % Area under peak.
        QuantArray ( i , 3 ) = max ( YjC ) ; % Peak height without baseline subtraction.

    end

    xlabel ( 'Frequency (Hz)' ) ;
    ylabel ( 'Intensity (fm^2/Hz)' ) ;

    % SensInputArrayCrop ( 2 ) = [] ; % Remove value of 1 that skews fit.
    % QuantArray ( 2 , : ) = [] ;
    SensHeightFit = polyfit ( SensInputArrayCrop , QuantArray ( : , 1 ) , 2 ) ; % Quadratic fit to the fitted peak heights.
    SensHeightFitArray =  ( SensHeightFit ( 1 ) .* ( HighResInputArray .^ 2 ) + ( SensHeightFit ( 2 ) .* HighResInputArray )...
        + SensHeightFit ( 3 ) ) ; % Fitting height data.
    figure ( 5 ) ;
    plot ( SensInputArrayCrop , QuantArray ( : , 1 ) , HighResInputArray , SensHeightFitArray ) ;
    
else
    
    SensHeightFit = [ 4.84134733788011 , 4.08962592567363 , -18.0846280945172 ] ; % Quadratic calibration parameters for fitted peak heights.
    SensHeightFitArray =  ( SensHeightFit ( 1 ) .* ( HighResInputArray .^ 2 ) + ( SensHeightFit ( 2 ) .* HighResInputArray )...
        + SensHeightFit ( 3 ) ) ; % Fitting height data.
    
end

iSet = 1 ; % Second run to input data files.
[ NFC , DataArray ] = CantileverThermalFileSelector ( iSet ) ;


iAccept = 0 ; % Whether or not to run fitting.

while iAccept == 0

    figure ( 1 ) ;

    % Defines data x and y values.
    X = DataArray ( : , 1 ) ;
    Y = DataArray ( : , 2 ) ./ NFC ; % Y data normalised by number of files for single peak height, 

    if CropSize > 0

        % Defines range of values to fit, if not using all peak data.
        [ MaxVal , MaxIndex ] = max ( Y ) ;
        X = X ( MaxIndex - CropSize : MaxIndex + CropSize ) ;
        Y = Y ( MaxIndex - CropSize : MaxIndex + CropSize ) ;

    end

    % Use peakfit.m to fit all or part of peak data.
    [ FitResults , GoodnessOfFit , Baseline , Coefficients , Residual , Xi , Yi ] = peakfit ( [ X Y ] , 0 , 0 , NumPeaks , PeakShape , Extra , ...
        NumTrials , 0 , Background , 0 , 1 , Bipolar , 0 , Variance ) ;
%     LinearB = ( Baseline ( 1 ) * Xi' ) + Baseline ( 2 ) ;
    if NumPeaks > 1
    
        Yn = sum ( Yi ) ; % Collapses components into single peak.
       
    else
        
        Yn = Yi ; % Preserve single component.
        
    end
    
    hold on
    
    % Extracts peak centre value.
    Yn = Yn - min ( Yn ) ;
    [ YMax , MIndex ] = max ( Yn ) ;
    fRMultiFit = Xi ( MIndex ) ;
    
    subplot ( 2 , 1 , 1 ) ;
    legend ( strcat ( 'Position:' ,  num2str ( fRMultiFit , 20 ) ) ) ;
    
    Accept = MFquestdlg ( [ 0.7 0.73 ] , 'Accept fit, change parameters or repeat?' , 'Check fit suitability' , 'Accept' , ...
        'Change parameters' , 'Repeat' , 'Repeat' ) ;

    switch Accept % Let the operator decide whether to use the fitted data,
                              % repeat the fit or change parameters. Script can be terminated by selecting 'Change arameters'
                              % and closing the dialog.

        case 'Accept'  % Keep data and move on to next peak.

            % Retains current fit parameters and saves a copy of the current figure..
            save ( strcat ( 'parameters' , '.txt' ) , 'NumPeaks', 'PeakShape' , 'NumTrials' , 'Extra' , 'Bipolar' , 'Variance' , 'CropSize' , '-ascii' ) ;
            savefig ( 1 , strcat ( 'thermalplot' , '.fig' ) ) ; % Saves current fit, overwriting previous.

            iAccept = 1 ; % Check for parameters file on next peak.

        case 'Change parameters' % Allow user to change parameters and repeat fit for same peak.
                                                       % User can terminate by closing this dialog.

            UserPrompt = { 'Numer of peaks' , 'Peak shape' , 'Numer of trials' , 'Broadening' , 'Unipolar (0) or bipolar (1)?' , 'Variance' , ...
                'No cropping (0) or define crop size (>0)?' } ;
            PromptTitle = 'Adjust fit paremeters' ;
            PromptLines = 1 ;
            DefaultAnswers = { ( num2str ( NumPeaks ) ) ( num2str ( PeakShape ) ) ( num2str ( NumTrials ) ) ( num2str ( Extra ) ) ...
                ( num2str ( Bipolar ) ) ( num2str ( Variance ) ) ( num2str ( CropSize ) ) } ;
            AlterProperties = inputdlg ( UserPrompt , PromptTitle , PromptLines , DefaultAnswers ) ;

            NumPeaks = str2double ( AlterProperties { 1 } ) ; % Number of peaks to fit to data.
            PeakShape = str2double ( AlterProperties { 2 } ) ; % 1 = Gaussian, 2 = Lorentzian, 5 = exponentially-broadened Gaussian,
                                                                                                   % 13 = Gaussian/Lorentzian blend, 18 = exponentially-broadened Lorentzian,
                                                                                                   % 20 = Voigt.
            NumTrials = str2double ( AlterProperties { 3 } ) ; % Number of trial fits.

            % Exponential/Voigt broadening factor. For Voigt profile, cannot be negative.
            if PeakShape == 20

                Extra = abs ( str2double ( AlterProperties { 4 } ) ) ;

            else

                Extra = str2double ( AlterProperties { 4 } ) ;

            end

            Bipolar = str2double ( AlterProperties { 5 } ) ;
            Variance = str2double ( AlterProperties { 6 } ) ;
            CropSize = str2double ( AlterProperties { 7 } ) ;
            
    end

end

ThermalArray = cell ( 2 , 2 ) ; % Output array for values of interest.

DeflPrompt = { 'Input deflection sensitivity of measurement in nm/V' } ;
PromptTitle = 'Input deflection sensitivity in nm/V' ;
PromptLines = 1 ;
DeflInput = inputdlg ( DeflPrompt , PromptTitle , PromptLines ) ;
DeflMeasure = str2double ( DeflInput ( 1 ) ) ; % Deflection sensitivity value of measurements.
[ MeasureClosest , MeasureIndex ] = min ( abs ( ( HighResInputArray - DeflMeasure ) ) ) ;
[ RealClosest , RealIndex ] = min ( abs ( ( HighResInputArray - 62.6119478149421 ) ) ) ;
MeasureRatio = YMax / SensHeightFitArray ( MeasureIndex ) ;
RealSignal = MeasureRatio * SensHeightFitArray ( RealIndex ) ;

% Calculates peak areas and spring constants assuming temperature of 21C.
XScale = ( Xi ( 2 ) - Xi ( 1 ) ) * 1E3 ; % Get scaling factor for summing along frequency axis.
YNorm = ( Yn ./ max ( Yn ) ) .* RealSignal ; % Scales measured signal according to real sensitivity.
ThermalArea = ( sum ( YNorm ) * 1E-30 ) * XScale ; % Sum of fitted peak data.
ScaleConst = 0.971 / ( 1.08 ^ 2 ) ;
Temp = 21 ;
kBoltz = 1.38064852E-23 ;
k = ( ScaleConst * kBoltz * ( 237.15 + Temp ) ) / ThermalArea ; % Spring constant for full peak.

ThermalArray { 1 } = 'Thermal k (N/m)' ;
ThermalArray { 2 } = 'Thermal k % error' ;
ThermalArray { 3 } = k ;
ThermalArray { 4 } = ( 0.27512599546642 / 62.6119478149421 ) * 100 ;

ThermLocation = cd ;

end
