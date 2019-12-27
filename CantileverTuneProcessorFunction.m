function [ TuneArray , TuneLocation ] = CantileverTuneProcessorFunction ( )

%% Imports and processes sets of .txt peak data files cotaining two columns.
% The first column should be a scale and the second an intensity. A single peak per
% file is assumed, though this may be a convolution of smaller peaks. The data is
% fitted using peak.m and findpeaksL.m, with the mean, standard deviation and
% standard error on the peak positions calculated. These values are also computed
% for the fitted curve in the region around the peak maximum, in order to assess
% the goodness of the fit. Designed for use in extracting precise values for the
% resonant frequency of an AFM cantilever when tuned in tapping mode.

 % Gets data from selected files.
addpath ( cd ) ;
[ NFC , DataArray ] = CantileverTuneFileSelector ( ) ;

figure ( 1 ) ; % Makes figure 1 current, so new fits will always be placed here.

fRArray = zeros ( NFC , 2 ) ; % Recipient array for resonant frequency values.

iNFC = 1 ; % File importation index.
iMaximise = 1 ; % Figure window maximisation index.
iFileCheck = 1 ; % File checking index.
    
Background = 1 ; % 1 = linear background, 2 = quadratic background.   
SamplePoints = 25 ; % Number of data points to use in smoothing for LFit.
SampleSize = 152 ; % Number of data points to use in fitting for LFit.

while iNFC <= NFC % Fit data.
    
    if iFileCheck == 1 % Only check on presence of parameters file once.
        
        % Inital fitting parameters.
        FileName = strcat ( 'parameters' , num2str ( iNFC ) , '.txt' ) ;
        LastFileName = strcat ( 'parameters' , num2str ( ( iNFC - 1 ) ) , '.txt' ) ;
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
            
        elseif exist ( LastFileName , 'file' ) > 0 % Checks to see if previous parameters file exists.
            
            ParametersFile = fopen ( LastFileName ) ; % If present, uses previous parameters.
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

            NumPeaks = 7 ;  % Number of peaks to fit to data.
            PeakShape = 20 ; % 1 = Gaussian, 2 = Lorentzian, 5 = exponentially-broadened Gaussian, 13 = Gaussian/Lorentzian blend,
                                           % 18 = exponentially-broadened Lorentzian, 20 = Voigt.
            NumTrials = 5 ; % Number of trial fits.
            Extra = 1.02 ; % Exponential/Voigt broadening factor.
            Bipolar = 0 ; % 0 = all peaks positive, 1 = can use negative peaks.
            Variance = 1 ; % Variation between fit attempts when NumPeaks > 1.
            CropSize = 0 ; % 0 = no data cropping, > 0 = crop data to this many points either side of maximum y value.

        end
    
    end
    
    clf % Clears  figure contents.
    
    % Defines data x and y values.
    X = DataArray ( : , 1 , iNFC ) ;
    Y = DataArray ( : , 2 , iNFC ) ;
    
    if CropSize > 0
        
        % Defines range of values to fit, if not using all peak data.
        [ MaxVal , MaxIndex ] = max ( Y ) ;
        X = X ( MaxIndex - CropSize : MaxIndex + CropSize ) ;
        Y = Y ( MaxIndex - CropSize : MaxIndex + CropSize ) ;
     
    end
    
    % Use peakfit.m to fit all or part of peak data.
%     MinWidth = ( X ( 2 ) - X ( 1 ) ) ; % Minimum peak width, here set to x-axis interval.
    [ FitResults , GoodnessOfFit , Baseline , Coefficients , Residual , Xi , Yi ] = peakfit ( [ X Y ] , 0 , 0 , NumPeaks , PeakShape , Extra , ...
        NumTrials , 0 , Background , 0 , 1 , Bipolar , 0 , Variance ) ;
    LinearB = ( Baseline ( 1 ) * Xi' ) + Baseline ( 2 ) ;                                                                                                                                      
    YFit = ( sum ( Yi , 1 ) )' +  LinearB ;
    XFit = Xi' ;
    
    % Computes statistics of fit, using area around peak maximum, and displays bounds on plot.
    StatWidth = round ( numel ( X ) / 3 ) ;
    StatData = Y - interp1 ( XFit , YFit , X , 'Spline' ) ; % Interpolating the fit to the same length as the data, calculates the residual.
    ResidualMean = mean ( StatData ( StatWidth : end - StatWidth ) ) ;
    ResidualStD = std ( StatData ( StatWidth : end - StatWidth ) ) ;
    ResidualStE = ResidualStD ./ sqrt ( numel ( StatData ) ) ;
    LowerBound = X ( StatWidth ) .* ones ( numel ( X ) , 1 ) ;
    UpperBound = X ( end - StatWidth ) .* ones ( numel ( X ) , 1 ) ;
    subplot ( 2 , 1 , 1 ) ;
    hold on
    plot ( LowerBound , Y , UpperBound , Y ) ;
    
    % Extracts peak centre value.
    [ YMax , MIndex ] = max ( YFit ) ;
    fRMultiFit = XFit ( MIndex ) ;
    
    % Use findpeaksL.m to fit peak maximum region.
    LorentzFinder = findpeaksL ( X , Y , 0 , 0.2 , SamplePoints , SampleSize , 3 ) ;
    fRLFit = LorentzFinder { ( 2 ) } ;
    LX = LorentzFinder { ( 6 ) } ;
%     LY = LorentzFinder { ( 7 ) } ;
    LZ = LorentzFinder { ( 8 ) } ;
    Coeffs = LorentzFinder { ( 9 ) } ;
    
    % Replots findpeaksL.m fit using second order polynomial fit coefficients.
    LFit = ( Coeffs ( 1 ) .* ( LZ .^ 2 ) ) + ( Coeffs ( 2 ) .* LZ ) + Coeffs ( 3 ) ;
    LFitAdj = ( LFit - min ( LFit ) ) ;
    LFitNorm = LFitAdj ./ max ( LFitAdj ) ;
    
    % Plots pair of fits and data, zoomed in to peak maximum.
    subplot ( 2 , 1 , 2 ) ;
    plot ( X , Y , '.' , XFit , YFit , LX , ( ( LFitNorm .* max ( Y ) ) + 0.02 ) , X , ( Y + 0.02 ) , '.' , LowerBound , Y , UpperBound , Y ) ;
    xlim ( [ ( LowerBound ( 1 )  ) ( UpperBound ( 1 ) ) ] ) ;
    ylim ( [ ( max ( Y ) - 0.25 ) ( max ( Y ) + 0.03 ) ] ) ;
    xlabel ( strcat ( 'Number of trials:' , { ' ' } , num2str ( NumTrials ) , ', Fit', { ' ' } , num2str ( iNFC ) , ' of' , { ' ' } , num2str ( NFC ) ) ) ;
    legend ( 'Data' ,  num2str ( fRMultiFit , 20 ) , num2str ( fRLFit , 20 ) ) ;
    Statistics = char ( strcat ( 'Residual mean:' , { ' ' } , num2str ( ResidualMean ) , { ' ' } , ...
        'Residual stdev:' , { ' ' } , num2str ( ResidualStD ) , { ' ' } , ...
        'Residual sterr:' , { ' ' } , num2str ( ResidualStE ) ) ) ;
    delete ( findall ( gcf , 'type' , 'annotation' ) ) ;
    annotation ( 'TextBox' , [ 0.15 0.44 0.15 0 ] , 'String' , 'Current' , 'FitBoxToText' , 'off' ) ;
    annotation ( 'TextBox' , [ 0.15 0.415 0.15 0 ] , 'String' , Statistics , 'FitBoxToText' , 'off' ) ;
    
    if iMaximise == 1 % Executes only on first run.
        
        Maximise = questdlg ( 'Maximise figure window?' , 'Maximise figure window?' , 'Yes' , 'No' , 'Yes' ) ;
    
        switch Maximise % User option to maximise figure window for ease of inspection.

            case 'Yes'

                set ( gcf , 'units' , 'normalized' , 'outerposition' , [ 0 0 1 1 ] ) ;

        end
    
        StatisticsLast = Statistics ; % Retain statistics for next iteration.

        iMaximise = iMaximise + 1 ;
        
    else
        
        % Print last statistics and retain current.
        annotation ( 'TextBox' , [ 0.15 0.34 0.15 0 ] , 'String' , 'Last' , 'FitBoxToText' , 'Off' ) ; 
        annotation ( 'TextBox' , [ 0.15 0.315 0.15 0 ] , 'String' , StatisticsLast , 'FitBoxToText' , 'Off' ) ;
        StatisticsLast = Statistics ;
            
    end
    
    Accept = MFquestdlg ( [ 0.7 0.73 ] , 'Accept fit, change parameters or repeat?' , 'Check fit suitability' , 'Accept' , ...
        'Change parameters' , 'Repeat' , 'Repeat' ) ;
        
    switch Accept % Let the operator decide whether to use the fitted data,
                              % repeat the fit or change parameters. Script can be terminated by selecting 'Change arameters'
                              % and closing the dialog.
    
        case 'Accept'  % Keep data and move on to next peak.

            fRArray ( iNFC , 1 ) = fRMultiFit ;
            fRArray ( iNFC , 2 ) = fRLFit ;
            
            % Retains current fit parameters and saves a copy of the current figure..
            save ( strcat ( 'parameters' , num2str ( iNFC ) , '.txt' ) , 'NumPeaks', 'PeakShape' , 'NumTrials' , 'Extra' , 'Bipolar' , 'Variance' , 'CropSize' , '-ascii' ) ;
            savefig ( 1 , strcat ( 'tuneplot' , num2str ( iNFC ) , '.fig' ) ) ; % Saves current fit, overwriting previous.
            
            iNFC = iNFC + 1 ; % Move to next peak.
            iFileCheck = 1 ; % Check for parameters file on next peak.
            
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

            % Retains current fit parameters.
            iFileCheck = 2 ; % Don't check parameters file on subsequent iteration of parameters.
            
    end
    
end

TuneArray = cell ( 6 , 2 ) ; % Output array for values of interest.

fRMeanMulti = mean ( fRArray ( : , 1 ) ) ; % Mean from peakfit.m. peak centre values.
fRMeanSingle = mean ( fRArray ( : , 2 ) ) ; % Mean from findpeaksL.m peak centre values.
StDevMulti = std ( fRArray ( : , 1 ) ) ; % Standard deviation of peakfit.m centre values.
StDevSingle = std ( fRArray ( : , 2 ) ) ; % Standard deviation of findpeaksL.m centre values.
StErrorMulti = StDevMulti / sqrt ( NFC ) ; % Standard error of peakfit.m centre values.
StErrorSingle = StDevSingle / sqrt ( NFC ) ; % Standard error of findpeaksL.m centre values.

TuneArray { 1 } = 'fRMeanMulti (kHz)' ;
TuneArray { 2 } = 'fRMeanSingle (kHz)' ;
TuneArray { 3 } = 'StDevMulti (kHz)' ;
TuneArray { 4 } = 'StDevSingle (kHz)' ;
TuneArray { 5 } = 'StErrorMulti (kHz)' ;
TuneArray { 6 } = 'StErrorSingle (kHz)' ;
TuneArray { 7 } = fRMeanMulti ;
TuneArray { 8 } = fRMeanSingle ;
TuneArray { 9 } = StDevMulti ;
TuneArray { 10 } = StDevSingle ;
TuneArray { 11 } = StErrorMulti ;
TuneArray { 12 } = StErrorSingle ;

TuneLocation = cd ;

end
