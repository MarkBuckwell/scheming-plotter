function [ EstimateArray ] = CantileverForceEstimatesFunction ( ResonantFrequency , DeflectionSensitivity )

%% Etimates DDESP-type effective probe mass or spring constant range of
% variation with resonant frequency using nominal and max/min quoted values for k, f0 and l.

ResPrompt = { 'Probe type - 1 = DDESP V2, 2 = DDESPFM V2, 3 = custom ranges.' } ;
PromptTitle = 'Input probe type' ;
PromptLines = 1 ;
ResInput = inputdlg ( ResPrompt , PromptTitle , PromptLines ) ;
DDESP = str2double ( ResInput ) ; % Probe type, use 1 for DDESP V2, 2 for DDESPFM V2

if DDESP == 1
    
    % DDESP V2 quoted parameters.
    fRRange = [ 300 , 450 , 600 ] ; % Range of resonant frequency values.
    lRange = [ 115 125 135 ] ; % Range of lengths.
    kRange = [ 30 , 80 , 180 ] ; % Range of spring constants.

elseif DDESP == 2 
    
    % DDESPFM V2 quoted parameters.
    fRRange = [ 80 , 105 , 130 ] ; % Range of resonant frequency values.
    lRange = [ 215 225 235 ] ; % Range of lengths.
    kRange = [ 3 , 6 , 12 ] ; % Range of spring constants.
    
elseif DDESP == 3
    
    % Allow user to input their own cantilever parameter ranges.
    
    UserPrompt = { 'Minimum fR (kHz)' , 'Nominal fR (kHz)' , 'Maximum fR (kHz)' ,...
        'Min length (um)' , 'Nominal length (um)' , 'Maximum length (um)' ,...
        'Minimum spring constant (N/m)' , 'Nominal spring constant (N/m)' ,...
        'Maximum spring constant (N/m)'} ;
            PromptTitle = 'Input cantilever parameter ranges:' ;
            PromptLines = 1 ;
            DefaultAnswers = { '50' , '70' , '90' , '225' , '240' , '255' ,...
                '0.6' , '2' , '3.5' } ; % Quoted values for OSCM-PT-R3
            InputRanges = inputdlg ( UserPrompt , PromptTitle , PromptLines , DefaultAnswers ) ;
            fRRange = zeros ( 1 , 3 ) ;
            lRange = fRRange ;
            kRange = fRRange ;
            
            for i = 1 : 3
                
                fRRange ( i ) = str2double ( InputRanges { i } ) ;
                lRange ( i ) = str2double ( InputRanges { i + 3 } ) ;
                kRange ( i ) = str2double ( InputRanges { i + 6 } ) ;
                
            end
            
else
    
    % Use DDESP V2 quoted parameters as default.
    fRRange = [ 300 , 450 , 600 ] ; % Range of resonant frequency values.
    lRange = [ 115 125 135 ] ; % Range of lengths.
    kRange = [ 30 , 80 , 180 ] ; % Range of spring constants.
    
end

% Value ranges for fitting.
klLowers = kRange .* ( ( lRange ( 2 ) / lRange ( 3 ) ) ^ 3 ) ; % Spring contants for longest cantilever.
klLowers ( 1 ) = kRange ( 1 ) ; % Pins min value to min quoted value.
klUppers = kRange .* ( ( lRange ( 2 ) / lRange ( 1 ) ) ^ 3 ) ; % Spring contants for shortest cantilever.
klUppers ( 3 ) = kRange ( 3 ) ; % Pins max value to max quoted value.
mRange = kRange ./ ( ( 2 * pi .* ( fRRange .* 1E3 ) ) .^ 2 ) ; % Mass values, assuming simple harmonic oscillator model.
fRModel = ( fRRange ( 1 ) : 0.00001 : fRRange ( end ) ) ; % Higher resolution resonant frequency values for modelling.

% Fits quoted parameters with quadratic functions.
RangeFit = polyfit ( fRRange , kRange , 2 ) ; % Fits spring constant to resonant frequency for nominal cantilever length.
kModel = ( RangeFit ( 1 )  .* ( fRModel .^ 2 ) ) + ( RangeFit ( 2 ) .* fRModel ) + RangeFit ( 3 ) ; % Produce high resolution fitted spring contant values.
lLowerParams = polyfit ( fRRange , klLowers , 2 ) ; % Fits spring constant to resonant frequency for longest cantilever.
lLowerModel = ( lLowerParams ( 1 ) .* ( fRModel .^ 2 ) ) + ( lLowerParams ( 2 ) .* fRModel ) +  lLowerParams ( 3 ) ;
lUpperParams = polyfit ( fRRange , klUppers , 2 ) ; % Fits spring constant to resonant frequency for shortest cantilever.
lUpperModel = ( lUpperParams ( 1 ) .* ( fRModel .^ 2 ) ) + ( lUpperParams ( 2 ) .* fRModel ) + lUpperParams ( 3 ) ;
mModel = kModel ./ ( ( 2 * pi .* ( fRModel .* 1E3 ) ) .^ 2 ) ; % Produce high resolution fitted mass values for nominal cantilever length.
mMinModel = lLowerModel ./ ( ( 2 * pi .* ( fRModel .* 1E3 ) ) .^ 2 ) ; % Produce high resolution fitted mass values for longest cantilever.
mMaxModel = lUpperModel ./ ( ( 2 * pi .* ( fRModel .* 1E3 ) ) .^ 2 ) ; % Produce high resolution fitted mass values for shortes cantilever.

iP = 0 ; % Whether or not to plot.

if iP == 1
    
    plot ( fRRange , kRange , 'o' , fRModel , kModel , fRModel , lLowerModel , fRModel , lUpperModel ) ;
    set ( gca , 'FontSize' , 18 ) ;
    xlabel ( 'Resonant frequency (kHz)' , 'FontSize' , 18 ) ;
    ylabel ( 'Spring constant (N/m)' , 'FontSize' , 18 ) ;
    legend ( 'Quoted values' , 'Nominal cantilever length fit' , 'Longest cantilever fit' , 'Shortest cantilver fit bound' , 'Location' , 'NorthWest' ) ;

end

%% Find mass and spring constant of probe from cantilever resonant frequency.

[ fRClosest , fRIndex ] = min ( abs ( ( fRModel - ResonantFrequency ) ) ) ;

% Indices of parameters closest to input values.
mClosest = mModel ( fRIndex ) ; % Estimated mass value.
mMax = mMaxModel ( fRIndex ) ; % Maximum mass estimate.
mMin = mMinModel ( fRIndex ) ; % Minimum mass estimate.
kClosest = kModel ( fRIndex ) ; % Estimated spring constant value.
kMax = lUpperModel ( fRIndex ) ; % Maximum spring constant estimate.
kMin = lLowerModel ( fRIndex ) ; % Minimum spring constant estimate.

% Convenient ouputs.
EstimateArray = cell ( 9 , 2 ) ;

EstimateArray { 1 } = 'Mass (pg)' ;
EstimateArray { 2 } = 'Max mass (pg)' ;
EstimateArray { 3 } = 'Min mass (pg)' ;
EstimateArray { 4 } = 'k (N/m)' ;
EstimateArray { 5 } = 'Max k (N/m)' ;
EstimateArray { 6 } = 'Min k (N/m)' ;
EstimateArray { 7 } = 'F/V (uN/V)' ;
EstimateArray { 8 } = 'Max F/V (uN/V)' ;
EstimateArray { 9 } = 'Min F/V (uN/V)' ;
EstimateArray { 10 } = mClosest * 1E12 ;
EstimateArray { 11 } = mMax * 1E12 ;
EstimateArray { 12 } = mMin * 1E12 ;
EstimateArray { 13 } = kClosest ;
EstimateArray { 14 } = kMax ;
EstimateArray { 15 } = kMin ;
EstimateArray { 16 } = kClosest * DeflectionSensitivity * 1E-3 / 0.71 ;
EstimateArray { 17 } = kMax * DeflectionSensitivity * 1E-3 / 0.71  ;
EstimateArray { 18 } = kMin * DeflectionSensitivity * 1E-3 / 0.71  ;

end
