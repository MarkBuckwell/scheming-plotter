clear
clc

%% This suite allows AFM probe characterisation data to be imported and processed.
% Sets of data files pertaining to probe resonant frequency, deflection
% sensitivity, thermal tune spectra... ADD MORE IN TIME... may be
% systematically imported and processed to output values for spring
% constant and effective mass, along with ranges and error values. In each
% case, one or more .txt data files may be input. Measurement scale should
% be in the first column and intensity data in the second. More than 2
% columns of data is ok, additional columns will (currently...) be ignored.

% Tapping-mode cantilver tuning data. The first column should be
% frequency data, the second column should be intensitiy.
[ TuneArray , TuneLocation ] = CantileverTuneProcessorFunction ( ) ;

ResonantFrequency = TuneArray { 7 } ;

% Contact-mode force-distance spectra. The first column
% should be z movement, the second column should be deflection intensity.
[ SensArray , SensLocation ] = CantileverDeflectionSensitivityProcessorFunction ( ) ;

DeflectionSensitivity = SensArray { 11 } ;
SensitivityError = SensArray { 14 } ;

% Contact-mode thermal tune data. The first column should be frequency
% data, the second column should be intensity. Requires value of deflection
% sensitivity used by software when aquiring data.
[ ThermalArray , ThermLocation ] = CantileverThermalProcessorFunction ( DeflectionSensitivity , SensitivityError ) ;

% Uses gathered data to estimate additional probe parameters.
[ EstimateArray ] = CantileverForceEstimatesFunction ( ResonantFrequency , DeflectionSensitivity ) ;

% Puts all information into a single array and show file locations to check
% all correspond to the same probe and measurement session.
FullArray = vertcat ( TuneArray , SensArray , EstimateArray , ThermalArray ) ;
openvar ( 'FullArray' ) ;
disp ( [ 'Resonant frequency data: ' , TuneLocation ] ) ;
disp ( [ 'Deflection sensitivity data: ' , SensLocation ] ) ;
disp ( [ 'Thermal tuning data: ' , ThermLocation ] ) ;
