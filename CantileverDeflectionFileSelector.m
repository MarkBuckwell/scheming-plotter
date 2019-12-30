function [ FlipArray , SensVals , kVals , DataHeaders , NFC , aH , aV , rH , rV ] = CantileverDeflectionFileSelector ( )

% Function to import AFM force-distance spectra as .txt files, exported
% from Nanoscope analysis. Could be generalised to any columnar data
% files. Function could be developed a little to read the column headers
% and provide options and clarity on what's done with each column. At
% the moment, it assumes the first 2 columns are the ramp data and subsequnt
% columns are deflection, friction and height sensor. Also assumes that
% alternating columns need inversion so all data starts from the same side
% of the plot.

for i = 1 : length ( findobj ( 'Type' , 'Figure' ) )
    
    figure ( i ) ;
    clf
    
end

% Select files.

[ FileGroup , DataPath ] = uigetfile ( '*.txt' , 'DialogTitle' , 'Select force-distance data files' , 'MultiSelect' , 'on' ) ; % Gets file names and location.

SingleFile = double ( ischar ( FileGroup ) ) ; % Used to account for the case where only a single file is selected.

if  SingleFile > 0

    NFC = 1 ;
    
    FileSet = strcat ( DataPath , FileGroup ) ;
    NFile = fopen ( FileSet , 'r' ) ; % Chosen file to determine number of columns, set to 1 to use first specified file.
 
else
    
    NFC = length ( FileGroup ) ;  % Number of files to import.
    
    FileSet = repmat ( { '' } , 1 , NFC ) ; % Generates cell array to place filenames into.

    %% Generates array of files for analysis.

    for i = 1 : NFC
    
        FileSet ( i ) = strcat ( DataPath , FileGroup ( i ) ) ; % Concatenates path and file strings and adds to output array,
    
    end
    
    NFile = fopen ( char ( FileSet ( : , 1 ) ) , 'r' ) ; % Chosen file to determine number of columns, set to 1 to use first specified file.

end

cd ( DataPath ) ; % Sets target directory to current folder.

%% Determines number of headerlines.

frewind ( NFile ) ;
HeaderArray = textscan ( NFile, '%s' , 'Delimiter' , '\n'  ) ;
HeaderArray = HeaderArray { 1 } ;
HeaderRow = find ( not ( cellfun ( 'isempty' , ( strfind ( HeaderArray ,...
    'Calc_Ramp' ) ) ) ) * 1 ) ; % Specifies the row containing headers.

%% Determines number of data columns.

Delimiter = '\t' ;
% Gets column headers. Assumes all files have same column and data format.
DataHeaders = textscan ( HeaderArray { HeaderRow } , '%s' , 'Delimiter' , Delimiter ) ;
DataHeaders = DataHeaders { 1 }  .' ;
frewind ( NFile ) ;
LArray = textscan ( NFile , '' , 'Delimiter' , Delimiter , 'EmptyValue' , NaN , 'HeaderLines' , HeaderRow , 'ReturnOnError' , false ) ;
N = size ( LArray , 2 ) - 1 ; % Specifies number of columns of data.
FormatSpec = repmat ( '%f' , 1 , N ) ; % Specifies import format; N cells of data columns.
Samples = length ( cell2mat ( LArray ( 1 ) ) ) ; % Number of samples for extend/retract.
fclose ( NFile) ;
    
%% Imports data.

DataArray = zeros ( Samples , N , NFC ) ; % Sets up recipient matrix for all data.
SensVals = zeros (  NFC , 1 ) ; % Matrix for deflection sensitivity values.
kVals = zeros ( NFC , 1 ) ; % Matrix for spring constant values.

for i = 1 : NFC
    
    if NFC == 1 
        
        FileChoice = char ( FileSet ) ; % Choose single file.
        
    else
        
        FileChoice = char ( FileSet ( : , i ) ) ; % Choose file from set.
        
    end
    
    % Opens and reads selected file data to temporary array.
    FileID = fopen ( FileChoice , 'r' ) ;
    
    SensLine = textscan ( FileID, '%s' , 'Delimiter' , '\n'  ) ;
    SensLine = SensLine { 1 } ;
    SensRow = find ( not ( cellfun ( 'isempty' , ( strfind ( SensLine ,...
    '@Sens. DeflSens' ) ) ) ) * 1 ) ; % Specifies the row containing deflection sensitivity.

    if isempty (SensRow ) == 0 % Only try to extraxt sensitivity and spring constant
        % values if those rows exist.
        
        SensVal = textscan ( SensLine { SensRow } , '%s' , 'Delimiter' , ' ' ) ;
        SensVals ( i ) = str2double ( SensVal { 1 }  ( 4 ) ) ;

        kRow = find ( not ( cellfun ( 'isempty' , ( strfind ( SensLine ,...
        'Spring Constant' ) ) ) ) * 1 ) ; % Specifies the row containing spring constant.
        kVal = textscan ( SensLine { kRow ( 1 ) } , '%s' ) ;
        kVal = kVal { 1 }  ;
        kVal{ 3 } ( end ) = [ ] ;
        kVals ( i ) = str2double ( kVal { 3 } ) ;
        
    end
    
    frewind ( FileID ) ;
    FileArray = textscan ( FileID , FormatSpec , 'Delimiter' , Delimiter , 'EmptyValue' , NaN , 'HeaderLines' , HeaderRow  , 'ReturnOnError' , false ) ;
    fclose ( FileID ) ;
    
    DataArray ( : , : , i ) = cell2mat ( FileArray ) ;
    
end

aH = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
    'Calc_Ramp_Ex_nm' ) ) ) ) * 1 ) ; % Specifies the column containing height extend data.
rH = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
    'Calc_Ramp_Rt_nm' ) ) ) ) * 1 ) ; % Specifies the column containing height retract data.
aV = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
    'Defl_V_Ex' ) ) ) ) * 1 ) ; % Specifies the column containing deflection extend data.
rV = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
    'Defl_V_Rt' ) ) ) ) * 1 ) ; % Specifies the column containing deflection retract data.
aF = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
    'Friction_V_Ex' ) ) ) ) * 1, 1 ) ; % Specifies the column containing deflection extend data.
rF = find ( not ( cellfun ( 'isempty' , ( strfind ( DataHeaders ,...
    'Friction_V_Rt' ) ) ) ) * 1 ) ; % Specifies the column containing deflection retract data.

FlipArray = DataArray ; % Creates array to process data, leaving raw data untouched, though it is not passed to script.
FlipArray ( : , rH , : ) = flipud ( DataArray ( : , rH , : ) ) ; % Inverts height retract data.
FlipArray ( : , rV , : ) = flipud ( DataArray ( : , rV , : ) ) ; % Inverts voltage retract data.

if isempty ( aF ) == 0
    
    FlipArray ( : , rF , : ) = flipud ( DataArray ( : , rF , : ) ) ; % Inverts frictionS retract data.

end

%% Option to baseline-correct data. Could be expanded/improved by adding a
% listdlg to give user option of which data types to correct, in addition
% to ensuring correct column indexing.

BackgroundSize = size ( FlipArray , 1  ) / 2 ;
CorrQuest = questdlg ( 'Non-contact approach and retract regions will be separately fitted and the corresponding linear background subtracted from each dataset.' , ...
    'Perform linear baseline correction?' , 'On extend data' , 'On retract data' , 'Do not correct' ,'On retract data' ) ;
    
if strcmp ( CorrQuest , 'On extend data' ) == 1
    
    for i = 1 : NFC
              
        % Deflection data.
        RetractDBackground = polyfit ( FlipArray ( 10 : BackgroundSize , aH , i ) , FlipArray ( 10 : BackgroundSize , aV , i ) , 1 ) ;
        BackgroundDFit = ( FlipArray ( : , aH , i ) * RetractDBackground ( 1 ) ) + RetractDBackground ( 2 ) ;
        
        figure ( 5 ) ;
        
        subplot ( 2 , 2 , 1 ) ;
        plot ( FlipArray ( : , aH , i ) , FlipArray ( : , aV , i ) , FlipArray ( : , aH , i ) , BackgroundDFit ) ;
        hold on
        
        FlipArray ( : , aV , i ) = FlipArray ( : , aV , i ) - BackgroundDFit ;
        FlipArray ( : , rV , i ) = FlipArray ( : , rV , i ) - BackgroundDFit ;
        
        subplot ( 2 , 2 , 3 ) ;
        plot ( FlipArray ( : , aH , i ) , FlipArray ( : , aV , i ) ) ;
        hold on
        
        if isempty ( aF ) == 0
            
            % Friction data.
            RetractFBackground = polyfit ( FlipArray ( 10 : BackgroundSize , aH , i ) , FlipArray ( 10 : BackgroundSize , aF , i ) , 1 ) ;
            BackgroundFFit = ( FlipArray ( : , aH , i ) * RetractFBackground ( 1 ) ) + RetractFBackground ( 2 ) ;

            subplot ( 2 , 2 , 2 ) ;
            plot ( FlipArray ( : , aH , i ) , FlipArray ( : , aF , i ) , FlipArray ( : , aH , i ) , BackgroundFFit ) ;
            hold on

            FlipArray ( : , aF , i ) = FlipArray ( : , aF , i ) - BackgroundFFit ;
            FlipArray ( : , rF , i ) = FlipArray ( : , rF , i ) - BackgroundFFit ;

            subplot ( 2 , 2 , 4 ) ;
            plot ( FlipArray ( : , aH , i ) , FlipArray ( : , aF , i ) ) ;
            
            hold on
        
        end
        
    end
    
elseif strcmp ( CorrQuest , 'On retract data' ) == 1
    
    for i = 1 : NFC
    
        % Deflection data.
        RetractDBackground = polyfit ( FlipArray ( 10 : BackgroundSize , rH , i ) , FlipArray ( 10 : BackgroundSize , rV , i ) , 1 ) ;
        BackgroundDFit = ( FlipArray ( : , rH , i ) * RetractDBackground ( 1 ) ) + RetractDBackground ( 2 ) ;
        
        figure ( 5 ) ;
        
        subplot ( 2 , 2 , 1 ) ;
        plot ( FlipArray ( : , rH , i ) , FlipArray ( : , rV , i ) , FlipArray ( : , rH , i ) , BackgroundDFit ) ;
        hold on
        
        FlipArray ( : , aV , i ) = FlipArray ( : , aV , i ) - BackgroundDFit ;
        FlipArray ( : , rV , i ) = FlipArray ( : , rV , i ) - BackgroundDFit ;
        
        subplot ( 2 , 2 , 3 ) ;
        plot ( FlipArray ( : , rH , i ) , FlipArray ( : , rV , i ) ) ;
        hold on
        
        if N == 6
            
            % Friction data.
            RetractFBackground = polyfit ( FlipArray ( 10 : BackgroundSize , rH , i ) , FlipArray ( 10 : BackgroundSize , rF , i ) , 1 ) ;
            BackgroundFFit = ( FlipArray ( : , rH , i ) * RetractFBackground ( 1 ) ) + RetractFBackground ( 2 ) ;

            subplot ( 2 , 2 , 2 ) ;
            plot ( FlipArray ( : , rH , i ) , FlipArray ( : , rF , i ) , FlipArray ( : , rH , i ) , BackgroundFFit ) ;
            hold on

            FlipArray ( : , aF , i ) = FlipArray ( : , aF , i ) - BackgroundFFit ;
            FlipArray ( : , rF , i ) = FlipArray ( : , rF , i ) - BackgroundFFit ;

            subplot ( 2 , 2 , 4 ) ;
            plot ( FlipArray ( : , rH , i ) , FlipArray ( : , rF , i ) ) ;
            hold on
    
        end
            
    end
    
else
    
    disp ( 'No baseline correction.' ) ;
    
end

%% To account for files in which the extend and retract data lengths are not
% equal (i.e. some data is missing, usually from the start of the extend
% curve, so this bit of code is designed for that case), the missing extend
% data is taken from the retract data, since that region is just the
% non-contact regime so is not relevant to analysis and should be flat
% (particuarly once baseline has been subtracted).
for i = 1 : NFC 
   
    % Find index of first NaN in extend ramp data.
    NanPos = find ( isnan ( FlipArray ( : , aH , i ) ) , 1 ) ;
    
    % Process if NaN elements exist.
    if isempty ( NanPos ) == 0
        
        % Replace extend ramp data with retract data, rather than trying to
        % mimick small variations in increments between datapoints.
        FlipArray ( : , aH , i ) = FlipArray ( : , rH , i ) ;
        
        % Shift extend deflection data to end of column.
        FlipArray ( : , aV , i ) = circshift ( FlipArray ( : , aV , i ) ,...
            size ( FlipArray , 1 ) - NanPos + 1 , 1 ) ;
        % Replace with retract data.
        FlipArray ( 1 : size ( FlipArray , 1 ) - NanPos + 1 , aV , i ) =...
            FlipArray( 1 : size ( FlipArray , 1 ) - NanPos + 1 , rV , i ) ;
        
        % Shift extend friction data to end of column.
        FlipArray ( : , aF , i ) = circshift ( FlipArray ( : , aF , i ) ,...
            size ( FlipArray , 1 ) - NanPos + 1 , 1 ) ;
        % Replace with retract data.
        FlipArray ( 1 : size ( FlipArray , 1 ) - NanPos + 1 , aF , i ) =...
            FlipArray ( 1 : size ( FlipArray , 1 ) - NanPos + 1 , rF , i ) ;        
        
    end
    
end
    
end
