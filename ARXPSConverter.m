clear
clc

%% Import and reformat ARXPS data files to split out all angles and depths
% so that they import more readily into CasaXPS.

%% Choose data files, which should be .avg files.

addpath ( cd ) ;

[ FileGroup , DataPath ] = uigetfile ( '*.avg' , 'DialogTitle' ,...
    'Select Data Files' , 'MultiSelect' , 'on' ) ;
% Gets file names and location.

SingleFile = double ( ischar ( FileGroup ) ) ;
% Used to account for the case where only a single file is selected.

if  SingleFile > 0

    NFC = 1 ;
    
else
    
    NFC = length ( FileGroup ) ;  % Number of files to import.
    
end

FileSet = repmat ( { '' } , 1 , NFC ) ; % Generates cell array to place filenames into.

cd ( DataPath ) ; % Sets target directory to current folder.

iFI = 1 ; % File importation index,

while iFI <= NFC
    
    if SingleFile > 0
        
        % Concatenates path and file strings and adds to output array,
        FileSet ( iFI ) = cellstr ( strcat ( DataPath , FileGroup ) ) ;
        
    else
        
        % Concatenates path and file strings and adds to output array,
        FileSet ( iFI ) = strcat ( DataPath , FileGroup ( iFI ) ) ;
        
    end
    
    iFI = iFI + 1 ; % Increments importation index.

end

%% Import cell arrays from data files, just taking in each row as a cell.

FormatSpec = '%q%[^\n\r]' ;
Delimiter = '' ; % No delimiting, just import each row individually.
ScanArray = repmat ( '' , NFC ) ;
SamplesArray = zeros ( 1 , NFC ) ; % Array for number of energy samples in spectrum.
AnglesArray =  SamplesArray ; % Array for number of angles per sample.
DepthArray = SamplesArray ;  % Array for number of depth levels sampled.

for i = 1 : NFC
   
    % Chooses file.
    FileChoice = char ( FileSet ( : , i ) ) ;
    
    FileID = fopen ( FileChoice , 'r' ) ;
    FileArray = textscan ( FileID , FormatSpec , 'Delimiter' , Delimiter ,...
        'EmptyValue' , NaN , 'ReturnOnError' , false ) ;
    fclose ( FileID ) ;
    
    ScanArray { i } = [ FileArray{ 1 : end - 1 } ] ;
    
    % Get sampling parameters; number of samples (energy channels), angles, depths.
    SamplesCheckCell = find ( not ( cellfun ( 'isempty' ,...
        ( strfind ( ScanArray {  i }  , '$SPACEAXES=4' ) ) ) ) ) ;
    
    SamplesCheck = textscan ( ScanArray { i } { SamplesCheckCell + 1 } ,...
        '%s' , 'Delimiter' , ',' ) ; % Reads sampling config line of data file.
    SamplesArray ( i ) =  str2double ( SamplesCheck { 1 } { 3 } ) ;
    AnglesCheck =  textscan ( ScanArray { i } { SamplesCheckCell + 2 } ,...
        '%s' , 'Delimiter' , ',' ) ; % Reads angles config line of data file.
    AnglesArray ( i ) =  str2double ( AnglesCheck { 1 } { 3 } ) ;
    DepthCheck =  textscan ( ScanArray { i } { SamplesCheckCell + 4 } ,...
        '%s' , 'Delimiter' , ',' ) ; % Reads depth config line of data file.
    DepthArray ( i ) =  str2double ( DepthCheck { 1 } { 3 } ) ;    
    
end

SeparateRegions = questdlg ( 'Export regions (files) to separate folders?' ,...
    'Region separation.' , 'Yes' , 'No' , 'Other' , 'No' ) ;

if strfind ( SeparateRegions , 'No' ) == 1
        
    cd ( DataPath ) ;
    mkdir ( 'Output files' ) ;
    NewDataPath = strcat ( DataPath , '/Output files' ) ;
    cd ( NewDataPath ) ;
    
end        

%% Go through each data file, condensing and separating into separate files
% and save pre-written vb sum funtion for use in CasaXPS. This will sum all
% angles at a given depth to give a sum spectrum, as if no ARXPS has been used.
for i = 1 : NFC % Step through elements/region spectra.
    
    DataStep = ceil ( SamplesArray ( i ) / 4 ) + 5 ;
    % Number of columns to step through for each spectrum.
    
    % Generates CasaXPS Sample Identifier prefix.
    SampleTitle = textscan ( ScanArray { i } { 2 } , '%s' , 'Delimiter' , '\' ) ;
    SampleTitle { 1 } { end - 1 } = strcat ( SampleTitle { 1 } { end }...
        ( 1 : end - 9 ) , ' DP' ) ;
    SampleIdentifier = SampleTitle { 1 } { 1 } ;
    
    for h = 2 : length ( SampleTitle { 1 } )
        
        SampleIdentifier = strcat ( SampleIdentifier , '\' ,...
            SampleTitle { 1 } { h } ) ;
        
    end
    
    if strfind ( SeparateRegions , 'Yes' ) == 1
        
        cd ( DataPath ) ;
        mkdir ( SampleTitle { 1 } { 8 } ) ;
        NewDataPath = strcat ( DataPath , '/' , SampleTitle { 1 } { 8 } ) ;
        cd ( NewDataPath ) ;
        
   end
    
    ScanArray { i } { 2 } = SampleIdentifier ;
    
    %% Removes depth increment from data label.
    DataIncrement = textscan ( ScanArray { i } { SamplesCheckCell + 8 +...
        ( AnglesArray ( i ) * DataStep ) } , '%s' , 'Delimiter' , '=' ) ;
    DataIncrement = DataIncrement { 1 } ( end ) ;
    DataIncrement = DataIncrement { 1 } ( 1 : 5 ) ;
    NoIncrement = strfind ( ScanArray { i } { SamplesCheckCell + 3 } , DataIncrement ) ;
    
    if DataIncrement ( 2 ) == '.'
        
        ScanArray { i } { SamplesCheckCell + 3 } ( NoIncrement : NoIncrement + 7 ) = '0.000000' ;
        
    elseif DataIncrement ( 3 ) == '.' 
        
        ScanArray { i } { SamplesCheckCell + 3 } ( NoIncrement : NoIncrement + 8 ) = ' 0.000000' ;
        
    elseif DataIncrement ( 4 ) == '.'
        
        ScanArray { i } { SamplesCheckCell + 3 } ( NoIncrement : NoIncrement + 9 ) = '  0.000000' ;        
        
    end
    
    for k = 1 : AnglesArray ( i ) % Step through angles, 1 per file.
        
        % Generates CasaXPS BlockIDs correlated with angle for each file.
        AnglesTitles = textscan ( ScanArray { i } { ( SamplesCheckCell + 7 +...
            ( ( k - 1 ) * DataStep ) ) } , '%s' , 'Delimiter' , '=' ) ;
        AngleName = ( AnglesTitles { 1 } { 7 } )  ;
        AngleTitle = num2str ( round ( str2num ( AngleName ) , 1 ) ) ;
        % Should stay as str2num, not str2double, as last character is ';'.
        SpectrumTitle = textscan ( ScanArray { i } { 9 } , '%s' , 'Delimiter' , '=' ) ;
        SpectrumTitle { 1 } { 2 } = strcat ( '''' , AngleTitle , '''' ) ;
        ScanArray { i } { 9 } = strcat ( SpectrumTitle { 1 } { 1 } , '=' ,...
            SpectrumTitle { 1 } { 2 } ) ;
        
        DataOutput = ScanArray { i } ;
        StartLine = 88 + ( DataStep * ( k - 1 ) ) ;
        
        % Remove data up to current angle at first depth (no removal for first angle).
        % An error will occur if the settings for any region were changed
        % during measurements.
        if k > 1 
                
            DataOutput ( SamplesCheckCell + 6 : StartLine - 1 ) = [] ;
            IncrementTransfer = textscan ( DataOutput { SamplesCheckCell + 10 } , '%s' , 'Delimiter' , ',' ) ;
            IncrementTransfer { 1 } { 2 } = '0' ;
            DataOutput { SamplesCheckCell + 10 } = strcat ( IncrementTransfer { 1 } { 1 } , ',' ,...
                IncrementTransfer { 1 } { 2 } , ',' , IncrementTransfer { 1 } { 3 } ) ;

        end
        
        for l = 1 : DepthArray ( i ) % Condense depths for each file.
            
            FromLine = SamplesCheckCell + 6 + ( DataStep * l ) ; % First line of data to remove.
            
            if FromLine + ( DataStep * ( AnglesArray ( i ) ) - 1 )...
                    <= length ( DataOutput ) 
                
                DataOutput ( FromLine : FromLine + ( DataStep *...
                    ( AnglesArray ( i ) - 1 ) - 1 ) ) = [ ] ;

                IncrementTransfer = textscan ( DataOutput { SamplesCheckCell + 10 +...
                    ( DataStep * l ) } , '%s' , 'Delimiter' , ',' ) ;
                IncrementTransfer { 1 } { 2 } = '0' ;
                DataOutput { SamplesCheckCell + 10 + ( DataStep * l ) } = strcat...
                    ( IncrementTransfer { 1 } { 1 } , ',' ,...
                    IncrementTransfer { 1 } { 2 } , ',' ,...
                    IncrementTransfer { 1 } { 3 } ) ;
 
            else
                
                DataOutput ( FromLine : end ) = [ ] ;
                                
            end
                             
        end

        %% Save individual .avg files.
        FileName = strcat ( SampleTitle { 1 } { 8 } , ',' , { ' ' } , AngleTitle ) ;
        dlmcell ( strcat ( FileName { 1 } , '.avg' ) , DataOutput ) ;
        
    end
    
    %% Create separate .txt files for CasaXPS VB summing function.
    if strfind ( SeparateRegions , 'Yes' ) == 1
    
    VBArray = repmat ( { '' } , DepthArray ( i ) , 1 ) ;
    
    for k = 1 : DepthArray ( i )
        
        SumLine = [ ] ;
        
        for l = 1 : AnglesArray ( i )
            
            SumLine = strcat ( SumLine , 'vb' , num2str ( ( ( l - 1 ) *...
                DepthArray ( i ) ) + k - 1 ) , '+' ) ;
            
        end
   
        SumLine ( end ) = [ ] ;
        VBArray { k } = SumLine ;
        
    end
    
    SumFile = cell2table ( VBArray ) ;
    writetable ( SumFile ,  strcat ( SampleTitle { 1 } { 8 } , ' sum' ) ) ;
    
    end
    
end

%% Create single, combined .txt file for CasaXPS VB summing function.
if strfind ( SeparateRegions , 'No' ) == 1
    
    cd ( NewDataPath ) ;
    
    VBArray = repmat ( { '' } , sum ( DepthArray ) , 1 ) ;
    VBInteger = 0 ;
    
    for i = 1 : size ( DepthArray , 2 )
        
        for k = 1 : DepthArray ( i )
            
            SumLine = [ ] ;
            
            for l = 1 : AnglesArray ( i )

                SumLine = strcat ( SumLine , 'vb' , num2str ( ( ( l - 1 ) .*...
                    DepthArray ( i ) ) + k - 1 + VBInteger ) , '+' ) ;

            end
            
            SumLine ( end ) = [ ] ;
            VBArray { ( ( i - 1 ) .* ( ( DepthArray ( i ) ) ) ) + k } = SumLine ;
        
        end
        
        VBInteger = VBInteger + ( DepthArray ( i ) .* AnglesArray ( i ) ) ;
   
    end
        
    SumFile = cell2table ( VBArray ) ;
    AllTitle = textscan ( DataPath , '%s' , 'Delimiter' , '\' ) ;
    AllTitle = AllTitle { 1 } ;
    VBTitle = strcat ( AllTitle { end - 1 } , { ' ' } , AllTitle { end } ,...
        { ' ' } , 'VB summation equation' ) ;    
    writetable ( SumFile , VBTitle { 1 } ) ;
    winopen ( strcat ( VBTitle { 1 } , '.txt' ) ) ;
         
end

disp ( DataPath ) ;
