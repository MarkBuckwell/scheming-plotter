function [ NFC , DataArray ] = CantileverThermalFileSelector ( iSet )

%% Opens user interface to select AFM cantilever tune data files to import for analysis.

if iSet == 0
    
    [ FileGroup , DataPath ] = uigetfile ( '*.txt' , 'DialogTitle' , 'Select thermal calibration data files' , 'MultiSelect' , 'on' ) ; % Gets file names and location.

elseif iSet == 1
    
    [ FileGroup , DataPath ] = uigetfile ( '*.txt' , 'DialogTitle' , 'Select thermal tune test data' , 'MultiSelect' , 'on' ) ; % Gets file names and location.

end

SingleFile = double ( ischar ( FileGroup ) ) ; % Used to account for the case where only a single file is selected.

if  SingleFile > 0

    NFC = 1 ;
    
else
    
    NFC = length ( FileGroup ) ;  % Number of files to import.
    
end

FileSet = repmat ( { '' } , 1 , NFC ) ; % Generates cell array to place filenames into.

cd ( DataPath ) ; % Sets target directory to current folder.


%% Generates array of files for analysis.

for i = 1 : NFC
    
    if SingleFile > 0
            
        FileSet ( i ) = cellstr ( strcat ( DataPath , FileGroup ) ) ; % Concatenates path and file strings and adds to output array,
    
    else
        
        FileSet ( i ) = strcat ( DataPath , FileGroup ( i ) ) ; % Concatenates path and file strings and adds to output array,
    
    end
    
end

Delimiter = '\t' ; % Could be specified in textscan to avoid variable creation.
FormatSpec = repmat ( '%f' , 1 , 2 ) ; % Specifies species labels import format based on number of species names in files.
StartRow = 2 ; % Specifies first row of data.

DataArray = zeros ( 41840 , 2 , NFC ) ; % Generates cell array to place data into.
SumArray = zeros ( 41840 , 2 , 1 ) ;

%% Opens, reads and closes data files.

for i = 1 : NFC

    % Chooses file.
    FileChoice = char ( FileSet ( : , i ) ) ;

    FileID = fopen ( FileChoice , 'r' ) ;
    FileArray = textscan ( FileID , FormatSpec , 'Delimiter' , Delimiter , 'EmptyValue' , NaN , 'HeaderLines' , StartRow - 1 , 'ReturnOnError' , false ) ;
    fclose ( FileID ) ;
    
    % Places array data into pre-defined recipient array.
    DataArray ( : , 1 , i ) = FileArray { 1 } ; % Column containing frequency data.
    DataArray ( : , 2 , i ) = FileArray { 2 } ; % Column containing intensity data
    SumArray ( : , 2 ) = SumArray ( : , 2 ) + FileArray { 2 } ; % Column containing summed intensity data.

end

if iSet == 1 
    
    SumArray ( : , 1 ) = FileArray { 1 } ; % Column containing summed corresponding frequency data.
    DataArray = SumArray ; % Switch output data to summed sampling data.
    
end

end
