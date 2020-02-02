clear
clc
set ( 0 , 'DefaultFigureWindowStyle' , 'Docked' ) ;

%% To model the electrical behaviour of a CAFM measurement.
% A WORK IN PROGRESS!!!
% This script can be used to asses the Cooulombic interaction between a CAFM
% probe and sample, as well as model the shape of a conductive feature that
% has been convoluted by the geometry of the scanning probe apex.

%% Firstly, looking at the Coulombic interaction between a CAFM probe and sample.
% This is used for thinking about the force that the probe will experience
% as a function of the applied voltage. Thus, considerations may be made on
% experimental observations, wherein there is some uncertainty on whether a
% feature is real, or the result of deflection/Coulombic forces.
% Assuming charges on each object are opposite (they are ideal conductors).
% Assuming circular contact area and sharp probe.
k = 1 / ( 4 * pi * 8.85418782E-12 ) ; % Coulomb constant.
d = ( 0.01 : 0.01 : 30 ) * ( 1E-9 ) ; % Separation in nm.
A = pi * ( 10E-9 ^ 2 ) ; % Contact area in nm^2.
Eair = 8.85418782E-12 ; % Permitivity of air.
ESiO2 = 3.9 * Eair ;
C = ESiO2 .* ( A ./ d ) ; % Capacitance of 'air gap'.

for i = 1 : 10
    
    V ( i ) = i ;
    
    Q = C * V ( i ) ;

    F = ( k * ( Q .^ 2 ) ) ./ ( d .^ 2 ) ; % Force between probe and sample.
    x { i } = F ./ 40 ; % Deflection of cantilever, given spring constant of 40 N/m.
    
    figure ( 2 ) ;
    semilogy ( ( d .* 1E9 ) , ( x { i } .* 1E6 ) , 'DisplayName' , strcat ( num2str ( V ( i ) ) , 'V' ) ) ;
    hold on
    
    figure ( 1 ) ;
    semilogy ( ( d .* 1E9 ) , F , 'DisplayName' , strcat ( num2str ( V ( i ) ) , 'V' ) ) ;
    hold on
    
end

xlabel ( 'Separation/nm' ) ;
ylabel ( 'Force/N' ) ;
ylim ( [ 1E-10 3E0 ] ) ;
legend show
nN = ones ( length ( d ) ) .* 1E-9 ;
halfnm = ones ( length ( F ) ) ./ 2 ;
plot ( ( d .* 1E9 ) , nN , halfnm , F ) ;

figure ( 2 ) ;
xlabel ( 'Separation/nm' ) ;
ylabel ( 'Deflection/um' ) ;
ylim ( [ 1E-6 1E6 ] ) ;
legend show
nm = ones ( length ( d ) ) .* 1E-3;
halfnm = ones ( length ( F ) ) ./ 2 ;
plot ( ( d .* 1E9 ) , nm , halfnm , x { i } ) ;

%% Now to model to the Lennard-Jones potential in order to convolute it with the Coulombic.

%% Now to model the normalised current (i.e. some arbitrary current) as a function
% of the tip-sample contact area. Initial, simple model is that the
% current scales linearly with the contact area (based on the equations in
% Mario's book).

D = 12 ; % Tip diameter in pixels. There might be a mismatch between this and
% the actual array size, as 1 must be added on to give a centre.
R = 16 ; % Feature width/diameter in pixels.

% Initially, assume tip and sample conductive feature are both flat and
% circular. Perform set of line scans across feature, wherein each pixel
% constitutes change in overlap between tip and sample, with the feature
% centred at the centre of the imaging area.

LineScale = 0.5 ; % Scalar for scan size. Set to 1 to give 1 nm/line and sample.
Samples = 256 ; % Resolution of image. 0 added in making array for definite
% mid point, along with additional tip diameter, in order for the tip to be
% able to scan up to having its centre at the edge of the image.
Resolution = 1 ; % Additional resolution scaler.

X = ( 0 : Resolution : Samples + D ) .* LineScale ; % x-axis lines.
Y = ( 0 : Resolution : Samples + D ) .* LineScale ; % y-axis samples.

% OverlapArray = zeros ( size ( X , 2 ) , size ( Y , 2 ) ) ; % Create an empty
% array for the overlap between tip and feature.

FeatureArray = zeros ( size ( X , 2 ) , size ( Y , 2 ) ) ; % Create an empty array for the sample feature.
FeatureArray ( ( ( Samples + D ) / ( 2 * Resolution ) ) + 1 ,...
    ( ( Samples + D ) / ( 2 * Resolution ) ) + 1 ) = 1 ; % Put a 1 at the centre of the image.
FeatureLogical = bwdist ( FeatureArray ) ; % Array of distances from centre of feature.
Feature = FeatureLogical <= R / 2 ; % Logical array where 1s denote the feature, whose radius is R / 2.

TipArray = zeros ( ( D / Resolution ) + 1 , ( D / Resolution ) + 1 ) ; % Create and empty array for the probe.
TipArray ( ( D / ( 2 * Resolution ) ) + 1 , ( D / ( 2 * Resolution ) ) + 1 ) = 1 ; % Put a 1 at the centre of the image.
TipLogical = bwdist ( TipArray ) ; % Array of distances from centre of tip.
Tip = TipLogical <= D / 2 ; % Logical array where 1s denote the tip, whose radius is D / 2.

PixelsSum = zeros ( Samples / Resolution , Samples / Resolution ) ; % Array for sum of nonzero elements
% when the tip is at each pixel in an image of size Samples x Samples.

for i = 1 : Samples / Resolution % Increment along slow axis (y).
   
    for j = 1 : Samples / Resolution % Increment along fast axis (x).
       
        % Map tip array onto feature arry, centred at pixel i , j, and compute
        % multiple, so any overlapping 1s will give 1, anywhere with a 0
        % gives 0.
        %% Uncomment this section to watch the CAFM probe scan the feature!
%         ScanArray = Feature ; % Array to check what the scan is doing.
%         ScanArray ( i : i + D , j : j + D ) =...
%             Feature ( i : i + D , j : j + D ) + Tip ;

%         imagesc ( ScanArray ) ;
%         pause ( 0.0001 ) ;

%%         OverlapArray = OverlapArray .* 0 ; % Needs to be reset, otherwise
        % the tip area is continuously written into the array as non-zero
        % elements and thus contributes to the overlap.
        
%         OverlapArray ( i : i + D , j : j + D ) =...
%             Feature ( i : i + D , j : j + D ) .* Tip ;

%         if sum ( sum ( OverlapArray ) ) > 0
%             
%             imagesc ( OverlapArray ) ;
%             pause ( 0.0001 ) ;
% 
%         end
 
%         PixelsSum ( i , j ) = sum ( sum ( OverlapArray ) ) ;

        % OverlapArray isn't needed, sum can be computed directly into
        % output instead.
        PixelsSum ( i , j ) = sum ( sum ( Feature ( i : i + ( D / Resolution ) ,...
            j : j + ( D / Resolution ) ) .* Tip ) ) ;
       
    end
    
        % To print out the location of the artefact when OvelapArray is not
        % reset.
%         if PixelsSum ( i , 1 ) >= 1
%             
%             disp ( 'artifact when i = ' ) ;
%             disp ( i )  ;
%             
%         end
        
end

% figure ;
% imagesc ( PixelsSum ) ;

hmax = D / 2 ;
% xpos = 1 - sin
CentreLine = ( PixelsSum ( ( Samples / 2 ) + 1 , : ) )' ;
Background = 0.016 ; % pA background current
BackgroundRatio = ( Background /12.28 ) * max ( CentreLine ) ;
CentreLineNorm = ( CentreLine + BackgroundRatio ) ./ max (CentreLine ) ;

xpos = CentreLine ;
figure ;
plot ( CentreLineNorm ) ;

open CentreLineNorm

% figure ( 2 ) ;
% 
% for i = 1 : 10 : Samples / ( 2 * Resolution )
%     
%     plot ( PixelsSum ( i , : ) ) ;
%     hold on
% 
% end
