clear
clc

%% Studies options for introducing a parallel resistor in order to enable switching RRAM with CAFM.
% Useful in determining what the maximum voltages and currents might be for
% conductive AFM measurements where the sample resistance is not
% significantly greater than the CAFM probe/contact resistance, and voltage
% is dropped.

Rprobe = 12600 ; % Probe/source resistance.
Rprist = 1E9 ; % Pristine device resistance.
Ron = 4E2 ; % On state device resistance.
Roff = 15E3 ; % Off state device resistance.
Rparallel = 100 : 10 : 1E7 ;
Vin = 10 ; % Input voltage.

%% Voltage/current across device as a function of parallel resistance, for fixed device state resistances.

Rdutprist = 1 ./ ( ( 1 ./ Rprist ) + ( 1 ./ Rparallel ) ) ;
Rduton = 1 ./ ( ( 1 ./ Ron ) + ( 1 ./ Rparallel ) ) ;
Rdutoff = 1 ./ ( ( 1 ./ Roff ) + ( 1 ./ Rparallel ) ) ;

Vprist = Vin .* ( Rdutprist ./ ( Rprobe + Rdutprist ) ) ;
Von = Vin .* ( Rduton ./ ( Rprobe + Rduton ) ) ;
Voff = Vin .* ( Rdutoff ./ ( Rprobe + Rdutoff ) ) ;

Iprist = Vprist ./ Rdutprist ;
Ion = Von ./ Rduton ;
Ioff = Voff ./ Rdutoff ;

Idutprist = Vprist ./ Rprist  ;
Iduton = Von ./ Ron ;
Idutoff = Voff ./ Roff ;

figure ( 1 ) ; 

subplot ( 2 , 2 , 1 ) ;
semilogx ( Rparallel , Vprist , Rparallel , Von , Rparallel , Voff ) ;
legend ( 'Pristine' , 'On' , 'Off' , 'Location' , 'Northwest' ) ;
xlabel ( 'Rparallel (Ohm)' ) ;
ylabel ( 'Device voltage (V)' );
set ( gca , 'FontSize' , 16 ) ;

subplot ( 2 , 2 , 2 ) ;
loglog ( Rparallel , Rdutprist , Rparallel , Rduton , Rparallel , Rdutoff ) ;
legend ( 'Pristine' , 'On' , 'Off' , 'Location' , 'Northwest' ) ;
xlabel ( 'Rparallel (Ohm)' ) ;
ylabel ( 'DUT resistance (Ohm)' );
set ( gca , 'FontSize' , 16 ) ;

subplot ( 2 , 2 , 3 ) ;
loglog ( Rparallel , Iprist , Rparallel , Ion , Rparallel , Ioff ) ;
legend ( 'Pristine' , 'On' , 'Off' , 'Location' , 'West' ) ;
xlabel ( 'Rparallel (?)' ) ;
ylabel ( 'Total current (A)' );
set ( gca , 'FontSize' , 16 ) ;

subplot ( 2 , 2 , 4 ) ;
loglog ( Rparallel , Idutprist , Rparallel , Iduton , Rparallel , Idutoff ) ;
legend ( 'Pristine' , 'On' , 'Off' , 'Location' , 'East' ) ;
xlabel ( 'Rparallel (?)' ) ;
ylabel ( 'Device current (A)' );
set ( gca , 'FontSize' , 16 ) ;

%% Voltage/current across device as a function of changing device resistance with fixed parallel resistor values.

Rparalleloptions = [ 1E2 , 1E4 , 1E6 ] ;
Rdevicechange = ( 100 : 100 : 1E9 ) ;

Rparallellow = ( Rparalleloptions ( 1 ) .* Rdevicechange ) ./ ( Rparalleloptions ( 1 ) + Rdevicechange ) ;
Rparallelmid = ( Rparalleloptions ( 2 ) .* Rdevicechange ) ./ ( Rparalleloptions ( 2 ) + Rdevicechange ) ;
Rparallelhigh = ( Rparalleloptions ( 3 ) .* Rdevicechange ) ./ ( Rparalleloptions ( 3 ) + Rdevicechange ) ;

Vlow = ( Vin .* Rparallellow ) ./ ( Rprobe + Rparallellow ) ;
Vmid = ( Vin .* Rparallelmid ) ./ ( Rprobe + Rparallelmid ) ;
Vhigh = ( Vin .* Rparallelhigh ) ./ ( Rprobe + Rparallelhigh ) ;
Vnoparallel = ( Vin .* Rdevicechange ) ./ ( Rprobe + Rdevicechange ) ;

Ilow = Vlow ./ Rdevicechange ;
Imid = Vmid ./ Rdevicechange ;
Ihigh = Vhigh ./ Rdevicechange ;
Inoparallel = Vnoparallel ./ Rdevicechange ;

figure ( 2 ) ;

subplot ( 2 , 1 , 1 ) ;
semilogx ( Rdevicechange , Vlow , Rdevicechange , Vmid , Rdevicechange , Vhigh , Rdevicechange , Vnoparallel ) ;
legend ( '100 Ohm parallel' , '10 kOhm parallel' , '1 MOhm parallel' , 'No parallel' , 'Location' , 'West' ) ;
xlabel ( 'Device resistance (?)' ) ;
ylabel ( 'Device voltage (V)' ) ;
set ( gca , 'FontSize' , 16 , 'XDir' , 'Reverse' ) ;

subplot ( 2 , 1 , 2 ) ;
loglog ( Rdevicechange , Ilow , Rdevicechange , Imid , Rdevicechange , Ihigh , Rdevicechange , Inoparallel) ;
legend ( '100 Ohm parallel' , '10 kOhm parallel' , '1 MOhm parallel' , 'No parallel' , 'Location' , 'Northwest' ) ;
xlabel ( 'Device resistance (?)' ) ;
ylabel ( 'Device current (I)' );
set ( gca , 'FontSize' , 16 , 'XDir' , 'Reverse' ) ;
