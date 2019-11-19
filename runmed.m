function gefiltert=runmed(ungefiltert,filterbreite,erster_wert,nan_handling);

% running median filter (MATLAB 5)
% gefiltert=runmed(ungefiltert,filterbreite,erster_wert,nan_handling)
% �bergabeparameter: ungefiltert:   ungefilterte Daten; Zeilen- od. Spaltenvektor; L�nge >= filterbreite
%                    filterbreite:  Breite des Medianfilters; zum Vektorrand hin wird die Breite bei Bedarf
%                                   verringert
%                    erster_wert:   erster zu berechnender Wert; die �brigen Randwerte werden mit NaNs 
%                                   aufgef�llt
%                    nan_handling:  0: befindet sich ein NaN innerhalb des Filterfensters, ist der gefilterte
%                                      Wert ebenfalls NaN
%                                   1: NaNs im Filterfenster werden ignoriert und die Berechnung mit den
%                                      verbleibenden Werten ausgef�hrt
% R�ckgabeparameter: gefiltert:     gefilterte Daten (Zeilenvektor)
%
% Jochen Ditterich, 9/97

clear gefiltert;

if ~isnumeric(ungefiltert)
   	error('Es k�nnen nur numerische Daten gefiltert werden!');
end;

if ~isreal(ungefiltert)
   error('Es k�nnen nur reelle Daten gefiltert werden!');
end;

[m n]=size(ungefiltert);

if n==1
   ungefiltert=ungefiltert';
end;

[m n]=size(ungefiltert);

if m~=1
   error('Es k�nnen nur Vektoren gefiltert werden!');
end;

if filterbreite<1
   error('Filterbreite zu klein!');
end;

if round(filterbreite)~=filterbreite
   error('Filterbreite mu� ganzzahlig sein!');
end;

if floor(filterbreite/2)*2==filterbreite
   error('Filterbreite mu� ungerade sein!');
end;

if n<filterbreite
   error('Vektor zu kurz!');
end;

if erster_wert>floor((n+1)/2)
   error('Erster Wert zu gro� gew�hlt!');
end;

if erster_wert<1
   error('Erster Wert zu klein gew�hlt!');
end;

if round(erster_wert)~=erster_wert
   error('Erster Wert mu� ganzzahlig sein!');
end;

if erster_wert>(filterbreite+1)/2
   error('Erster Wert zu gro� gew�hlt!');
end;

letzter_wert=n-erster_wert+1;

if nan_handling~=0 & nan_handling~=1
   error('Unzul�ssiger Wert f�r NaN-Handling!');
end;

% Randbereiche
if erster_wert>1
   gefiltert(1:erster_wert-1)=nan;
   gefiltert(letzter_wert+1:n)=nan;
end;

start_voll=(filterbreite+1)/2;
halbe_breite=(filterbreite-1)/2;

if nan_handling==0 % keine NaN-Korrektur

	% Bereich mit voller Filterbreite
	for i=start_voll:n-start_voll+1
   	gefiltert(i)=median(ungefiltert(i-halbe_breite:i+halbe_breite));
   end;
   
   % Zwischenbereiche
   if erster_wert<start_voll
      for i=erster_wert:start_voll-1
         gefiltert(i)=median(ungefiltert(1:2*i-1));
      end;
      for i=n-start_voll+2:letzter_wert
         gefiltert(i)=median(ungefiltert(2*i-n:n));
      end;
   end;
   
else % NaN-Korrektur
   
	% Bereich mit voller Filterbreite
	for i=start_voll:n-start_voll+1
      vektor=ungefiltert(i-halbe_breite:i+halbe_breite);
      nutzdaten=vektor(find(~isnan(vektor)));
      if length(nutzdaten)==0
         gefiltert(i)=nan;
      else
         gefiltert(i)=median(nutzdaten);
      end;
   end;
   
   % Zwischenbereiche
   if erster_wert<start_voll
      for i=erster_wert:start_voll-1
         vektor=ungefiltert(1:2*i-1);
         nutzdaten=vektor(find(~isnan(vektor)));
         if length(nutzdaten)==0
            gefiltert(i)=nan;
         else
            gefiltert(i)=median(nutzdaten);
         end;
      end;
      for i=n-start_voll+2:letzter_wert
         vektor=ungefiltert(2*i-n:n);
         nutzdaten=vektor(find(~isnan(vektor)));
         if length(nutzdaten)==0
            gefiltert(i)=nan;
         else
            gefiltert(i)=median(nutzdaten);
         end;
      end;
   end;
   
end;
