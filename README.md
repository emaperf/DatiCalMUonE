La cartella AngleCompare contiene:
il foglio Excel Angle.ods in cui sono riportati tutti i dati relativi agli angoli valutati considerando una matrice 3x3, 5x5 e 7x7 attorno al cristallo col massimo e i grafici di delta e sigma delta riferiti alle matrici 3x3 e 5x5 per un confronto
 
La cartella ShowerSize1 in cui ci sono tutti i plots, divisi per target (inseriti per confrontare i grafici degli angoli considerando diverse matrici attorno al cristallo col massimo)

La cartella ShowerSize2_s1+s2 contiene:
1) i Plots, divisi in base al target urtato (T1, T20, T40), in cui sono riportati i plots di theta, delta theta, E, Erel; 
theta e delta theta sono ottenuti considerando una matrice 5x5 attorno al cristallo con massima energia; 
E, Erel presentano la sovrapposizione di tre grafici:
- in nero, la somma delle energie su tutti i cristalli, senza smearing
- in verde, l'energia dei cristalli interni con lo smearing dell'elettronica e del termine fisso (s1+s2)
- in rosso, l'energia di cristalli con lo smearing nel blocco 5x5 attorno al massimo;
Inoltre nelle cartelle T1,T20,T40 sono riportati, nei file CrD.., il numero dei cristalli interni illuminati.

2) Nella cartella Grafici i file:
- Delta raffiguranti i grafici di Delta e sigmaDelta in funzione di E
- Erel raffiguranti i grafici di Erel in funzione di E e di ErelSmeared in funzione di ESmeared (quest'ultimo considera       l'energia con lo smear dei grafici verdi dei plot)
- RMS raffiguranti la sigma di Erel e la sigma di ErelSmeared rispettivamente in funzione di E e di ESmeared (sempre relativa al grafico verde)

3) Foglio Excel MUonE in cui Smeared (5x5) è riferito al grafico rosso mentre Smeared (whole) al grafico verde delle energie

La cartella ShowerSize3 contiene:
1) Le cartelle T1 T20 con i plot di theta, delta theta, E, Erel in cui sia gli angoli che le energie sono stati valutati considerando una matrice 7x7 attorno al cristallo col massimo di energia. Lo smearing (grafico rosso) non considera il termine fisso s2
2) Il foglio Excel MUonE con i dati relativi ai plot

La cartella ShowerSize2 contiene:
1) Le cartelle T1 T20 T40 con i plot di theta, delta theta, E, Erel in cui sia gli angoli che le energie sono stati valutati considerando una matrice 5x5 attorno al cristallo col massimo di energia. Lo smearing (grafico rosso) non considera il termine fisso s2

Smear.C è il codice utilizzato per i plots di ShowerSize2 e ShowerSize3

Smear2.C è il codice utilizzato per i plots di ShowerSize2_s1+s2
