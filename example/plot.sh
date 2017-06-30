R=1600.000000/1800.000000/50.000000/150.000000
J=X9i/6i
PS=1600.000000~1800.000000.ps
PDF=1600.000000~1800.000000.pdf
gmt surface out.txt -R$R -I0.250000/0.500000 -Gout.txt.grd
gmt grd2cpt out.txt.grd -Cjet -Z>tmp.cpt
gmt psxy -R$R -J$J -K -T>$PS
gmt grdimage out.txt.grd -R -J -K -O -Bx20+l"Time(sec)" -By15.000000+l"Backazimuth(deg)" -BWSen -Ctmp.cpt>>$PS
gmt psscale -Ctmp.cpt -R -J -K -O -D9.4i/3i/12/0.8 -Ba0.4:"Spectral Power":>>$PS
gmt psxy -R -J -O -T>>$PS
ps2pdf $PS $PDF
gmt psconvert -Tg -A -P $PS
rm tmp.* gmt.*
evince $PDF
