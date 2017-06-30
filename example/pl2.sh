R=600/800/50/150
J=X9i/7i
PS=plot.ps
PDF=plot.pdf

gmt surface out.txt -R$R -I1/0.2 -Gtmp.grd
gmt grd2cpt tmp.grd -Cjet >tmp.cpt


gmt psxy -R$R -J$J -K -T>$PS
gmt grdimage tmp.grd -R -J -K -O -Ctmp.cpt -Bx20+l"Time(sec)" -By10+l"Backazimuth(deg)" -BWSen >>$PS
#gmt grdcontour tmp.grd -R -J -K -O -C0.2 -Bx100+l"Time(sec)" -By2+l"Slowness(sec/deg)" -BWSen >>$PS

gmt psxy -R -J -O -T>>$PS

ps2pdf $PS $PDF
rm tmp.* gmt.*
evince $PDF
