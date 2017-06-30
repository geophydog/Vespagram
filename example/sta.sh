R=-170/-80/8/68
J=M7i
PS=sta.ps
PDF=sta.pdf

gmt grdcut etopo5.grd -R$R -Gcut.grd
gmt makecpt -T-7200/4500/10 -Z -Cwiki-scotland.cpt >tmp.cpt
gmt grdgradient cut.grd -A0/90 -Gcut_grad.grd -Ne0.7

gmt psxy -R$R -J$J -K -T>$PS
gmt pscoast -R -J -K -O -Bx10 -By10 -BWSen -A10000 -Gdarkyellow -Sskyblue >>$PS
gmt grdimage cut.grd -R -J -K -O -Icut_grad.grd -Ctmp.cpt>>$PS
gmt psxy -R -J -K -O -W2p,blue>>$PS<<EOF
-90.9488 13.7527
-158.628403 62.445450
EOF

saclst stlo stla f *.SAC | awk '{print $2,$3}' | gmt psxy -R -J -K -O -St0.7c -W1.5p>>$PS
echo -90.9488 13.7527 | gmt psxy -R -J -K -O -Sa0.5c -W2p,red>>$PS
gmt psxy -R -J -O -T>>$PS
gmt psconvert -Tg -A -P $PS
ps2pdf $PS $PDF
rm gmt.* tmp.*
evince $PDF
