/*--------------------------------------------------------------------------------
 * This C programm is to execute Velocity Spectral Analysis (Vespagram or Vespa) *
 * 2017-6-29  Initially coded by Xuping Feng @NJU                                *
 * ------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "sacio.h"

#define PI 3.14159265

/*---------------------------------------------------to exclude space or new line chracter of a string-------------------------------------------------------*/
void no_spa(char *ps) {
    char *pt = ps;
    while ( *ps != '\0' ) {
        if ( *ps != ' ' && *ps != '\n' ) {
            *pt++ = *ps;
        }
        ++ps;
    }
    *pt = '\0';
}

float sign( float x ) {
    if ( x >= 0. ) return 1.;
    else return -1.;
}

/*----------------------------------------------------------------------------MAIN PROGRAM-------------------------------------------------------------------*/
int main( int argc, char *argv[] ) {
    int i, j, size = 256, sac_npts, shift_index, count = 0, sta_index = 0, begin_index, end_index, beam_npts;
    float *data, t1, t2, f1, f2, Nth_root, time, center_lon = 0., center_lat = 0., dx, dy, shift_time, delta,\
          time_start, time_end, **allamp, **coordi, *beam, peak = 0., tmp1, tmp2;
    SACHEAD hd;
    char *ss, *slow = {"slow"}, *backazimuth = {"baz"}, ch[20];
    FILE *fin, *fout, *fbp, *fp;

    if ( argc != 13 ) {
        fprintf(stderr,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        fprintf(stderr,"++    USAGE: vespa <filelist> <t1> <t2> <f1> <f2> <ID:slow/baz> <slow/baz> <slow_low/baz_low> <slow_high/baz_high> <step_len> <Nth_root> <output>    ++\n");
        fprintf(stderr,"++            return output_results of 3 columns: col1: time; col2: slowness/backazimuth col3: amplitude                                             ++\n");
        fprintf(stderr,"++            <filelist>           [1]  file containing every SAC format file                                                                        ++\n");
        fprintf(stderr,"++            <t1>                 [2]  beginning time of executing Velocity Spectral Analisys                                                       ++\n");
        fprintf(stderr,"++            <t2>                 [3]  ebding time of executing Velocity Spectal Analysis                                                           ++\n");
        fprintf(stderr,"++            <f1>                 [4]  low limitation of corner frequency of bandpass filter                                                        ++\n");
        fprintf(stderr,"++            <f2>                 [5]  high limitation of corner frequency of bandpass filter                                                       ++\n");
        fprintf(stderr,"++            <ID:slow/baz>        [6]  \"slow\" or \"baz\", which means fix slowness/baz and scan baz/slowness                                          ++\n");
        fprintf(stderr,"++            <slow/baz>           [7]  fixed slowness/backazimuth                                                                                   ++\n");
        fprintf(stderr,"++            <slow_low/baz_low>   [8]  low limitation of scanning sloness or backazimuth                                                            ++\n");
        fprintf(stderr,"++            <slow_high/baz_high> [9]  high limitation of scanning slowness or backazimuth                                                          ++\n");
        fprintf(stderr,"++            <step_len>           [10] step length of scanning slowness/backazimuth                                                                 ++\n");
        fprintf(stderr,"++            <Nth_root>           [11] Nth root slant-stacking, here specially, N consistant with 1 means linear stacking                           ++\n");
        fprintf(stderr,"++            <output>             [12] file saving outputting results: col1: time col2: slowness/baz col3: amplitude                                ++\n");
        fprintf(stderr,"++                                 ATTENTION!!! Executing VESPA requires SAC(Seismic Analysis Code) to filter with bandpass                          ++\n");
        fprintf(stderr,"++                                              Executing plot.sh requires GMT(the Generic Mapping Tools) with primary version 5,here 5.3.1          ++\n");
        fprintf(stderr,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        exit(1);
    }

    fin = fopen(argv[1],"r"); fout = fopen(argv[12],"w"); fp = fopen("plot.sh","w");
    t1 = atof(argv[2]); t2 = atof(argv[3]); f1 = atof(argv[4]); f2 = atof(argv[5]); ss = (char*) malloc(size); Nth_root = atof(argv[11]);

    time_start = clock();
    while ( fgets(ss, size, fin) ) {
        no_spa(ss);
        data = read_sac(ss, &hd);
        delta = hd.delta; sac_npts = hd.npts; center_lon += hd.stlo; center_lat += hd.stla; count += 1;
    }
    free(data); fclose(fin);
    center_lon /= count; center_lat /= count; begin_index = (int) (t1/delta); end_index = (int) (t2/delta); beam_npts = (int) ((t2-t1)/delta);

/*--------------------------------------------------------allocate dynamic memories to array "allamp", "coordi" and "amp"---------------------------*/
    allamp = (float**) malloc(sizeof(float*) * count); coordi = (float**) malloc(sizeof(float*) * count);
    beam = (float*) malloc(sizeof(float) * sac_npts);
    for ( i = 0; i < count; i ++ ) {
        allamp[i] = (float*) malloc(sizeof(float) * sac_npts);
        coordi[i] = (float*) malloc(sizeof(float) * 2);
    }

/*--------------------------------------------------------filtering ans save data in corresponding arrays-------------------------------------------*/
    fin = fopen(argv[1], "r");
    while ( fgets( ss, size, fin ) ) {
        no_spa(ss);
        fbp = fopen("bandpassfilter.sh","w");
        fprintf(fbp,"SAC_DISPLAY_COPYRIGHT=0\n");
        fprintf(fbp,"sac<<END\n");
        fprintf(fbp,"r %s\n", ss);
        fprintf(fbp,"bp c %f %f n 4 p 2\n", f1, f2);
        strcat(ss,".bp");
        fprintf(fbp,"w over %s\n", ss);
        fprintf(fbp,"q\n");
        fprintf(fbp,"END\n");
        fclose(fbp); system("sh bandpassfilter.sh"); system("rm bandpassfilter.sh");
        data = read_sac(ss, &hd);
        coordi[sta_index][0] = hd.stlo; coordi[sta_index][1] = hd.stla;
        for ( i = 0; i < sac_npts; i ++ ) allamp[sta_index][i] = data[i];
        sta_index += 1;
    }
    system("rm *.bp"); free(data); fclose(fin);

/*--------------------------------------------------------fixed slowness and scanning backazimuth---------------------------------------------------*/
    if ( !strcmp(slow, argv[6]) ) {
        float slowness, baz_low, baz_high, baz_step, baz_scan, slow_x, slow_y;
        slowness = atof(argv[7]); baz_low = atof(argv[8]); baz_high = atof(argv[9]); baz_step = atof(argv[10]);
        baz_scan = baz_low;
        while ( baz_scan <= baz_high ) {
            for ( i = 0; i < beam_npts; i ++ ) beam[i] = 0.;
            for ( i = 0; i < count ; i++ ) {
                slow_x = slowness * cos( (90.-baz_scan)/180.*PI ); slow_y = slowness * sin( (90.-baz_scan)/180.*PI );
                dx = center_lon - coordi[i][0]; dy = center_lat - coordi[i][1];
                shift_time = slow_x * dx + slow_y * dy; shift_index = (int) (shift_time/delta);
                for ( j = 0; j < beam_npts; j ++ ) {
                    if ( (begin_index + j + shift_index) >= 0 && (begin_index + j + shift_index) < sac_npts ) {
                        tmp1 = sign(allamp[i][begin_index+j+shift_index])*pow(fabs(allamp[i][begin_index+j+shift_index]),1./Nth_root)/count;
                        beam[j] += tmp1;
                    }
                    else beam[j] += 0.;
                }
            }
            for ( j = 0; j < beam_npts; j ++ ) {
                time = t1 + j * delta;
                tmp2 = sign(beam[j])*pow(fabs(beam[j]),Nth_root);
                beam[j] = tmp2;
                if ( peak < fabs(beam[j]) ) peak = fabs(beam[j]);
                fprintf(fout, "%f %f %f\n", time, baz_scan, beam[j]);
            }
            baz_scan += baz_step;
        }
    /*----------------------------------------------------Scanning backazimuth:save plot script in file "plot.sh"----------------------- -----------*/
    fprintf(fp,"R=%f/%f/%f/%f\n", t1, t2, baz_low, baz_high);
    fprintf(fp,"J=X9i/6i\nPS=%f~%f.ps\nPDF=%f~%f.pdf\n", t1, t2, t1, t2);
    fprintf(fp,"gmt gmtset FONT_LABEL 25,27,black\n");
    fprintf(fp,"gmt gmtset FONT_TITLE 30,27,black\n");
    fprintf(fp,"awk '{print $1,$2,$3/%f}' %s >tmp.txt\n", peak, argv[12]);
    fprintf(fp,"gmt surface tmp.txt -R$R -I%f/%f -G%s.grd\n", delta*10, baz_step/2, argv[12]);
    fprintf(fp,"gmt makecpt -Cpolar.cpt -T-1/1.0/0.01 -Z>tmp.cpt\n");
    fprintf(fp,"gmt psxy -R$R -J$J -K -T>$PS\n");
    fprintf(fp,"gmt grdimage %s.grd -R -J -K -O -Bx%d+l\"Time(sec)\" -By%f+l\"Backazimuth(deg)\" -BWSen+t\"Fixed slowness sec/deg%.3f\" -Ctmp.cpt>>$PS\n", argv[12], (int)((t2-t1)/10), (baz_high-baz_low)/10.,slowness);
    fprintf(fp,"gmt psscale -Ctmp.cpt -R -J -K -O -D9.4i/3i/15.1/1 -Bx0.2+l\"Normalized Spectral Power\">>$PS\n");
    fprintf(fp,"gmt psxy -R -J -O -T>>$PS\nps2pdf $PS $PDF\n");
    fprintf(fp,"gmt psconvert -Tg -A -P $PS\n");
    fprintf(fp,"rm tmp.* gmt.*\n"); fprintf(fp,"evince $PDF\n");
    fclose(fp);

    }

/*--------------------------------------------------------fixed backazimuth and scanning slowness---------------------------------------------------*/
    else if ( !strcmp(backazimuth, argv[6]) ) {
        float backazimuth, slow_low, slow_high, slow_step, slow_scan, slow_x, slow_y;
        backazimuth = atof(argv[7]); slow_low = atof(argv[8]); slow_high = atof(argv[9]); slow_step = atof(argv[10]);
        slow_scan = slow_low;
        while ( slow_scan <= slow_high ) {
            for ( i = 0; i < beam_npts; i ++ ) beam[i] = 0.;
            for ( i = 0; i < count; i ++ ) {
                slow_x = slow_scan * cos( (90.-backazimuth)/180.*PI ); slow_y = slow_scan * sin( (90.-backazimuth)/180.*PI );
                dx = center_lon - coordi[i][0]; dy = center_lat - coordi[i][1];
                shift_time = slow_x * dx + slow_y * dy; shift_index = (int) (shift_time/delta);
                for ( j = 0; j < beam_npts; j ++ ) {
                    if ( (begin_index + j + shift_index >= 0) && (begin_index + j + shift_index < sac_npts) ) {
                        tmp1 = sign(allamp[i][begin_index+j+shift_index])*pow(fabs(allamp[i][begin_index+j+shift_index]),1./Nth_root)/count;
                        beam[j] += tmp1;
                    }
                    else beam[j] += 0.;
                }
            }
            for ( j = 0; j < beam_npts; j ++ ) {
                time = t1 + j * delta;
                tmp2 = sign(beam[j])*pow(fabs(beam[j]), Nth_root);
                beam[j] = tmp2;
                if ( peak < fabs(beam[j]) ) peak = fabs(beam[j]);
                fprintf(fout,"%f %f %f\n", time, slow_scan, beam[j]);
            }
            slow_scan += slow_step;
        }
    /*----------------------------------------------------Scanning slowness:save plot script in file "plot.sh"--------------------------------------*/
    fprintf(fp,"R=%f/%f/%f/%f\n", t1, t2, slow_low, slow_high);
    fprintf(fp,"J=X9i/6i\nPS=%f~%f.ps\nPDF=%f~%f.pdf\n", t1, t2, t1, t2);
    frpintf(fp,"gmt gmtset FONT_LABEL 25,27,black\n");
    fprintf(fp,"gmt gmtset FONT_TITLE 30,27,black\n");
    fprintf(fp,"awk '{print $1,$2,$3/%f}' %s >tmp.txt\n", peak, argv[12]);
    fprintf(fp,"gmt surface tmp.txt -R$R -I%f/%f -G%s.grd\n", delta*10, slow_step/2, argv[12]);
    fprintf(fp,"gmt makecpt -Cpolar.cpt -T-1/1./0.01 -Z>tmp.cpt\n");
    fprintf(fp,"gmt psxy -R$R -J$J -K -T>$PS\n");
    fprintf(fp,"gmt grdimage %s.grd -R -J -K -O -Bx%d+l\"Time(sec)\" -By%f+l\"Slowness(sec/deg)\" -BWSen+t\"Fixed back-azimuth %.2f deg\" -Ctmp.cpt>>$PS\n", argv[12], (int)((t2-t1)/10), (slow_high-slow_low)/10., backazimuth);
    fprintf(fp,"gmt psscale -Ctmp.cpt -R -J -K -O -D9.4i/3i/15.1/1 -Bx0.2+l\"Normalized Spectral Power\">>$PS\n");
    fprintf(fp,"gmt psxy -R -J -O -T>>$PS\n"); fprintf(fp,"ps2pdf $PS $PDF\n");
    fprintf(fp,"gmt psconvert -Tg -A -P $PS\n");
    fprintf(fp,"rm tmp.* gmt.*\n"); fprintf(fp,"evince $PDF\n");
    fclose(fp);
    }

/*--------------------------------------------------------inputting ID (slow/baz) is not proper and exit--------------------------------------------*/
    else {
        fprintf(stderr,"Inputting unproper ID and now exiting ...\n");
        exit(1);
    }

/*--------------------------------------------------------output array coordinates------------------------------------------------------------------*/
    for ( i = 0; i < count; i ++ ) {
        printf("  Station %3d longitude ---> %10.6f  latitude ---> %10.6f\n", i+1, coordi[i][0], coordi[i][1]);;
    }
    printf("Array center: longitude ---> %10.6f  latitude ---> %10.6f\n", center_lon, center_lat);

 /*--------------------------------------------------------release dynamic memories------------------------------------------------------------------*/
    for ( i = 0; i < count; i++ ) {
        free(allamp[i]); free(coordi[i]);
    }
    free(allamp); free(coordi); free(beam); fclose(fout);



    time_end = clock();
    printf("Time used: %.6f seconds!!!\n", (time_end-time_start)/CLOCKS_PER_SEC);
    return 0;
}
