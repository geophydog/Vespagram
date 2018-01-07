:hotel:[Return to Home Page](https://github.com/geophydog/geophydog.github.io/blob/master/README.md)

***

# Vespagram
- To execute velocity spectrum analysis (VESPAGRAM or VESPA)

***

# Basic principles
_Rost S, Thomas C (2002) Array seismology: methods and applications. Rev Geophys 40:1008. doi:
10.1029/2000RG000100_
# Dependencies
- Linux or Mac OS platform
- SAC(Seismic Analysis Code, here, version 101.6a)
    - [Request for SAC](http://ds.iris.edu/ds/nodes/dmc/forms/sac/)
- GMT(The Generic Mapping Tools, here version 5.3.1)
    - [Download GMT](http://gmt.soest.hawaii.edu/projects/gmt/wiki/Download)

***

# Installation
- After downloading the C programming source code just run "make" to compile it, and you will get the executable binary file "vespa", and copy it in your bin directory.

***

# Usage
vespa sacfile.lst t1 t2 fre_low fre_high ID<baz/slow> baz/slow slow_low/baz_low slow_high/baz_high slow_step/baz_step Nth_root output_file

- NOTICE!!!
 ```
 Here, slowness actually means horizontal slowness or ray parameter.
 ```

|   parameter   | mean  |
| ------------- | ----- |
| sacfile.lst   | file containing names of SAC format files |
|      t1       | begin time of doing vespagram   |
|      t2       | end time of doing vespagram     |
|      fre_low  | low limitation of corner frequency of SAC files |
|      fre_high | high limitation of corner frequency of SAC files |
|  ID<baz/slow> | "baz": fixed back-azimuth and scanning slowness; "slow": fixed slowness and scanning back-azimuth  |
|    baz/slow   | if back-azimuth fixed, its specific fixed value; or slowness fixed, its specific fixed value  |
| slow_low/baz_low | if back-azimuth fixed, low limitation of scanning slowness; or slowness fixed, low limitation of scanning back-azimuth|
| slow_high/baz_high| if back-azimuth fixed, high limitation of scanning slowness; or slowness fixed, high limitation of scanning back-azimuth |
| slow_step/baz_step | if back-azimuth fixed, step length of scanning slowness; or slowness foxed, step length of scanning back-azimuth |
| Nth_root | stacking method, specially, it means linear stacking when N is 1 |
| output_file | file name of outputting results |

***

# Example
### :one: The seismic stations mapping of this array.
   ![STATIONS](https://github.com/geophydog/Vespagram/blob/master/images/STA.jpg)
    
### :two: SAC format data recorded by this seismic array.
   ![SAC](https://github.com/geophydog/Vespagram/blob/master/images/SAC.jpg)
    
### :three: Find the back-azimuth with beam-forming method
   - [Beam-forming](https://github.com/geophydog/Beamforming_in_time_domain)
   ![BAZ](https://github.com/geophydog/Vespagram/blob/master/images/BEAM-FORMING.jpg)
  _Its back-azimuth is about 116 degrees_

### :four: Fix the back-azimuth with 116 degrees and scan horizontal slowness or ray parameters.
   ![RESULTS](https://github.com/geophydog/Vespagram/blob/master/images/results.jpg)

***

# Contribution
- Athor: Xuping Feng
- Email: geophydogvon@gmail.com.
