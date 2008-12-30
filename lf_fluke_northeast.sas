* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  This program uses MRFSS intercept survey data to generate          *
*  a representative length-frequency distribution for a particular    *
*  species.  The program does the following:                          *
*     1) outputs a table showing the distribution of measured fish    *
*        among state/mode/area/wave strata,                           *
*     2) outputs a table showing the distribution of estimated catch  *
*        among state/mode/area/wave strata,                           *
*     3) calculates and outputs relative length-class frequencies for *
*        each state/mode/area/wave stratum,                           *
*     4) calculates appropriate relative weighting factors to be      *
*        applied to the length-class frequencies for each             *
*        state/mode/area/wave stratum prior to pooling among strata,  *
*     5) produces unbiased estimates of length-class frequencies for  *
*        more than one state/mode/area/wave strata by summing         *
*        respectively weighted relative length-class frequencies      *
*        across strata.                                               *
*  The program also calculates relative length-class frequencies by   *
*  pooling the same data without properly weighting to account for    *
*  sampling biases toward particular strata.  This allows you to      *
*  directly compare the unbiased and biased estimates to assess       *
*  the impact of intended and/or unintended sampling biases on the    *
*  estimated size distribution of fish landed.                        *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *;

**********************************************
** CHANGE THE FOLLOWING DIRECTORY LOCATIONS **
**********************************************;

*** LOCATION OF INTERCEPT DATASETS ***;
libname agi "m:\intercept\ag";

*** LOCATION OF CATCH ESTIMATE DATASETS ***;
libname est "m:\products\estimates\catch\ag";
libname cmb "M:\programs\www_tools\comb_charter_est2000\combined_data";

*** LOCATION OF FORMATS LIBRARY ***;
libname library 'm:\programs\formats';

*** LOCATION OF DATASETS CREATED BY PROGRAM ***;
libname out "c:\mrfss\products\freq\trial";

options ls=80 ps=80;

*****************************************
** CHANGE THE FOLLOWING SPECIFICATIONS **
*****************************************;

*** SPECIES CODE ***;
%let species = '8835430101'; ** ScUp **;
%let label = SCUP;           ** for output length frequency file name **;

*** YEAR ***;
%let yr = 2007;

*** STATES TO BE INCLUDED IN POOLING OF LENGTH DATA ***;
%let states = (st=25);

*** CATCH ESTIMATE DATASETS TO BE INCLUDED IN ANALYSIS ***;
*%let estput=EST.AG_&yr.1 EST.AG_&yr.2 EST.AG_&yr.3 EST.AG_&yr.4 EST.AG_&yr.5 EST.AG_&yr.6;
%let estput=cmb.AG_&yr.1 cmb.AG_&yr.2 cmb.AG_&yr.3 cmb.AG_&yr.4 cmb.AG_&yr.5 cmb.AG_&yr.6;

*** TITLES FOR TABLES ***;
%let header = "Length Frequencies for Landed SCUP - MASSACHUSETTS";

******************************************************************************
* GO TO HIGHLIGHTED SECTIONS TO EDIT LENGTH CLASSES AND LENGTH CLASS FORMATS *
******************************************************************************;

run;

** DATA ON AVAILABLE OR "A" CATCH;

data lng_dat;
   set agi.i3_&yr.1 agi.i3_&yr.2 agi.i3_&yr.3
       agi.i3_&yr.4 agi.i3_&yr.5 agi.i3_&yr.6
       (keep=id_code sp_code sub_reg st
        mode_f mode_fx area_x wave year lngth);

   if year=2004 and sub_reg lt 6 then do;
      if mode_f='6' or mode_f='7' then mode_fx='6';
   end;

   if (year=2005 or year=2006) and sub_reg lt 6 then do;
      if mode_f='6' then mode_fx='4';
      if mode_f='7' then mode_fx='5';
   end;

* Subset species for which length frequencies are to be estimated;

   if sp_code = &species;

* Keep only fish records with length measurements;

   if lngth ne 9999 and lngth ne . ;

* Subset state/mode/area/wave cells to include in pooling of data;

   if &states;

* Convert lengths in millimeters to lengths in inches;

   length=(lngth/25.4);

*** Assign records to length classes for length frequency analyses;

*******************************
** EDIT LENGTH CLASSES BELOW **
*******************************;

   if length < 40 then lngcat = 40;
   if length < 39 then lngcat = 39;
   if length < 38 then lngcat = 38;
   if length < 37 then lngcat = 37;
   if length < 36 then lngcat = 36;
   if length < 35 then lngcat = 35;
   if length < 34 then lngcat = 34;
   if length < 33 then lngcat = 33;
   if length < 32 then lngcat = 32;
   if length < 31 then lngcat = 31;
   if length < 30 then lngcat = 30;
   if length < 29 then lngcat = 29;
   if length < 28 then lngcat = 28;
   if length < 27 then lngcat = 27;
   if length < 26 then lngcat = 26;
   if length < 25 then lngcat = 25;
   if length < 24 then lngcat = 24;
   if length < 23 then lngcat = 23;
   if length < 22 then lngcat = 22;
   if length < 21 then lngcat = 21;
   if length < 20 then lngcat = 20;
   if length < 19 then lngcat = 19;
   if length < 18 then lngcat = 18;
   if length < 17 then lngcat = 17;
   if length < 16 then lngcat = 16;
   if length < 15 then lngcat = 15;
   if length < 14 then lngcat = 14;
   if length < 13 then lngcat = 13;
   if length < 12 then lngcat = 12;
   if length < 11 then lngcat = 11;
   if length < 10 then lngcat = 10;
   if length <  9 then lngcat =  9;
   if length <  8 then lngcat =  8;
   if length <  7 then lngcat =  7;
   if length <  6 then lngcat =  6;
   if length <  5 then lngcat =  5;
   if length <  4 then lngcat =  4;
   if length <  3 then lngcat =  3;
   if length <  2 then lngcat =  2;
   if length <  1 then lngcat =  1;

* Format length category variable for later use in output tables;


*************************************
** EDIT LENGTH CLASS FORMATS BELOW **
*************************************;

proc format;
   value lngthcat
                  1 = 'Less than 1 inch'
                  2 = '1-1.99 inches'
                  3 = '2-2.99 inches'
                  4 = '3-3.99 inches'
                  5 = '4-4.99 inches'
                  6 = '5-5.99 inches'
                  7 = '6-6.99 inches'
                  8 = '7-7.99 inches'
                  9 = '8-8.99 inches'
                  10= '9-9.99 inches'
                  11= '10-10.99 inches'
                  12= '11-11.99 inches'
                  13= '12-12.99 inches'
                  14= '13-13.99 inches'
                  15= '14-14.99 inches'
                  16= '15-15.99 inches'
                  17= '16-16.99 inches'
                  18= '17-17.99 inches'
                  19= '18-18.99 inches'
                  20= '19-19.99 inches'
                  21= '20-20.99 inches'
                  22= '21-21.99 inches'
                  23= '22-22.99 inches'
                  24= '23-23.99 inches'
                  25= '24-24.99 inches'
                  26= '25-25.99 inches'
                  27= '26-26.99 inches'
                  28= '27-27.99 inches'
                  29= '28-28.99 inches'
                  30= '29-29.99 inches'
                  31= '30-30.99 inches'
                  32= '31-31.99 inches'
                  33= '32-32.99 inches'
                  34= '33-33.99 inches'
                  35= '34-34.99 inches'
                  36= '35-35.99 inches'
                  37= '36-36.99 inches'
                  38= '37-37.99 inches'
                  39= '38-38.99 inches'
                  40= '39-39.99 inches'
                 ;

* Calculate and save total measured fish by mode/area/wave cell;

proc freq data=lng_dat formchar='|----|+|---';
   title 'Distribution of measured fish among mode/area/wave strata';
   title2 &header;
   format st astatfor. mode_fx $tripfor. area_x $areafor. wave bimnth.;
   table st*mode_fx*area_x*wave / out=meas_fsh LIST;
run;

* Change name of count variable for total measured fish by m/a/w cell;

data meas_fsh; set meas_fsh;
   meas_fsh = count;
   drop count percent;
run;


* Obtain total estimated harvest by mode/area/wave cell;

data estland;

* Input Catch Estimate Files for Waves of Interest;

   set &estput. (keep=sp_code sub_reg st
       wave mode_fx area_x estclaim estharv);

* Designate species for which length frequencies are to be estimated;

   if sp_code = &species.;

* Define mode/area/wave cells to include for pooling of data;

   if &states;
   year=&yr;
   estland = estclaim + estharv;
run;

* Calculate total estimated harvest by mode/area/wave cell;

proc summary data=estland nway;
   class year st mode_fx area_x wave;
   var estland;
   output out=cell_lnd sum=cell_lnd;
   PROC PRINT;
    title "Landings of species";
run;


* Calculate total estimated harvest for all cells combined;

proc summary data=estland nway;
   class year;
   var estland;
   output out=all_land sum=tot_lnd;

run;

* Create file of fish weighting factors based on estimates of fish
   landed by mode/area/wave;

data cell_wgt;
   merge cell_lnd all_land; by year;
   wgt_fact = cell_lnd/tot_lnd;
proc sort; by st mode_fx area_x wave;
run;

* Calculate length frequencies by mode/area/wave;

proc freq data=lng_dat formchar='|----|+|---';
   title 'Length frequencies by mode/area/wave';
   title2 "&yr. - " &header;
   format st astatfor. mode_fx $tripfor. wave bimnth. area_x $areafor.
      lngcat lngthcat.;
   table st*mode_fx*area_x*lngcat*wave / out=lng_freq noprint;
run;

* Rename variable for total measured fish by
  mode/area/wave and length category;

data lng_freq;
   set lng_freq;
   lng_fish=count;
   sp_code=&species;
   drop count percent;
   proc sort;
    by st mode_fx area_x wave;
run;


* Merge length frequencies by mode/area/wave with weighting factors;

data lng_frq1;

   merge lng_freq meas_fsh cell_wgt; by st mode_fx area_x wave;

* Calculate weighted length class frequencies by mode/area/wave;
   if lng_fish ne .;
   lng_prop = (lng_fish/meas_fsh) * wgt_fact;

* Sum weighted frequencies to obtain representative
  length class frequencies;

proc summary data=lng_frq1 nway;
   class lngcat;
   id sp_code tot_lnd;
   var lng_prop;
   output out = lng_frq2 sum=tot_prop;

* Adjust length class proportions to account for cells with estimated
  catch but no measured lengths;

proc summary data=lng_frq2 nway;
   class sp_code;
   id tot_lnd;
   var tot_prop;
   output out = lng_frq3 sum=sum_prop;

data out.lfrq&label.;
   merge lng_frq2 lng_frq3 (keep=sp_code sum_prop);
      by sp_code;
   end_prop = tot_prop / sum_prop;
   end_lnd = end_prop*tot_lnd;
   end_pct = 100*end_prop;

/*
   proc dbload dbms=xls;
    path="C:\mrfss\products\freq\rdgrpr\rgprlfrq_&yr..xls";
    version=5;
    putnames=yes;
    limit=0;
    load;
   run;
*/

* Print list of relative frequencies of fish landed by length category
   for species of interest;

proc tabulate data=out.lfrq&label. formchar='|----|+|---';
   title 'Relative Frequencies of Landed Fish by Length Categories';
   title2 "&yr. - " &header;
   title3 'Data Reweighted prior to Pooling';
   class lngcat;
   var end_pct end_lnd;
   format lngcat lngthcat.;
   label end_pct='Estimated Relative Frequency (%)'
         end_lnd='Estimated Number of Landed Fish in Length Class';
   tables lngcat all='Total', end_pct*f=12.2 end_lnd*f=comma14.;
run;
