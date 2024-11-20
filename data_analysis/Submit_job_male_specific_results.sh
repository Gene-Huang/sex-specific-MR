#!/bin/bash

for exposure in OSA Insomnia Sleepiness ShortSleep LongSleep 
do 
  for bmi_type in bmi_unadj bmi_adj
  do
    for outcome in AF CAD CKD HF HTN T2DM
    do
      Rscript /Male_specific_results_code.R ${exposure} ${bmi_type} ${outcome}
    done
  done
done
          
          