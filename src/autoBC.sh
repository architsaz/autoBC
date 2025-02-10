#!/bin/bash

# Get the name of aneu case from argument
casename=$1

#finding mdl directory in aneu database 
mdl_dir=$(find /dagon1/jcebral/aneuDB/db.R01/db.*/ -maxdepth 1 -name $casename)/mdl

#cp zhemo.info file:
   if [ -e "$mdl_dir/zhemo.info" ]; then
      cp $mdl_dir/zhemo.info  .
   else
      cp $mdl_dir/zhemo.info.gz  .
      gzip -df zhemo.info.gz 
   fi
   if [ ! -e zhemo.info ]; then
      echo "ERROR: zhemo.info file is not availabe for $casename" 
      exit 1 # Failure signal
   fi
#cp Open surface file :
   surf=$(head -n 1 zhemo.info | sed "s/INPUT: //")
   if [ -e $mdl_dir/$surf ]; then
      cp $mdl_dir/$surf  . 2>/dev/null
   else
      if [ -e mdl/$surf.gz ]; then
          cp $mdl_dir/$surf.gz  . 2>/dev/null
      else
          cp $(find /dagon1/jcebral/aneuDB/db.R01/db.*/ -maxdepth 1 -name $casename)/bld/$surf  . 2>/dev/null
          cp $(find /dagon1/jcebral/aneuDB/db.R01/db.*/ -maxdepth 1 -name $casename)/bld/$surf.gz  . 2>/dev/null
          modi_surf=$(basename $surf .zfem)
          cp $(find /dagon1/jcebral/aneuDB/db.R01/db.*/ -maxdepth 1 -name $casename)/bld/$modi_surf.gz  ./$surf.gz 2>/dev/null
      fi
      gzip -df $surf.gz 
   fi

   if [ ! -e $surf ]; then
      echo "ERROR: $surf file is not availabe for $casename" 
      exit 1 # Failure signal 
   fi
# generate the corrdinate for centeriod of inlet/outlet 
    rm -r inlet.txt outlet.txt >>/dev/null
    zhemo_cmdl --gen-centers >>/dev/null
    if [ ! -e centers.txt ]; then
      echo "ERROR: centers.txt file is not availabe for $casename" 
      exit 1 # Failure signal 
   fi
    mapfile -t inlet_id < <(grep "INFLOW" zhemo.info | awk '{print $2}')
    mapfile -t outlet_id < <(grep "OUTFLOW" zhemo.info | awk '{print $2}')
    FILENAME="centers.txt"
    line_count=0
    while read -r line; do
        for id in "${inlet_id[@]}"; do
            if [[ "$id" == "$line_count" ]]; then
                #echo "$id inlet: $line"
                echo "$line" >> inlet.txt
                break
            fi
        done
        for id in "${outlet_id[@]}"; do
            if [[ "$id" == "$line_count" ]]; then
                #echo "$id outlet; $line"  
                echo "$line" >> outlet.txt
                break
            fi
        done
        ((line_count++))  # Increment line count
    done < "$FILENAME"
    if [ ! -e inlet.txt ] || [ ! -e outlet.txt ]; then
        echo "ERROR: inlet or outlet files did not generated for $casename"
        exit 1  # Failure signal
    fi
# creat mask file
/dagon1/achitsaz/mylib/EXECs/EXEC_autoBC $casename
