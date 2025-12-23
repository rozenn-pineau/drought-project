#!/bin/bash
  
#goal : to separate the vcf file into chromosomes for faster coding

for chr in {1..16}; do
        echo "processing: $chr "
        awk -v awk_var="$chr" ' $1== awk_var ' two_pulse_flexible_prop_2_values_ID.vcf > ${chr}.vcf
done
