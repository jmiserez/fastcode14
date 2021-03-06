#!/bin/bash
echo "======================================================="
echo " Runtime"
echo "======================================================="
z=$((0)); for i in *.run.log; do z=$(($z+1)); (echo "#$i"); echo "yvals$z = ("; x=$(cat $i | grep "CPU_CYCLE" | tail -n6 | tr -d " " | cut -d ":" -f8| sed 's/$/,/'); y=$(echo "$x" | tr '\n' ' ');  echo "$y)"; done
echo "======================================================="
echo " Performance"
echo "======================================================="
z=$((0)); for i in *.run.log; do z=$(($z+1)); (echo "#$i"); echo "yvals$z = ("; x=$(cat $i | grep "Performance  " | tail -n6 | tr -d " " | cut -d ":" -f2| sed 's/$/,/'); y=$(echo "$x" | tr '\n' ' ');  echo "$y)"; done
echo "======================================================="
echo " Flops"
echo "======================================================="
z=$((0)); for i in *.run.log; do z=$(($z+1)); (echo "#$i"); echo "yvals$z = ("; x=$(cat $i | grep "Flops" | tail -n6 | tr -d " " | cut -d ":" -f2| sed 's/$/,/'); y=$(echo "$x" | tr '\n' ' ');  echo "$y)"; done
echo "======================================================="
echo " Operational intensity"
echo "======================================================="
z=$((0)); for i in *.run.log; do z=$(($z+1)); (echo "#$i"); echo "yvals$z = ("; x=$(cat $i | grep "Operational" | tail -n6 | tr -d " " | cut -d ":" -f2 | sed 's/$/,/' ); y=$(echo "$x" | tr '\n' ' ');  echo "$y)"; done
echo "======================================================="
echo " Memory bandwidth in bytes/cycle"
echo "======================================================="
z=$((0)); for i in *.run.log; do z=$(($z+1)); (echo "#$i"); echo "yvals$z = ("; x=$(cat $i | grep "bytes" | tail -n6 | tr -d " " | cut -d ":" -f2 | sed 's/$/,/' ); y=$(echo "$x" | tr '\n' ' ');  echo "$y)"; done
echo "======================================================="
echo " Memory bandwidth in MB/s"
echo "======================================================="
z=$((0)); for i in *.run.log; do z=$(($z+1)); (echo "#$i"); echo "yvals$z = ("; x=$(cat $i | grep "MiB" | tail -n6 | tr -d " " | cut -d ":" -f2 | sed 's/$/,/' ); y=$(echo "$x" | tr '\n' ' ');  echo "$y)"; done
echo "======================================================="
echo " Cache miss rate"
echo "======================================================="
z=$((0)); for i in *.run.log; do z=$(($z+1)); (echo "#$i"); echo "yvals$z = ("; x=$(cat $i | grep "Rate" | tail -n6 | tr -d " " | cut -d ":" -f2| sed 's/$/,/'); y=$(echo "$x" | tr '\n' ' ');  echo "$y)"; done

