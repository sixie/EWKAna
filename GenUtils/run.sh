RESBOSDIR=/data/blue/ksung/ResBos/resbos_cp/output

#
# d-dbar process
#
DATASETS=(
zee_ddb_resumm_ct10
zee_ddb_resumm_nlo_ct10
zee_ddb_resumm_ct66
zee_ddb_resumm_nlo_ct66
zee_ddb_resumm_nnlo_ct66
)

for i in "${DATASETS[@]}"
do
echo "Processing ${i}..."
XSEC=`grep "TOTAL CROSS SECTION" $RESBOSDIR/$i.out | awk '{print $5}'`
NEVENTS=`grep -c "^ 90" $RESBOSDIR/$i.hep`
root -l -q -b ResBosNtupler.C+\(\"${RESBOSDIR}/${i}.hep\",1,${XSEC},${NEVENTS},\"${RESBOSDIR}/${i}_ntuple.root\"\)
done

#
# u-ubar process
#
DATASETS=(
zee_uub_resumm_ct10
zee_uub_resumm_nlo_ct10
zee_uub_resumm_ct66
zee_uub_resumm_nlo_ct66
zee_uub_resumm_nnlo_ct66
)

for i in "${DATASETS[@]}"
do
echo "Processing ${i}..."
XSEC=`grep "TOTAL CROSS SECTION" $RESBOSDIR/$i.out | awk '{print $5}'`
NEVENTS=`grep -c "^ 90" $RESBOSDIR/$i.hep`
root -l -q -b ResBosNtupler.C+\(\"${RESBOSDIR}/${i}.hep\",2,${XSEC},${NEVENTS},\"${RESBOSDIR}/${i}_ntuple.root\"\)
done
