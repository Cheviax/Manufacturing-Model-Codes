#!/bin/sh

XPLAMMPSDIR="XPLAMMPS"
NBXP=16
INPUTSLURRYFILENAME="inliqreal.run"
INPUTVISCOSITYFILENAME="viscosity_slurry.run"



echo ${XPLAMMPSDIR}
rm -R ${XPLAMMPSDIR}*


for numxp in `seq 1 ${NBXP}`
do
    XPDIR=${XPLAMMPSDIR}${numxp}
    mkdir ${XPDIR}
    cp script.template ${XPDIR}/script.sh

    cp ${INPUTSLURRYFILENAME} ${XPDIR}/
    
    while [ ! -f ${XPDIR}/${INPUTSLURRYFILENAME} ]
    do
	echo "the file "${INPUTSLURRYFILENAME}" does not exit. Wait 1 second"
	sleep 1
    done

    cd ${XPDIR}
    sbatch --reservation midmem script.sh
    cd ..
done



# JE SAIS QUE X EST TERMINE

