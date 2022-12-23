for i in {1..421}
do
sed s/NUMBER1/${i}/ /projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/rsetlist_template.r | sed s/NUMBER1/${i}/  > /projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/files/makersetlist${i}
sed s/INDEXNUMBER/${i}/ /projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/exe_temp > /projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/files/exeR_${i}
qsub /projects/lu-lab/BenStuff/research/projects/oldhome/scripts/networks/rset/rsetlists/files/exeR_${i}
done
