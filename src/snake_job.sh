

files=$(snakemake -n --rerun-incomplete | grep -e 'rule all' -A1 | tail -n1 | sed 's/,/\n/g' | grep -v 'input'  |  grep -Pe "$1" | xargs echo )
if $files ; then
	echo "no files matching"
	exit
fi
echo 'runnning'
echo $files

snakemake -j 50  -k -p --restart-times 1 --max-jobs-per-second 5 -s Snakefile --cluster-config ../src/config_pipeline.json  --rerun-incomplete --use-conda --drmaa=" -cwd -V -l h_vmem={cluster.h_vmem} -l h_rt={cluster.h_rt} -pe {cluster.pe} -j yes -P medium -o sge_log" $files
