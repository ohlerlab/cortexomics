set -uev

files=''

pattern=${1:-''}

if [ "$pattern" ]; then 
	echo "filtering for files matching $pattern"
	files=$(snakemake --summary  | cut -f 1 | tail -n+2 |   grep -Pe "$pattern" | xargs echo ) 
	
	if [ -z "$files" ]; then
		echo "no files matching"
		exit
	fi
fi

mkdir -p cluster_log

set -x 
snakemake -j 50 -k -p --restart-times 0 --max-jobs-per-second 5 -s Snakefile --cluster-config ../src/config_pipeline.json --rerun-incomplete --use-conda '--drmaa=  --mem={cluster.mem}  --time={cluster.time}   -p {cluster.A} -e {cluster.log_e} -o {cluster.log_o}' $files
