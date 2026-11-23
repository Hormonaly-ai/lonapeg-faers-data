#!/bin/bash
# Download each quarterly zip (resumable, verified by zipfile open), reduce to parquet, delete zip. Logs to data/raw/download_log.txt.
# Usage: 00_pipeline.sh START END   (line numbers of data/raw/urls_2010Q4_2026Q3.txt)
T=$(cd "$(dirname "$0")/.." && pwd); OUT=$T/data/processed/faers_reduced; mkdir -p $OUT /tmp/tfaers
for i in $(seq $1 $2); do
  U=$(sed -n "${i}p" $T/data/raw/urls_2010Q4_2026Q3.txt); F=/tmp/tfaers/$(basename $U)
  Q=$(basename $U | grep -oiE '[0-9]{4}q[0-9]' | tr a-z A-Z)
  grep -q "\"quarter\": \"$Q\"" $OUT/process_log.jsonl 2>/dev/null && continue
  for try in 1 2 3 4 5 6; do
    S=$(python3 -c "import datetime;print(datetime.datetime.now(datetime.timezone.utc).strftime('%FT%TZ'))"); curl -s -A "Mozilla/5.0" -C - --retry 3 -o "$F" "$U"
    if python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).namelist()" "$F" 2>/dev/null; then break; fi
  done
  echo "$Q $(stat -c %s $F) start=$S end=$(python3 -c "import datetime;print(datetime.datetime.now(datetime.timezone.utc).strftime('%FT%TZ'))")" >> $T/data/raw/download_log.txt
  python3 $T/code/01_process_quarter.py "$F" $OUT >> /tmp/tfaers/proc_$1.log 2>&1 && rm -f "$F"
done
echo DONE $1 $2 >> /tmp/tfaers/proc_$1.log
