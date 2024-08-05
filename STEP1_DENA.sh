#!/bin/bash
echo Start time is `date +%Y/%m/%d--%H:%M`

work_dir=/data1/m6A_analysis/UC/Patient
cd ${work_dir}
mkdir DENA_N
mkdir ${work_dir}/DENA_N/single_read_fast5
mkdir ${work_dir}/DENA_N/basecalls
sample_name=Patient_N

/data1/miniconda2/envs/DENA/bin/python /data1/DENA-release/step4_predict/LSTM_extract.py get_pos --fasta /data1/Reference/Homo_sapiens_assembly38.fasta  --motif 'RRACH' --output ${work_dir}/DENA_N/candidate_predict_pos.txt

/data1/miniconda2/envs/DENA/bin/multi_to_single_fast5 -t 10 -i /data1/m6A_analysis/UC/Patient/pangguang_N/fast5_pass -s ${work_dir}/DENA_N/single_read_fast5 --recursive

/data1/ont-guppy-cpu/bin/guppy_basecaller -i ${work_dir}/DENA_N/single_read_fast5 -s ${work_dir}/DENA_N/basecalls --flowcell FLO-MIN106 --kit SQK-RNA002 --cpu_threads_per_caller 5 --qscore_filtering --fast5_out --records_per_fastq 0 --recursive

cat ${work_dir}/DENA_N/basecalls/pass/*.fastq > ${work_dir}/DENA_N/basecalls/basecalls.fastq

/data1/miniconda2/envs/DENA/bin/tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir ${work_dir}/DENA_N/single_read_fast5 --fastq-filenames ${work_dir}/DENA_N/basecalls/basecalls.fastq --processes 4 --sequencing-summary-filenames ${work_dir}/DENA_N/basecalls/sequencing_summary.txt

/data1/miniconda2/envs/DENA/bin/tombo resquiggle --rna --processes 4 --corrected-group RawGenomeCorrected_001 --basecall-group Basecall_1D_000 --include-event-stdev --overwrite --ignore-read-locks ${work_dir}/DENA_N/single_read_fast5 /data1/Reference/Homo_sapiens_assembly38.fasta

/data1/miniconda2/envs/py39/bin/minimap2 -ax map-ont -L --secondary=no /data1/Reference/Homo_sapiens_assembly38.fasta ${work_dir}/DENA_N/basecalls/basecalls.fastq | samtools view -bh -F 2324 | samtools sort -O bam > ${work_dir}/DENA_N/basecalls.bam
samtools index ${work_dir}/DENA_N/basecalls.bam

cd ${work_dir}/DENA_N/

/data1/miniconda2/envs/DENA/bin/python /data1/DENA-release/step4_predict/LSTM_extract.py predict --processes 4 --fast5 ${work_dir}/DENA_N/single_read_fast5  --corr_grp RawGenomeCorrected_001 --bam ${work_dir}/DENA_N/basecalls.bam --sites ${work_dir}/DENA_N/candidate_predict_pos.txt --label ${sample_name} --windows 2 2

mkdir ${work_dir}/DENA_N/${sample_name}

## Predict
/data1/miniconda2/envs/DENA/bin/python /data1/DENA-release/step4_predict/LSTM_predict.py -i ${work_dir}/DENA_N -m /data1/DENA-release/DENA_LSTM_Model -o ${work_dir}/DENA_N/${sample_name} -p ${sample_name} -d

echo Finish time is `date +%Y/%m/%d--%H:%M`

