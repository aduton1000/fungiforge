// Stage 04b — mitochondrial genome (W2.7). Annotates the mitochondrial contigs stage 04 separated:
// the 15 core genes by tblastn of the bundled reference proteins (presence, extra copies, introns
// from HSP gaps), rnl/rns by blastn, circularity from terminal repeats; then maps a read subsample
// to the mitogenome for depth (copy number relative to the nuclear genome) and heteroplasmy
// (minor alleles >= 10 % at >= 20x). Emits organelle.json and a simple GFF. An isolate without a
// separated mitogenome gets mito_present=false, never a failure.
process ORGANELLE {
  tag { meta.id }
  label 'decontam'
  publishDir { "${params.outdir}/${meta.id}/04_decontam" }, mode: 'copy', pattern: "*.{json,gff}"
  input:  tuple val(meta), path(mito), path(nuclear), path(ont), path(r1), path(r2)
  output: tuple val(meta), path("${meta.id}.organelle.json"), emit: json
          tuple val(meta), path("${meta.id}.mito.gff"),       emit: gff
  script:
  def has_illumina = r1.name != 'NO_R1'
  def sub_cmd = has_illumina ? "seqkit head -n ${params.organelle_reads} ${r1} > sub_1.fq && seqkit head -n ${params.organelle_reads} ${r2} > sub_2.fq"
                             : "seqkit head -n ${params.organelle_reads_ont} ${ont} > sub.fq"
  def reads   = has_illumina ? "sub_1.fq sub_2.fq" : "sub.fq"
  def preset  = has_illumina ? "sr" : "map-ont"
  def het_frac = has_illumina ? "0.10" : "0.20"       // ONT substitution noise needs a higher minor-allele floor
  """
  source ff_status.sh; ff_init organelle "${meta.id}" ${meta.id}.organelle.json --best-effort
  ff_version blastn -- blastn -version
  ff_version minimap2 -- minimap2 --version
  ff_version samtools -- samtools --version
  if [ ! -s ${mito} ]; then
    ff_skip annotate "no mitochondrial contig was separated from the assembly"
    mito_annotate.py --sample "${meta.id}" --mito ${mito} --out-json ${meta.id}.organelle.json --out-gff ${meta.id}.mito.gff
  else
    ff_run makeblastdb -- bash -c "makeblastdb -in ${mito} -dbtype nucl -out mito_db > makeblastdb.log 2>&1"
    ff_run tblastn -- bash -c "tblastn -query '${params.mito_proteins}' -db mito_db -evalue 1e-5 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen' > genes.tsv 2>tblastn.log"
    ff_run blastn_rrna -- bash -c "blastn -query '${params.mito_rrna}' -db mito_db -evalue 1e-10 -num_threads ${task.cpus} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' > rrna.tsv 2>blastn.log"
    MITO_DEPTH=""; NUC_DEPTH=""; PILEUP=""
    if [ "${params.organelle_reads_check}" = "true" ]; then
      ff_run subsample --optional -- bash -o pipefail -c "${sub_cmd}"
      if [ "\$FF_RC" -eq 0 ]; then
        cat ${mito} ${nuclear} > ref.fa
        ff_run map_reads --optional -- bash -o pipefail -c "minimap2 -ax ${preset} -t ${task.cpus} --secondary=no ref.fa ${reads} 2>minimap2.log | samtools sort -@ ${task.cpus} -m 1G -o reads.bam - && samtools index reads.bam"
        if [ "\$FF_RC" -eq 0 ]; then
          samtools coverage reads.bam > coverage.tsv
          MITO_DEPTH=\$(awk -v ids="\$(grep '>' ${mito} | sed 's/>//; s/ .*//' | tr '\\n' ' ')" 'BEGIN{n=split(ids,a," "); for(i=1;i<=n;i++) m[a[i]]=1} !/^#/ && (\$1 in m) {bp+=\$3-\$2+1; d+=(\$3-\$2+1)*\$7} END{if(bp>0) printf "%.1f", d/bp}' coverage.tsv)
          NUC_DEPTH=\$(awk -v ids="\$(grep '>' ${mito} | sed 's/>//; s/ .*//' | tr '\\n' ' ')" 'BEGIN{n=split(ids,a," "); for(i=1;i<=n;i++) m[a[i]]=1} !/^#/ && !(\$1 in m) {bp+=\$3-\$2+1; d+=(\$3-\$2+1)*\$7} END{if(bp>0) printf "%.1f", d/bp}' coverage.tsv)
          # pileup on the mitochondrial contigs from a BAM downsampled to ~300x; -B (no BAQ realignment:
          # on kilobase ONT reads BAQ exhausts memory and yields nothing); MAPQ 0 is kept because
          # overlapping fragments of the circle make every mitochondrial read a multi-mapper
          samtools faidx ${mito}
          FRAC=\$(python3 -c "d=float('\${MITO_DEPTH:-0}' or 0); print(min(1.0, 300.0/d) if d > 0 else 1.0)")
          samtools view -b -s "\$FRAC" -o mito_sub.bam reads.bam \$(cut -f1 ${mito}.fai | tr '\\n' ' ') && samtools index mito_sub.bam
          : > mito.pileup
          for c in \$(cut -f1 ${mito}.fai); do samtools mpileup -B -f ${mito} -q 0 -Q 15 -d 500 -r "\$c" mito_sub.bam 2>/dev/null >> mito.pileup || true; done
          PILEUP="--mpileup mito.pileup --het-min-frac ${het_frac}"
        fi
      fi
    else
      ff_skip map_reads "read check disabled (--organelle_reads_check false)"
    fi
    ff_run mito_annotate -- mito_annotate.py --sample "${meta.id}" --mito ${mito} --tblastn genes.tsv --rrna-blastn rrna.tsv \\
        \${PILEUP} \${MITO_DEPTH:+--mito-depth \$MITO_DEPTH} \${NUC_DEPTH:+--nuclear-depth \$NUC_DEPTH} \\
        --out-json ${meta.id}.organelle.json --out-gff ${meta.id}.mito.gff
  fi
  ff_finalize
  """
  stub:
  "echo '{\"sample\":\"${meta.id}\",\"stage\":\"organelle\",\"mito_present\":false}' > ${meta.id}.organelle.json; echo '##gff-version 3' > ${meta.id}.mito.gff"
}
