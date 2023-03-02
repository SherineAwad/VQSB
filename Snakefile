

rule all: 
      input: 
          expand("{cohort}_excesshet.vcf", cohort=config['COHORT']), 
          expand("{cohort}_siteonly.vcf", cohort=config['COHORT']),
          expand("{cohort}_indels.recal", cohort= config['COHORT']),
          expand("{cohort}_snps.recal", cohort= config['COHORT']),
          expand("{cohort}.indel.recalibrated.vcf.gz",cohort= config['COHORT']), 
          expand("{cohort}.snps.recalbrated.vcf.gz", cohort=config['COHORT'])
rule ExcessHet: 
        input: 
             "{cohort}.vcf"
        output: 
             "{cohort}_excesshet.vcf" 
        params: 
            config['ExcessHet']
	shell: 
           """ 
	   gatk --java-options "-Xmx3g -Xms3g" VariantFiltration -V {input} --filter-expression "ExcessHet > {params}" --filter-name ExcessHet -O {output}
	   """


rule MakeSitesOnlyVcf: 
     input: 
          "{cohort}_excesshet.vcf"
     output: 
          "{cohort}_siteonly.vcf" 
     shell: 
          """
           gatk MakeSitesOnlyVcf -I {input} -O {output}
          """

rule tranches_indels: 
     input: 
         "{cohort}_siteonly.vcf" 
     output: 
         "{cohort}_indels.recal", 
         "{cohort}_indels.tranches" 
     params: 
         gold_standard = config['GOLD_STANDARD'],
         axiom_exome = config['Axiom_Exome'], 
         dbsnp = config['DBSNP'],
     shell: 
        """
        gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
        -V {input} \
        --trust-all-polymorphic \
        -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
        -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \      
        -mode INDEL \
        --max-gaussians 4 \
        -resource:mills,known=false,training=true,truth=true,prior=12:{params.gold_standard}\
        -resource:axiomPoly,known=false,training=true,truth=false,prior=10:{params.axiom_exome} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2:{params.dbsnp} \
        -O {output[0]} 
        --tranches-file {output[1]}
       """ 

rule tranches_snps: 
     input:
          "{cohort}_siteonly.vcf" 
     output:
          "{cohort}_snps.recal", 
          "{cohort}_snps.tranches"
     params:
         hapmap_3 = config['hapmap_3'], 
         G_omni2 = config['G_omni2'],
         G_phase1 = config['G_phase1'],  
         dbsnp = config['DBSNP'],
     shell: 
          """ 
          gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
          -V {input} \
          --trust-all-polymorphic \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
          -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
          -mode SNP \
          --max-gaussians 6 \
          -resource:hapmap,known=false,training=true,truth=true,prior=15:{params.hapmap_3}\
          -resource:omni,known=false,training=true,truth=true,prior=12:{params.G_omni2} \
          -resource:1000G,known=false,training=true,truth=false,prior=10:{params.G_phase1} \
          -resource:dbsnp,known=true,training=false,truth=false,prior=7:{params.dbsnp} \
          -O {output[0]} \
          --tranches-file {output[1]} 
       """ 


rule apply_VQSR_indels: 
     input: 
          "{cohort}_excesshet.vcf",
          "{cohort}_indels.recal",
          "{cohort}_indels.tranches" #in previous rules this need to be output instead of param 
     output: 
          "{cohort}.indel.recalibrated.vcf.gz"
     params: 
     shell: 
          """
          gatk --java-options "-Xmx5g -Xms5g" \
   	  ApplyVQSR \
          -V {input[0]}\
          --recal-file {input[1]} \
          --tranches-file {input[2]} \
          --truth-sensitivity-filter-level 99.7 \
          --create-output-variant-index true \
          -mode INDEL \
          -O {output} 
          """ 

rule apply_VQSR_SNPs: 
     input: 
         "{cohort}.indel.recalibrated.vcf.gz", 
         "{cohort}_snps.recal",
         "{cohort}_snps.tranches", 
     output: 
         "{cohort}.snps.recalbrated.vcf.gz" 
     shell: 
         """ 
         gatk  --java-options "-Xmx5g -Xms5g" \
         -V {input[0]} \ 
         --recal-file {input[1]}  \
         --tranches-file {input[2]} \
         --truth-sensitivity-filter-level 99.7 \ 
         --create-output-variant-index true  \
         -mode SNP \
         -O {output} 
         """ 
 
