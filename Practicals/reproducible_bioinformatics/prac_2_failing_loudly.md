# Prac 2 Failing loudly - handling error in clinical screening

#### By Evelyn Collen


## **1. Introduction**

### What happens when we can't afford to fail - DPYD pharmacogenomic screening 


In last week's practical, we started to get the hang of handling failing scripts and the concept that failure, in general, isn't a bad thing! Here we are going to focus on how we can best handle errors in diagnostic reporting. 

The DPYD gene is responsible for generating the dihydropyrimidine dehydrogenase (DPD) enzyme, which plays a key role in the metabolism of toxic compounds. Deficiency in this enzyme can cause fatal toxicity to fluoropyrimidine chemotherapy treatments (e.g., 5-fluorouracil, capecitabine), which are widely used in the treatment of solid tumours such as colorectal cancer, breast cancer and gastrointestinal cancers. Variants in the DPYD can cause different functionality in the DPD enzyme's function, categorised into zero, decreased, or normal function. 

It's important to optimise the variant-tailored dosage, as the standard dose that is effective for some people can be fatal for certain variant carriers. Severe toxicity to these drugs occurs in about 10% to 40% of patients, and around 7% of Europeans carry some variant that impairs function. In South Australia, there has unfortunately been a recorded case of patient fatality due to DPD enzyme deficiency and an incorrectly tailored dosage. 


![alt text](images_and_refs/dpyd2.jpg)

Figure 1. The relative risk of toxicity in those with versus without the
specified variant when treated with full or individualised dose (taken from [Sonic Genetics](https://www.sonicgenetics.com.au/))
 
 
Depending on the configuration of alleles, and whether those alleles have zero, decreased or normal function, patients will receive a metabolism rating. Have a look at how the rating is worked out in this [DPYD metabolism rating table](images_and_refs/DPYD_metabolism_rating_and_recommendations.xlsx)


Today, we are going to determine whether four cancer patients are carrying particular variants in a DPYD gene, and what phenotype and metabolism rating can be deduced from their allele types. It's absolutely critical we get this diagnostically right in the end, as the oncologist will use this information to determine the best course of treatment that can aggressively handle the cancer, whilst minimising the effect on the patient as much as possible. 
The patient is waiting on this information - so not only we have to get it right, we want do it quickly. Last week we made sure our errors were noisy in our scripting; this week we are going to make our pipeline noisy about the whole process. 

### 1.1 Reminder about virtual Machines

As usual we will be connecting the virtual machines: 

**Please [go here](../../Course_materials/vm_login_instructions.md) for instructions on connecting to your VM.**

## 1.2 Learning Outcomes

1. Understand clinical failure in a clinical bioinformatics and how to minimise it
2. Learn troubleshooting errors from bioinformatic tools
3. Learn why QC at every crucial step is important to ensure errors are mitigated
4. Learn about how clinical screening results are curated and reported out 
5. Learn about variant reporting in the DPYD gene

## 1.3 About the dataset

Today we are looking at four patients, who I have anonymised to Patient A, Patient B, Patient C and Patient D. There is a bonus Patient E who has not been reported out, and you will soon find out why. 

Patients A,B,C and D have had real diagnostic reports issued out for them, and you can have a look at the anonymised versions here. [DPYD patient reports](images_and_refs/Patient_reports.docx)

**Questions:**
1.  Referring to the DPYD metabolism table, what score would be given to patient who has 1 normal funtion allele and 1 decreased function? Would you classify their phenotype as normal, intermediate or poor?
2. Referring to the DYPD patient reports, which patient actually has 2 decreased function alleles that are both classified as poor metabolisers? If you had to guess, what dosage would likely be given to this patient?
 
<details>
<summary>Answers</summary>
<ul><li>1. Intermediate metaboliser </li>
<ul><li>2. Probably a minimum dosage or an alternative therapy </li> </ul>
</details>


## **2. Errors with sample integrity**

### 2.1 Getting scripts and data ready

load software
```bash
source activate bioinf
```

create all directories (will do nothing if they already exist) and move into project directory
```bash
mkdir -p ~/Practical_Failing_Loudly/{0_scripts,1_vcfs,2_bam,3_reports,4_refs}
mkdir -p ~/Practical_Failing_Loudly/3_reports/{1_fingerprint_check,2_contam_check,3_varcall_check}
cd ~/Practical_Failing_Loudly
```

copy scripts and data and also make symlinks, if not existing from last week's prac

```bash
[ -f /dest/file ] || cp ~/data/failing_loudly/0_scripts/* 0_scripts/
[ -e /path/to/link ] || ln -s ~/data/failing_loudly/1_vcfs/*.vcf 1_vcfs/
[ -e /path/to/link ] || ln -s ~/data/failing_loudly/2_bam/*.bam 2_bam/
[ -e /path/to/link ] || ln -s ~/data/failing_loudly/4_refs/* 4_refs/
```

copy new scripts and data for this week's prac

```bash
ln -s ~/data/failing_loudly/1_vcfs/*.vcf 1_vcfs/
```

### 2.2 Checking a sample's integrity with genetic fingerprinting

The throughput of NGS samples going through clinical laboratories is really high, and getting higher every year. Some labs are seeing more than 20,000 NGS samples processed each year. This is pushing a lot of automation both in lab and in silico, but even now some lab steps are manual. With so many samples being manually handled at certain steps, how can we guarantee that no sample has been swapped or contaminated with another? 

One way is to separate the sample into two, right at the beginning when the lab first receives the sample. The first part of the sample goes through the normal testing process, and we generate data for it. The second part goes through a completely independent workflow, where we target just a handful of common SNPs.


```mermaid
flowchart TD
    A[Sample Received at Lab]

    A --> B[DNA extraction split into two]

    B --> E1[Primary Sample A undergoes sequencing for diagnostic variants]
    B --> E2[Counterpart Sample A undergoes sequencing for fingerprinting SNPs]

    E1 --> F[Check SNPs match with GATK check fingerprint]
    E2 --> F

    F --> G{Fingerprints Match?}

    G -->|Yes| H[Confirm Identity]
    G -->|No| I[Flag Discrepancy / Investigate]
```

This is just one of many sanity checks we can do, including pedigree check relatedness between family members, population ancestry, etc.

Let's start by running a fingerprint check for Patient A. The counterpart vcfs are in a separate folder, as you can see here, and have been produced completely indepedently from our diagnostic vcfs (so they won't have any variants in the DPYD gene). They contain mostly germline SNPs, that occur with faily high frequency (> 5%) in the most common population databases.

```bash
ls 1_vcfs/counterpart_vcfs/
```
Patient_A_counterpart.gatk.hg38.vcf	
Patient_B_counterpart.gatk.hg38.vcf
Patient_C_counterpart.gatk.hg38.vcf	
Patient_D_counterpart.gatk.hg38.vcf
Patient_E_counterpart.gatk.hg38.vcf


We're going to check the concordance of genotypes using GATK. We need two input ref files for the fingerprint command: 

```bash
ls 4_refs/HaplotypeMap.vcf
ls 4_refs/Homo_sapiens_assembly38.fasta
```

4_refs/HaplotypeMap.vcf is provided by the argument --HAPLOTYPE_MAP, while 4_refs/Homo_sapiens_assembly38.fasta is provided to the tool by the argument -R.

Now let's run the fingerprint command:

```bash
gatk CheckFingerprint -R 4_refs/Homo_sapiens_assembly38.fasta -I 1_vcfs/Patient_A.vcf.gz --GENOTYPES 1_vcfs/counterpart_vcfs/Patient_A_counterpart.gatk.hg38.vcf.gz --HAPLOTYPE_MAP --GENOTYPE_LOD_THRESHOLD 0 --SUMMARY_OUTPUT 3_reports/Patient_A.fingerprint_summary.tsv --DETAIL_OUTPUT 3_reports/Patient_A.fingerprint_detailMetrics.tsv
```

You should get this output:

```
Illegal argument value: Positional arguments were provided ',0}' but no positional argument is defined for this tool.
Tool returned:
1
```

Aha! An error! Lucky for us, we eat errors for breakfast. Notice, the GATK tool hasn't even said it's thrown an error - it's just listed all the required arguments (to nudge you to use the right ones) and told us at the bottom that we are doing something illegal. It's also told us the exit code is 1, which is developer speak for something is not right.

'Illegal argument value' tells us something might be wrong with our arguments or inputs. Can you figure out the issue now that you've had some experience digging into errors? Have a careful look at the inputs and give it a go. If you truly get stuck, the right command to run is hidden below. 

<details>
<summary>Fixed_command</summary>
<ul>```bash
gatk CheckFingerprint -R 4_refs/Homo_sapiens_assembly38.fasta -I 1_vcfs/Patient_A.vcf.gz --GENOTYPES 1_vcfs/counterpart_vcfs/Patient_A_counterpart.gatk.hg38.vcf.gz --HAPLOTYPE_MAP 4_refs/HaplotypeMap.vcf --GENOTYPE_LOD_THRESHOLD 0 --SUMMARY_OUTPUT 3_reports/1_fingerprint_check/Patient_A.fingerprint_summary.tsv --DETAIL_OUTPUT 3_reports/1_fingerprint_check/Patient_A.fingerprint_detailMetrics.tsv
``` </ul>
<details>

Well done! If the command worked, you should see in your ouput Patient_A: LOD = 19.548716624813856. 

You can also see this number if you look at one of the output files: 

```bash
cat 3_reports/1_fingerprint_check/Patient_A.fingerprint_summary.tsv
```

The LOD score, or LL_EXPECTED_SAMPLE (log-likelihood) in the output, is the core metric in this output. It represents the base-10 logarithm of the likelihood that, based on genotype similarity of the SNPs, the counterpart sample is an identical match to the primary sample, versus a random sample.
Positive Value: The counterpart sample matches the primary sample (e.g., a LOD of 6 means it is \(10^{6}\) or 1,000,000 times more likely to be a match than not).
Negative Value: The counterpart sample does not match the primary sample, indicating a potential swap or contamination.
Near Zero: Inconclusive result, usually due to low coverage or non-informative genotypes


### 2.2 Checking all samples for integrity

Let's now run this on all the samples with a simple for loop: 

```bash

for vcf in 1_vcfs/Patient_*.vcf.gz; do \
sample=$(basename $vcf .vcf.gz); \
gatk CheckFingerprint -R 4_refs/Homo_sapiens_assembly38.fasta -I ${vcf} --GENOTYPES 1_vcfs/counterpart_vcfs/${sample}_counterpart.gatk.hg38.vcf.gz --HAPLOTYPE_MAP 4_refs/HaplotypeMap.vcf --GENOTYPE_LOD_THRESHOLD 0 --SUMMARY_OUTPUT 3_reports/1_fingerprint_check/${sample}.fingerprint_summary.tsv --DETAIL_OUTPUT 3_reports/1_fingerprint_check/${sample}.fingerprint_detailMetrics.tsv; done

```
You can scroll through the output or have a look at the lod scores: 

```bash
cat 3_reports/1_fingerprint_check/Patient_*.fingerprint_summary.tsv
```

You're looking for the value for the LL_EXPECTED_SAMPLE column, which should be the number at the bottom, fourth from the left. 


**Questions:**
1. Have all patient samples passed the sample integrity check, based on LOD score?
2. What is the likelihood that Patient A's counterpart sample is a true match, versus that it was swapped with another patient's?
3. I haven't mentioned another really common sanity check to interrogate sample integrity. Can you guess what it is? 
4. What would happen to the LOD score of Patient A if the sample swap occured *prior* to the lab receiving the sample? Would it be negative?

<details>
<summary>Answers</summary>
<ul><li> 1. No, all of them passed except Patient E:
Patient_A 19.548717
Patient_B 18.665838
Patient_C 15.562386
Patient_D 21.008243
Patient_E -216.23926626114772 </li>
<li>2. Patient A's lod score is 19.548717, so the chance that it is a true match would be roughly \(10^{20}\) times more likely than a random swap  </li>
<li>3. Doing a sex check to see that the genetic sex matches expected patient sex </li>
<li>4. Nothing - the LOD score would still be positive, as the counterpart sample would have all the same genotypes as the main sample, seeing as the swap occured prior to the counterpart sample being split off. </li> </ul>
</details>


## **3. Errors with sample contamination**

### 3.1 Checking for any evidence of contamination

Now that we have verified that all the samples are correctly identified and that no swaps have occured, another thing we can check is whether there is any evidence of contamination, especially low-level contamination that could influence heterozygous calls. Remember that germline SNPs should usually only ever sit around 1 or 0.5 in frequency, for homozygous and heterozygous SNPs respectively? Well, if you see variants with any VAFs different to this, it could be because 1) it's a somatic variant that has accumulated over the patient's lifespan and not present in all cells, 2) there are reads from other patients contaminating the sample and making up the other portion of the allel fraction.

Since our counterpart vcfs mostly have SNPs that are common in population databases, they're more likely to be germline than somatic, so if we see widespread homozygous alt variants at these positions, it's more likele case 2). We can check this, by counting up the number of SNPs in the counterpart vcf that have a vaf less than 0.9, at homozygous alt positions.

```bash
python ./0_scripts/check_vaf.py 1_vcfs/counterpart_vcfs/Patient_A_counterpart.gatk.hg38.vcf.gz
```

The output should show you there is 1 such SNP for Patient A:

```
Warning: allele frequency for homozygous variant is 0.854, < 0.9
Sample_name	number_homozygous_SNPs_low_VAF	Contamination_status
Patient_A	1	OK
```
Since there's only 1, and it's not too far off from 0.9, we're not concerned about this and the sample passes. Genetic sequencing data can be noisy, and with real contamination we'd expect more sites to be affected. If there were perhaps more than 3 sites like this, we might start to worry. 

Let's run a for loop on the other patients, and store the outputs:

```bash
for vcf in 1_vcfs/counterpart_vcfs/Patient_*_counterpart.gatk.hg38.vcf.gz; \
do a=$(basename $vcf .gatk.hg38.vcf.gz); ./0_scripts/check_vaf.py $vcf 2>&1 | tee 3_reports/2_contam_check/${a}_contam_check; done

```
This script classifies sample contamination status as "OK" if the number of low vaf SNPs is < 3. 

**Question:**
1.  Do all the samples pass contamination check?

<details>
<summary>Answer</summary>
<ul><li>1. Yes! </li>
 </ul>
</details>



## **3. Errors with variant calling**

### 3.2 Run bcftools as a second variant caller check 

No variant caller is perfect, and even the best and more robust programs can make mistakes or have inconsistencies. The caller we used for the DPYD variants was Vardict, which is known to be a very good caller for amplicon data that almost never misses. But even it can have problems - if you're interested, take a quick glance at one of Vardict's issue pages on Github:

[text](https://github.com/AstraZeneca-NGS/VarDictJava/issues/81)

As I've been emphasising constantly, making the right call is absolutely crucial for the patients, and once we've established the sample is the right one, the next crucial step is to ascertain that the right call was made. We can go ahead and double-check it with the aid of a totally different variant caller: bcftools.

Note our inputs for this command:

```
./2_bam/Patient_A.bam
./4_refs/Homo_sapiens_assembly38.fasta
```


```bash
bcftools mpileup --count-orphans --no-BAQ \
  --max-depth 12345 \
  --min-MQ 10 \
  --skip-indels \
  --annotate AD \
  -f ./4_refs/Homo_sapiens_assembly38.fasta \
  -r chr1:97573863,chr1:97450058,chr1:97515787,chr1:97082391 \
   ./2_bam/Patient_AB.bam \
| bcftools call -c -a GQ \
| bcftools view -e 'GT="0/0"' \
| bcftools +fill-tags -O v -o 3_reports/3_varcall_check/Patient_A_bcftools_check.vcf -- -t FORMAT/VAF
```

Well well, do I smell another input error? Have a look at the error message from bcftools. 
Do you think you can work it out and fix the typo in the above command? 


<details>
<summary>Fixed_command</summary>
<ul>```bash
bcftools mpileup --count-orphans --no-BAQ \
  --max-depth 12345 --min-MQ 10 --skip-indels --annotate AD \
  -f ./4_refs/Homo_sapiens_assembly38.fasta \
  -r chr1:97573863,chr1:97450058,chr1:97515787,chr1:97082391 \
   ./2_bam/Patient_A.bam \
| bcftools call -c -a GQ | bcftools view -e 'GT="0/0"' \
| bcftools +fill-tags -O v -o 3_reports/3_varcall_check/Patient_A_bcftools_check.vcf -- -t FORMAT/VAF
``` </ul>
<details>


Great job! Run the following to compare the call made by Vardict to the one we just made with bcftools:

```bash
grep 97573863 1_vcfs/Patient_A.vcf
grep 97573863 3_reports/3_varcall_check/Patient_A_bcftools_check.vcf
```
Did they make the same call?

### 3.2 Run bcftools for all the samples


```bash
for i in ./2_bam/*.bam; do a=$(basename $i .bam); bcftools mpileup --count-orphans --no-BAQ \
  --max-depth 12345 --min-MQ 10 --skip-indels --annotate AD \
  -f ./4_refs/Homo_sapiens_assembly38.fasta \
  -r chr1:97573863,chr1:97450058,chr1:97515787,chr1:97082391 \
   ${i} \
| bcftools call -c -a GQ \
| bcftools view -e 'GT="0/0"' \
| bcftools +fill-tags -O v -o 3_reports/3_varcall_check/${a}_bcftools_check.vcf -- -t FORMAT/VAF; done
```

## **4. Outputting QC reports to summarise pass/fail**

Well done for making it this far! That was quite a lot, and it gets very tedious to check each sample one by one. Let's finally run our nice pipeline, which will do all this hard work for us and give us a nice summary file, with pass/fail info for the steps that matter. 



```bash

```




## **5. What does it look like when things go wrong?**


Remember last week when we actively broke some python scripts? You may remember the error message from the validator script was relatively straightforward. What will happen if we run a dodgy vcf through our mini pipeline script? Let's give it a try. 

```bash
cp ~/data/failing_loudly/patient_1_dodgy.vcfs ./1_vcfs/patient_1_dodgy.vcf
```






### Bonus tasks if time permitting 
1. Take a look at our bash script. "DPYD_mini_pipeline.sh". In that bash scripts, are the paths to the python and awk scripts absolute or relative? Could this cause issues? Can you change the path to be absolute instead? (hint: to get the path of a script or file, you could run):

```bash
 ls -d "$PWD"/{script_or_file_name} )
```

2. In the reports, the following caveat is given: "For the HapB3 genotype “decreased function” is inferred by detecting the exonic tag SNP (c.1236G>A). Recent studies indicate that in rare cases, the causal decreased function variant c.1129-5923C>G may not be present despite having this tag SNP". What mechanism could cause the causal variant not be present in a patient, when the tag SNP itself is?


## Concluding remarks

Hopefully through doing this prac, you will see that even when the stakes for not failing are high, by being loud about everything that could go wrong in as many key places we could, we can assure our clinicians and our patients that everything that the reports have been painstakingly checked and are accurate to the best of our ability. Perhaps the old adage is true, and failing loud really is the key to success! 