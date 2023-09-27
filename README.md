# WES-analysis 
Project showing the basic steps for Whole Exome Sequence analysis to identify and discuss a single rare Non-synonymous variant
<p>(Please cite this work if you take help from it)<p>

### INTRODUCTION

The objective of this project was to analyze the Whole Exome Sequence(WES) of an Indian Telugu female living in the UK having South Asian ancestry, with an aim to identify and study a rare non-synonymous variant. This was achieved by aligning the exome sequence to the <a rel="noreferrer noopener" href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/" target="_blank">human reference genome</a>(GRCh38) obtained from the 1000 Genomes Project.</p>

#### Individual's sequence ID and information
<p>a) Sequence ID: HG03973</p>

<p>b) Biosample ID: SAME1839728</p>

<p>c) Sample Run ID: ERR250841</p>

<p>d) Population (ethnicity and geographical location): Indian Telugu in the UK, South Asian Ancestry</p>

<p>e) Gender: Female</p>

<p>f) Cell line source: HG03973 at Coriell</p>

### WORKFLOW
<!-- wp:image {"id":11268,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/final-analytical-pipeline.png"><img src="https://gtbinf.files.wordpress.com/2022/05/final-analytical-pipeline.png?w=322" alt="" class="wp-image-11268"/></a><figcaption><em>Figure 1- Analytical pipeline for WES data analysis</em></figcaption></figure>
<!-- /wp:image -->

<p></p>

### METHODOLOGY

<p>- Firstly, the GRCh38 no alt analysis set reference genome sequence file was downloaded from the following website: </p>

<p><code>https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz</code></p>

<p>- The bwa indexed files were downloaded from the same website: </p>

<p><code>https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz</code></p>

<p>- The Whole Exome Sequence of the individual with sequence id HG03973 was downloaded from the 1000Genomes website: </p>

<p><code>https://www.internationalgenome.org/data-portal/sample</code></p>

<p>- The reads were mapped using the bwa mem tool:</p>

<p><code>bwa mem -M -R '@RG\tID:flowcell\tSM: HG03973' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ERR250841_1.fastq ERR250841_2.fastq &gt; HG03973bwamem.sam</code></p>

<p>- Samtools was used to clean up the generated sam file:</p>

<p><code>samtools fixmate -O bam HG03973bwamem.sam HG03973bwamemfixmate.bam</code></p>

<p>- Samtools was then used to sort the bam file:</p>

<p><code>samtools sort -O bam -o HG03973sorted.bam -T /tmp/HG03973temp HG03973bwamemfixmate.bam</code></p>

<p>- In the improvement step, Samtools was used to index sorted bam file:</p>

<p><code>samtools index HG03973sorted.bam</code></p>

<p>- Variant calling was done using the bcftools to generate the vcf file:</p>

<p><code>bcftools mpileup -Ou -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna HG03973sorted.bam | bcftools call -vmO z -o HG03973.vcf.gz</code></p>

<p>- Finally, the vcf file was indexed using tabix:</p>

<p><code>tabix -p vcf HG03973.vcf.gz</code></p>

<p>- The generated unindexed gunzipped vcf file was uploaded in the <a rel="noreferrer noopener" href="https://wannovar.wglab.org/" target="_blank">wANNOVAR</a>[5] online tool for annotation. The following parameter settings were used:</p>

<p>&gt; Reference genome:<em> hg38</em></p>

<p>&gt; Gene Definition: <em>RefSeq Gene</em></p>

<p>&gt; Analysis: <em>Individual</em></p>

<p>- The list of annotated variants obtained through the wANNOVAR tool was first filtered to select variants having exonic function: “<em>nonsynonymous SNV</em>” which left 11,672 variants out of the total 24,608 variants. These leftout variants were filtered for the variants having ClinVar Significance: “<em>Pathogenic/Likely Pathogenic</em>” which left 21 variants. These 21 variants were filtered for the variants having 1000G_ALL and 1000G_SAS frequency &lt; <em>0.01</em>, which left only 2 variants and those 2 variants had 1000Genomes other subpopulation frequencies, ExAC frequency, and ExAC all subpopulation frequencies &lt; <em>0.01</em>. Among these two, the final variant was selected based on the CADD_phred score &gt; <em>30</em>, polyphen-2 score near to <em>1</em> and SIFT score &lt; <em>0.05</em>.</p>

### RESULTS AND DISCUSSION

***Summary table of variants:***

<!-- wp:list {"ordered":true,"type":"a"} -->
<ol type="a"><li>Total number of variants: <em>24,608</em></li></ol>
<!-- /wp:list -->

<!-- wp:image {"id":11275,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image.png?w=975" alt="" class="wp-image-11275"/></a><figcaption><em>
Table 1- Table showing the top 25 variants from all the 24,608 variants ordered with chromosome number </em></figcaption></figure>
<!-- /wp:image -->

<!-- wp:list {"ordered":true,"start":2} -->

<br /><ol start="2"><li>Total number of synonymous variants: <em>12,050</em></li></ol>
<!-- /wp:list -->

<!-- wp:image {"id":11281,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-1.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-1.png?w=975" alt="" class="wp-image-11281"/></a><figcaption><em>
Table 2- Table showing the top 25 variants from all the 12,050 synonymous variants</em></figcaption></figure>
<!-- /wp:image -->

<!-- wp:list {"ordered":true,"start":3} -->
<br /><ol start="3"><li>Total number of non-synonymous variants: <em>11,667</em></li></ol>
<!-- /wp:list -->

<!-- wp:image {"id":11283,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-2.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-2.png?w=975" alt="" class="wp-image-11283"/></a><figcaption><em>
Table 3- Table showing the top 25 variants from all the 11,667 nonsynonymous variants</em></figcaption></figure>
<!-- /wp:image -->

<!-- wp:list {"ordered":true,"start":4} -->
<br /><ol start="4"><li>Number of protein-truncating variants: <em>269</em></li></ol>
<!-- /wp:list -->

<!-- wp:image {"id":11285,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-3.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-3.png?w=975" alt="" class="wp-image-11285"/></a><figcaption><em>
Table 4- Table showing the top 25 variants from all the 269 protein-truncating variants</em></figcaption></figure></p><br />


***Table of potentially damaging or pathogenic nonsynonymous variants:***

<!-- wp:image {"id":11286,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-4.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-4.png?w=975" alt="" class="wp-image-11286"/></a><figcaption><em>
Table 5- List of potentially damaging/pathogenic nonsynonymous variants having ClinVar Significance Pathogenic/Likely Pathogenic along with their SIFT, Polyphen2 HDIV, Polyphen2 HVAR, CADD_raw, CADD_phred scores. Data obtained through the wANNOVAR annotation tool</em></figcaption></figure>
<p></p>

<!-- wp:image {"id":11287,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-5.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-5.png?w=975" alt="" class="wp-image-11287"/></a><figcaption><em>
Table 6- Remaining columns from Table 5. List of potentially damaging/pathogenic nonsynonymous variants having ClinVar Significance Pathogenic/Likely Pathogenic along with their allele frequency in the population of the individual(1000G_SAS) and popmax allele frequency in gnomAD (shown in <strong>Bold</strong>). Data obtained through the wANNOVAR annotation tool</em></figcaption></figure>
<p></p>

<!-- wp:image {"id":11288,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-6.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-6.png?w=975" alt="" class="wp-image-11288"/></a><figcaption><em>Table 7- List of potentially damaging/pathogenic nonsynonymous variants having ClinVar Significance Pathogenic/Likely Pathogenic along with their allele frequency in 1000G_ALL and 1000G different sub-populations</em></figcaption></figure>
<!-- /wp:image -->

<!-- wp:paragraph -->
<p></p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
***Damaging variant information:***
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>Selected single nonsynonymous variant that is predicted to be pathogenic, with low allele frequency (allele frequency less than 0.1%):</p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true,"type":"a"} -->
<ol type="a"><li>Gene: NUP93</li><li>Protein Name: Nuclear pore complex protein Nup93</li><li>Variant ID (rs number): rs145146218</li><li>ClinVar ID: RCV000210641.1</li><li>Exonic function: nonsynonymous SNV<ul><li>Variant genotype (DNA change and amino acid change):</li><li>DNA change: C to T (start and stop position 56831918)</li></ul></li><li>Amino acid change: R(Arginine) to W(Tryptophan) at position 388</li><li>Variant frequency in the overall human population and in different ethnic/regional populations – </li></ol>
<!-- /wp:list -->

<!-- wp:paragraph -->
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; According to 1000Genomes:</p>
<!-- /wp:paragraph -->

<!-- wp:list -->
<ul><li>1000G_All (variant frequency in overall population): 0.001</li><li>1000G_SAS (variant frequency in South Asian population): 0.003</li><li>1000G_EUR (variant frequency in European population): 0.002</li><li>Other populations (1000G_AFRICAN, AMERICAN, EAST ASIAN): 0.0</li><li>The highest frequency is in the South Asian population (T = 0.003)</li></ul>
<!-- /wp:list -->

<!-- wp:image {"id":11291,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-7.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-7.png?w=975" alt="" class="wp-image-11291"/></a><figcaption><em>Table 8– Variant frequency in overall and different human sub-populations according to 1000Genomes</em>[1]</figcaption></figure>
<!-- /wp:image -->

<!-- wp:list {"ordered":true,"start":8} -->
<ol start="8"><li>SIFT score: 0.001</li><li>Polyphen2_HDIV_score: 0.999</li><li>Polyphen2_HVAR_score: 0.94</li><li>CADD_raw: 7.534</li><li>CADD_phred: 34</li></ol>
<!-- /wp:list -->

<!-- wp:paragraph -->
***Damaging variant analysis:***
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>The nuclear pore complex is a massive structure that extends across the nuclear envelope, forming a gateway that regulates the flow of macromolecules between the nucleus and the cytoplasm. Nucleoporins are the main components of the nuclear pore complex in eukaryotic cells. The gene NUP93 encodes a nucleoporin protein that localizes both to the basket of the pore and to the nuclear entry of the central gated channel of the pore. The encoded protein is a target of caspase cysteine proteases that play a central role in programmed cell death by apoptosis. Alternative splicing results in multiple transcript variants encoding different isoforms[2].</p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>The likely effect of the variant is that the mutant protein is unable to constitute nuclear pore complex which was reported through a study on <em>Xenopus Laevis</em> egg extract by Braun et al. (2016)[3] and is responsible for the Steroid-resistant nephrotic syndrome (SRNS) disease of the renal glomerular filter in <em>Homo sapiens </em>(second most frequent cause of end-stage kidney disease (ESKD) in the first 3 decades of life)[3]. This variant gene has been reported in a Serbian girl with nephrotic syndrome type 12 disease who had compound heterozygous mutations in NUP93 gene and each unaffected parent was heterozygous for 1 of the mutations. Different variants of the NUP93 gene have been studied in vitro in the human podocytes by Braun et al. (2016) and some of the mutant proteins were not able to localize properly to the nuclear envelope. R388W mutant protein localized properly to the nuclear envelope in human podocyte cells but was unable to restore nuclear envelope and nuclear pore complex assembly in NUP93-depleted Xenopus egg extracts. The mutation also impaired BMP7 /SMAD4 -dependent gene transcription [3], [4].</p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
***Previous citations on the R388W variant and any other variant on the same gene:***
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>This variant(R388W) has been cited in one paper (PMID: 26878725) which summarizes its predicted role in the Steroid-resistant nephrotic syndrome through the interference with the BMP7 dependent SMAD signaling and the likely inability of Nucleoporin to form nuclear pore complex around the nucleus which is responsible for regulating the flow of macromolecules between the nucleus and the cytoplasm [3].</p>

<p>The same paper (PMID: 26878725) by Braun et al. (2016) reported four more damaging variants (GLY591VAL, TYR629CYS, 1-BP DEL 1326G, IVS13DS G-A +1) on the same gene.</p>

<p>The mutations <span style="text-decoration:underline;">GLY591VAL </span>(reported in 2 siblings, born of consanguineous Turkish parents, with nephrotic syndrome type 12) and <span style="text-decoration:underline;">TYR629CYS </span>(reported in 2 unrelated boys, each born of consanguineous Turkish parents, with nephrotic syndrome type 12) were shown, through in vitro functional studies, to have a role in abrogation of the normal interaction of NUP93 with the phosphorylated, activated forms of SMAD1 and SMAD5 and with the nuclear import receptor IPO7, and impair BMP7 /SMAD4-dependent gene transcription. [3]</p>

<p>The mutation <span style="text-decoration:underline;">1-BP DEL, 1326G</span> (reported in a German girl with nephrotic syndrome type 12), leads to a 1-bp deletion (c.1326delG, NM_014669.4) in exon 12, resulting in a frameshift and premature termination (Lys442AsnfsTer14) forming a truncated protein and the mutation <span style="text-decoration:underline;">IVS13DS, G-A, +1</span> (reported in a German girl with nephrotic syndrome type 12) lead to a G-to-A transition (c.1537+1G-A, NM_014669.4), resulting in a splice site mutation and the in-frame skipping of exon 13. Through in vitro studies, Braun et al. (2016) identified that in both the cases the mutant protein failed to properly localize to the nuclear envelope in human podocytes and was unable to restore nuclear envelope and nuclear pore-complex assembly in NUP93 depleted Xenopus egg extracts. It also abrogated and impaired the normal interaction of NUP93 with the phosphorylated, activated forms of SMAD1 and SMAD5 and with the nuclear import receptor IPO7, and impaired BMP7 /SMAD4-dependent gene transcription. [3]</p>

***Structural model of the protein:***

<!-- wp:image {"id":11294,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-8.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-8.png?w=658" alt="" class="wp-image-11294"/></a><figcaption><em>
</p>Figure 2- Nuclear pore complex protein Nup93. The colors indicate the Model confidence (Navy Blue-Very high (pLDDT* &gt; 90), Sky Blue-Confident (90 &gt; pLDDT &gt; 70), Yellow-Low (70 &gt; pLDDT &gt; 50), Orange-Very low (pLDDT &lt; 50))</em>. *pLDDT- per-residue confidence score</figcaption></figure>
<!-- /wp:image -->

<!-- wp:paragraph -->
<br />***Visualization of altered amino acid on the structural model:***
<!-- /wp:paragraph -->

<!-- wp:image {"id":11296,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-9.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-9.png?w=717" alt="" class="wp-image-11296"/></a><figcaption><em>
</p>Figure 3- Nup93 protein highlighting (red circle) the variant amino acid substitution position</em></figcaption></figure>
<!-- /wp:image -->

<!-- wp:paragraph -->
<p><br />The structures are taken from the AlphaFold Protein Structure Database[6]. AlphaFold produces a per-residue confidence score (pLDDT) between 0 and 100[6]. This position (388) has a score of 89.92 which according to AlphaFold lies in the category of “Confident” score and is considered a good score. The regions surrounding the site of the variant amino acid (position 388) position have a pLDDT of around 90 which is considered a very high confidence score according to the AlphaFold Protein Structure Database[6]. Therefore, this model can overall be considered reliable.</p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>The altered amino acid position is not exactly on the surface but slightly inwards. The variant protein has a single base substitution in the codon at position 388, which is not leading to a protein truncating amino acid, the structure of the protein is less likely to be altered but the substitution of the positively charged amino acid Arginine with the neutral amino acid Tryptophan (in the variant protein) is likely to affect amino acid chain interactions and compromise or alter the secondary and tertiary structures of the protein, which may lead to its abnormal interactions and eventually diseases. It has already been reported[4] that this changed amino acid does change the functioning of the Nup93 protein and is likely involved in the abnormal interaction of NUP93 with the phosphorylated, activated forms of SMAD1 and SMAD5 and with the nuclear import receptor IPO7 in the Steroid-resistant nephrotic syndrome, and impairs BMP7 /SMAD4-dependent gene transcription[3] .</p>
<!-- /wp:paragraph -->

<!-- wp:image {"id":11297,"sizeSlug":"large","linkDestination":"media"} -->
<figure class="wp-block-image size-large"><a href="https://gtbinf.files.wordpress.com/2022/05/image-10.png"><img src="https://gtbinf.files.wordpress.com/2022/05/image-10.png?w=733" alt="" class="wp-image-11297"/></a><figcaption><em>
</p>Figure 4- Slightly rotated view of the protein structure in figure2, for better visualization of the variant amino acid position</em></figcaption></figure>
<!-- /wp:image -->

<br />

### REFERENCES
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>[1]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; “rs145146218 RefSNP Report - dbSNP - NCBI.” https://www.ncbi.nlm.nih.gov/snp/rs145146218/#frequency_tab (accessed Nov. 29, 2021).</p>

<p>[2]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; “NUP93 nucleoporin 93 [Homo sapiens (human)] - Gene - NCBI.” https://www.ncbi.nlm.nih.gov/gene/9688 (accessed Nov. 23, 2021).</p>

<p>[3]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; D. A. Braun <em>et al.</em>, “Mutations in nuclear pore genes NUP93, NUP205 and XPO5 cause steroid-resistant nephrotic syndrome,” <em>Nature Genetics</em>, vol. 48, no. 4, pp. 457–465, Mar. 2016, doi: 10.1038/NG.3512.</p>

<p>[4]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; “OMIM Entry - * 614351 - NUCLEOPORIN, 93-KD; NUP93.” https://www.omim.org/entry/614351?search=nup93&amp;highlight=nup93 (accessed Nov. 23, 2021).</p>

<p>[5]&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; “wANNOVAR.” https://wannovar.wglab.org/ (accessed Nov. 29, 2021).</p>

<p>[6]          “AlphaFold Protein Structure Database.” https://alphafold.ebi.ac.uk/entry/Q8N1F7 (accessed Nov. 27, 2021).</p>
