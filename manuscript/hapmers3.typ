#import "@preview/arkheion:0.1.0": arkheion, arkheion-appendices

#show: arkheion.with(
  title: "Using readmers and hapmers in assessing phase switching after read error correction of Oxford Nanopore Sequences",
  authors: (
    (name: "Jean P. Elbers", email: "jean.elbers@gmail.com", affiliation: "Institute of Medical Genetics, Center for Pathobiochemistry and Genetics, Medical University of Vienna"),
    (name: "David Horner", email: "david.horner@meduniwien.ac.at", affiliation: "Institute of Medical Genetics, Center for Pathobiochemistry and Genetics, Medical University of Vienna"),
    (name: "Tamara Löwenstern", email: "tamara.loewenstern\n@meduniwien.ac.at", affiliation: "Institute of Medical Genetics, Center for Pathobiochemistry and Genetics, Medical University of Vienna"),
    (name: "Franco Laccone", email: "franco.laccone@meduniwien.ac.at", affiliation: "Institute of Medical Genetics, Center for Pathobiochemistry and Genetics, Medical University of Vienna")
  ),
  abstract: [Methods for sequence error correction can improve sequence accuracy; however, there can be unintended errors added during error correction. One such example is phase switching, whereby sequences derived from genomes containing more than one parental copy have contributions from more than one parental haplotype. Such switches are mistakes that can confound downstream analyses especially _de novo_ genome assembly. While DNA sequences do not possess words such as linguistic languages, one can partition a DNA sequence into word-like objects called k-mers. K-mers include pieces of DNA sequence of k length that can describe various properties of DNA sequences. With regard to phase switching, there are k-mers present in one parental haplotype not found in other(s). These so-called hapmers can represent inaccuracies in that parental haplotype's assembly but also correct, unique DNA variation. Here we investigated the effect of DNA sequence error correction on phase switching at the sequence/read level. Using several error-correction methods, we find all methods tested are similar to raw, presumably, phase-switch-free Oxford Nanopore Technologies (ONT) sequences in the percentage of readmers (k-mers from the ONT sequences) matching one parental haplotype's hapmers. This work demonstrates an efficient method to assess if an error-correction method has introduced phase switching implemented in the Julia programming language.],
  keywords: ("k-mers", "hapmers", "Oxford Nanopore Technologies", "ONT", "phase switches"),
  date: "October 18, 2024",
)

#set cite(style: "chicago-author-date")
#show link: underline

= Introduction
The advancement of sequencing technologies, such as Oxford Nanopore Technologies (ONT) sequencing, has revolutionized genomic research by enabling rapid and cost-effective long-read DNA sequencing. However, ONT sequencing reads have a tendency to have insertion and deletion (indels) errors, especially when the template sequence has homopolymers @delahaye2021sequencing@huang2021homopolish. Indel sequencing errors can be corrected by a variety of computational methods, but it is not always clear whether non-intentional errors are introduced during the indel error-correction process as no method is perfect @huang2021homopolish@wang2020nanoreviser. One such type of unintentionally-introduced error is phase switching, a phenomenon where reads derived from multiple parental genome copies contain a mixture of more than one parental copy's DNA nucleotides @merqury. Phase switching may not be detrimental to all analyses, but it can confound downstream applications such as _de novo_ genome assembly, whereby there could be loss of variation such as SNPs for phasing a genomic region resulting in that region being unphaseable (Konstantinos Kyriakidis, _pers. comm._).
#parbreak()
While DNA sequences do not possess words such as linguistic languages, one can partition a DNA sequence into word-like objects called k-mers @waterman198852. K-mers include pieces of DNA sequence of k length that can describe various properties of DNA sequences. With regard to phase switching, there are k-mers present in one parental haplotype not found in other(s). These so-called hapmers can represent inaccuracies in that parental haplotype's assembly but also correct, unique DNA variation (Heng Li, _pers. comm._). Further, during phase switching, there is a mixture of hapmers from heterozygous genomics regions that do not coexist on the same homologous chromosome.
#parbreak()
While there are computational methods for estimating phase switching at the assembled contig or chromosome level @yak@merqury, here we investigate the effect of DNA sequence error correction on phase switching at the sequence/read level by comparing the k-mers in ONT reads (readmers) and percentage of hapmers from one or the other parental haplotype with human cell-line DNA. Using a computational approach implemented in the Julia programming language, we present a systematic evaluation of several error-correction methods.

= Materials Methods
== Cell Culture
We cultured HG002/RM8391 cells (NIST, Coriell), following the manufacturer's culturing and expanding guidelines.

== DNA Extraction, Library Prep, Sequencing
We extracted DNA from the HG002-cultured cells using automated extraction with EZ1 Advanced XL Robot (Qiagen) and EZ1 DNA Tissue Kit (Qiagen) and then performed left-size selection using Short Read Eliminator (SRE) Kit (Circulomics) to enrich reads greater than 25,000 bases. We incubated the samples with the SRE-reagent in a ratio 1:1 for 1 hour at 50°C to increase DNA yields. Next, we used the ONT SQK-LSK114 kit for ONT library preparation and the ONT-protocol for “Ligation sequencing 30kb human gDNA” (Version: GDH_9174_v114, revD, 10.11.2022) starting from the step “DNA repair and end-prep”. Finally, we sequenced the libraries on two R10.4.1 ONT PromethION flowcells using an ONT P2 Solo device connected to an ONT GridION.

== Read Processing
We basecalled the resulting pod5 squiggles (i.e., the current raw output files from ONT sequencing) with Dorado 0.6.0 @ONT2024 using the dna\_r10.4.1\_e8.2\_400bps\_sup\@v4.3.0 basecalling model to initially generate SUP, raw ONT reads.

== Error correction
While we are well aware there are numerous raw ONT sequence error-correction tools, we focused on the following as these were used in a prior analysis and were assessed to see if they might negatively impact phase switching at the sequence/read level. During error correction, we made a Snakemake workflow @snakemake for processing and also used BBMap/BBTools's scripts @Bushnell such as partition.sh version 39.01, shred.sh version 39.01, and filterbyname.sh version 39.06 for various operations such as subsetting the reads, ensuring reads > 1,000,000 bases were split for effective downstream processing by one of the error-correction tools, Peregrine 2021 @Peregrine_2021, and downstream comparison of the same reads, respectively.

=== Herro
#cite(<Stanojevic>, form: "prose") developed the deep-learning, haplotype-aware ONT error-correction tool Herro. We used Herro commit \# 5e7b70f and the version 0.1 R10.4.1 read-correction model. We first used seqkit 2.5.1 @Shen to obtain read identifiers and minimap2 version 2.26-r1175 @Li1 to perform read overlapping, running the Herro developers' create_batched_alignments.sh script. We then ran Herro inference (i.e., error correction) using a batch size of 64 and a single Nvidia A100 40GB graphical processing unit.

=== Brutal Rewrite
Brutal Rewrite is an error-correction method that uses k-mers for error correction @Br. We used Brutal Rewrite commit \# ad87f92, the graph algorithm for de Bruijn graph @pevzner1989tuple based error correction, and a k-mer length of 19 on Herro-corrected reads.

=== Peregrine 2021
Peregrine 2021 @Peregrine_2021 is genome assembler similar to the Peregrine genome assembler @Peregrine with both assemblers using sparse-hierarchical minimizers. We used Peregrine 2021 version 0.4.13 performing one round of overlapped-based error correction using the Brutal Rewrite-corrected reads as input and 6, 56, and 80 as values for the reduction factor, k-mer size, and window size for Peregrine 2021 settings, respectively.

=== DeChat
DeChat is an error-correction method that we used on corrected reads but is designed for raw ONT reads @Li2024.05.09.593079. We used default settings for version 1.0.0, using the Peregrine 2021-corrected reads as input and retaining only the _recorrected.fa_ read DeChat output.

== Phase switching analyses
=== Read alignment
We used minimap2 version 2.28-r1209 with the following options to map raw, Herro, Brutal Rewrite, Peregrine 2021, and DeChat reads against either the maternal or paternal haplotypes separately from _hg002v1.0.1.fasta_ @HG002: output "X" and "=" in cigar strings, no secondary reads, soft-clipping of supplementary alignments, base-level alignment, output SAM format, and the high-quality, long-read preset. We then used SAMtools 1.9 @samtools to retain only primary alignments and those with mapping qualities equal to 60 and converted those SAM alignments into pairwise-mapping format (PAF) alignments with minimap2's paftools.js sam2paf subprogram. Next, we retained the longest, primary alignment from each read and excluded those reads that Herro split into parts for Herro, Brutal Rewrite, Peregrine 2021, and DeChat reads.

=== K-mers and hapmers
We developed a Julia programming language script @julia that takes read alignment information from minimap2 PAF files, FASTA files, and FASTA-index files of the reads and extracts ONT read k-mers (hereafter readmers) and associated k-mers unique to each haplotype (hereafter hapmers). We required that for a read alignment to be used for readmer and hapmer analysis, that query start and query stop had to be identical for maternal and paternal alignments (i.e., the same region of the read aligned to both haplotypes). Next, we iterated over each alignment meeting our criteria to extract the readmers and noted the total hapmers and matching hapmers for both maternal and paternal haplotypes. We processed the resulting tables with a Julia programming language script, converting the number of matches to percentages. For each read alignment, we aggregated the difference in the percentages of matching hapmers by taking the absolute value of the difference between the maternal and paternal matching percentages. For example, if we had a read alignment with 75% maternal hapmer matches and 25% paternal hapmer matches, then we would have 50% as the resulting value for that read alignment and likewise 50% if percentages were flipped between haplotypes (i.e., we would not have -50%). Finally, we generated histograms with a bin width of 1 from 1% to 100% with a cut-off of 6,000 read alignments per bin on y-axis or no cut-off. For all analyses we used a k-mer length of 15.

=== Alignment statistics
Using the same reads for hapmer analysis, we realigned these sequences to both _hg002v1.0.1.fasta_ assembly haplotypes simultaneously using minimap2 version 2.28-r1209 with the following options: output "X" and "=" in cigar strings, no secondary reads, soft-clipping of supplementary alignments, base-level alignment, output SAM format, and the high-quality, long-read preset back to an assembly preset (i.e., lr:hqae and not lr:hq). We then used the program Bam Error Stats Tools ("best", commit \#c1c69bb, @best) to analyze errors in the aligned reads.

= Results
We generated 15,048,314 raw, super (i.e., SUP) accuracy reads that we error corrected with Herro, resulting in 4,578,144 reads of which 582,208 were split parts of reads and not analyzed for phase-switching analysis. Brutal Rewrite does not filter out reads, but there were 4,490,689 reads remaining after one round of overlapping with Peregrine 2021 of which 532,757 were split parts and not analyzed for phase-switching analysis. DeChat did not filter out any reads from Peregrine 2021. @table1 shows results for the number of read alignments meeting criteria for each correction method, prior to requiring the same read identifier for all 5 datasets for downstream hapmer analysis.

In general, the percentage of read alignments passing our filtering criteria of mapping quality 60 and only the longest, primary alignment increased with each successive correction method, except there were lower percentages for Brutal Rewrite than Herro (@table1).

#show figure.caption: strong

#show figure.where(
  kind: table
): set figure.caption(position: top)

#figure(
  caption: [Number of read alignments and percent out of total reads after filtering for alignments with mapping quality 60 and the longest, primary alignments. Herro, Brutal Rewrite, Peregrine 2021, and DeChat read alignments were further filtered for any reads split by the Herro algorithm. Coverage is estimated based on a haploid genome size of 3.1 billion bases],
  table(
    align: center,
    columns: (auto, auto, auto, auto, auto, auto),
    row-gutter: (2pt, auto, auto, auto, auto, auto),
    stroke: 0.5pt,
    inset: 5pt,
    [], [Total Reads], [Maternal Alignments], [Paternal Alignments], [N50], [Coverage],
    [Raw], [15,048,314], [4,029,514 (26.7772%)], [3,962,042 (26.3288%)], [21,268], [40.3x],
    [Herro], [4,578,144], [1,733,156 (37.8572%)], [1,705,413 (37.2512%)], [24,437], [26.3x],
    [Brutal Rewrite], [4,578,144], [1,733,152 (37.8571%)], [1,705,385 (37.2506%)], [24,437], [26.3x],
    [Peregrine 2021], [4,490,689], [1,724,488 (38.4014%)], [1,697,130 (37.7922%)], [24,152], [25.7x],
    [DeChat], [4,490,689], [1,724,727 (38.4067%)], [1,697,341 (37.7969%)], [24,153], [25.7x],
  )
) <table1>

There were about 1,145,743 reads alignments that had the same query start and query stop for maternal and paternal haplotypes and were present in the raw, Herro, Brutal Rewrite, Peregrine 2021, and DeChat datasets. This 5-way intersection allowed us to keep reads identical but change merely what correction if any was applied to examine readmers matching hapmers.

Looking at @figure1, one can see that the majority of the read alignments had close to 100% hapmers from a single parental haplotype suggesting little to no phase switching. The Peregrine 2021 method introduced a few taller, lower-percentage peaks, which could indicate a slight increase in phase switching at the read level. That being said, this was only a tiny increase relative to the almost 1 million read alignments in the rest of the histogram (@figure2).

#show figure.caption: strong

#figure(
  image("hapmers-histogram-Newsame_reads2.svg", width: 110%),
  caption: [Histograms of absolute value of difference in matching hapmer percentage between maternal and paternal HG002 haplotypes. The y-axis was cut-off at 6,000 read alignments for visualization of low frequency data. The vast majority of occurrences were close to 100% suggesting that most of the read alignments possessed hapmers from a single haplotype and little to no phase switching.]
) <figure1>


#show figure.caption: strong

#figure(
  image("hapmers-histogram-Newsame_reads-full2.svg", width: 110%),
  caption: [Histograms of absolute value of difference in matching hapmer percentage between maternal and paternal HG002 haplotypes. The y-axis is not cut-off and shows the staggering number of read alignments with almost 100% hapmers from a single haplotype.]
) <figure2>

@figure2 shows in greater detail such a high proportion of the total read alignments coming from a single haplotype (i.e., how many read alignments had approximately 100% matching hapmers to a single haplotype).

There were some interesting scenarios whereby almost 12% analyzed, read alignments aligned in regions where 0 hapmers (i.e., homozygous genomic regions) were present (@table2). Keeping in mind the distinction between a matching hapmer (i.e., a hapmer in a heterozygous genomic region matching a readmer in a read alignment) and a hapmer (i.e., a hapmer in a heterozygous region), we also observed read alignments with 0 matching hapmers but >0 total hapmers in the genomic region (@table2).

#show figure.caption: strong

#show figure.where(
  kind: table
): set figure.caption(position: top)

#figure(
  caption: [Number of read alignments (out of 1,145,743 total) based on whether the reads were aligned in heterozygous or homozygous regions.],
  table(
    align: center,
    columns: (auto, auto, auto, auto),
    row-gutter: (2pt, 2pt, auto, auto, auto, auto, auto),
    stroke: 0.5pt,
    inset: 5pt,
    [], table.cell(colspan: 2)[Read Alignments in Heterozygous Regions #parbreak() (>0 Total Hapmers)], [Read Alignments in Homozygous Regions #parbreak() (0 Total Hapmers)],
    [], [Readmers Matching Hapmers], table.cell(colspan: 2)[Readmers Matching 0 Hapmers],
    [Raw], [1,009,418], [7,483], [128,842], 
    [Herro], [1,009,411], [2,232], [134,100],
    [Brutal Rewrite], [1,009,411], [2,232], [134,100], 
    [Peregrine 2021], [1,007,110], [2,084], [136,549],
    [DeChat], [1,007,150], [2,044], [136,549],
  )
) <table2>

#show figure.caption: strong

#show figure.where(
  kind: table
): set figure.caption(position: top)

#figure(
  caption: [Alignment statistics from the program Best (Bam Error Stats Tools) of reads used in hapmer analyses. Abbreviations: Compres for compressed, Kbp for Kilobase pair, Inser for Insertion, Del for Deletion, Hp for homopolymer.],
  table(
    align: center,
    columns: (auto, auto, auto, auto),
    row-gutter: (2pt, auto, auto, auto, auto, 2pt, 2pt, auto, auto, auto, auto, 2pt, 2pt, auto, auto, auto, auto, auto),
    stroke: 0.5pt,
    inset: 5pt,
    [], [Identity], [Identity Quality Value], [Gap Compres Identity],
    [Raw], [0.980335], [17.062956], [0.984901],
    [Herro], [0.999637], [34.396470], [0.999770],
    [Brutal Rewrite], [0.999637], [34.398948], [0.999772],
    [Peregrine 2021], [0.999726], [35.619356], [0.999830],
    [DeChat], [0.999741], [35.862487], [0.999842],
    [], [Matches per Kbp], [Mismatches per Kbp], [Non-hp Inser per Kbp],
    [Raw], [985.021819], [7.253553], [2.826602],
    [Herro], [999.725438], [0.052275], [0.033748],
    [Brutal Rewrite], [999.725298], [0.051433], [0.033581],
    [Peregrine 2021], [999.790755], [0.031596], [0.024639],
    [DeChat], [999.803664], [0.027125], [0.023584],
    [], [Non-hp Del per Kbp], [Hp Inser per Kbp], [Hp Del per Kbp],
    [Raw], [3.413819], [1.954716], [4.310809],
    [Herro], [0.044347], [0.055095], [0.177940],
    [Brutal Rewrite], [0.044946], [0.054916], [0.178323],
    [Peregrine 2021], [0.026350], [0.040332], [0.151300],
    [DeChat], [0.022643], [0.039365], [0.146568],
  )
) <table3>

Using the same reads used in the hapmer analysis but looking at the alignment statistics such as alignment identity, matches/mismatches per kilobase, and non-homopolymer/homopolymer insertion/deletions per kilobase: we saw that the raw reads had much lower identities, matches per kilobase, and higher values for other categories compared to all other corrected reads. Each successive correction method seemed to improve the alignment statistics, but one notable thing was that Brutal Rewrite had a lower non-homopolymer insertion rate per kilobase than Herro, but a higher non-homopolymer deletion rate per kilobase. The same difference between insertion and deletion rates occurred for homopolymers between Brutal Rewrite and Herro.

= Discussion
#parbreak()
== General
The findings presented here provide significant insights into the effects of various error-correction methods on phase switching at the read level in ONT sequencing. One of the key observations is that all tested error-correction methods—Herro, Brutal Rewrite, Peregrine 2021, and DeChat—exhibited similar percentages of readmers matching one parental haplotype's hapmers compared to raw ONT sequences. This suggests that while these methods improve overall sequence accuracy, they may not significantly exacerbate phase switching issues.
== Phase Switching Insights
The results indicate that the majority of read alignments primarily matched hapmers from a single parental haplotype, which suggests minimal phase switching occurred. However, the presence of taller, lower-percentage peaks in the Peregrine 2021 reads (@figure1) raises questions about Peregrine 2021's haplotype awareness. Although this increase was minor, it highlights the need for further investigation into how certain error-correction tools can inadvertently introduce errors during the correction process. Peregrine 2021 is very fast in read overlapping as it avoids all-versus-all read overlapping (i.e., the computationally expensive process of overlapping each and every single read) through its unique, sparse-hierarchical minimizer index that helps find overlap candidates instead all-versus-all overlapping. Perhaps not doing all-versus-all overlapping increases increases phase-switching ever so slightly. We experimented with alternative values for the reduction factor and window size when running Peregrine 2021, we used 4 and 64, respectively, but this did not change the overall pattern observed in the initial Peregrine 2021 hapmer analysis. As Peregrine 2021 uses hand-coded heuristics in error correction, it may not have performed as well as Herro, which is a deep learning based method, especially given that deep learning methods outperform hand-coded heuristic methods on low coverage data used here (Jason Chin, _pers. comm._).
#parbreak()
Interestingly, the analysis revealed a subset of read alignments that had no matching hapmers to readmers in regions where hapmers existed. This could be indicative of sequencing errors from those genomic regions, misalignments, or perhaps even low-levels of contamination. While we did not further investigate these situations, understanding these discrepancies could help refine error-correction strategies.
== Error Rates in Different Correction Methods
An intriguing finding is the observed difference in non-homopolymer insertion and deletion rates between the Brutal Rewrite and Herro methods. Specifically, Brutal Rewrite demonstrated a lower non-homopolymer insertion rate per kilobase but a higher non-homopolymer deletion rate compared to Herro. This raises questions about the underlying mechanics of these error-correction algorithms. It may suggest that Brutal Rewrite is more conservative in adding new bases but potentially more aggressive in removing bases that it deems erroneous. As our estimated coverage was only 26.3x for the reads used in Brutal Rewrite correction, perhaps higher coverage would have been beneficial. Further in regions of low-complexity repeats, Brutal Rewrite tends to select the shortest repetitive sequence (Antoine Limasset and Pierre Marjion, _pers. comm._).
#parbreak()
Similar patterns were observed in terms of homopolymer regions, where the insertion and deletion rates differed between the two methods. These variations highlight the importance of understanding the specific behaviors of each correction tool, as they may influence downstream analyses, such as variant calling and _de novo_ genome assembly.
== Future Directions
Future studies should focus on further characterizing the mechanisms behind phase switching and error introduction across different correction methods. A more comprehensive understanding of how each tool interacts with specific genomic features—such as homozygosity and heterozygosity—could lead to improved error-correction algorithms that maintain haplotype integrity while enhancing sequence accuracy.
== Conclusion
In conclusion, while current methods show promise in maintaining haplotype fidelity and reducing sequencing errors, ongoing research is essential to refine these tools and their applications in genomics. Addressing the discrepancies noted in this study will be crucial for advancing our understanding of error-correction for ONT sequences.

= Acknowledgements
We would like to acknowledge the Medical University of Vienna's High Performance Computing cluster for computing resources, Konstantinos Kyriakidis for thoughts regarding the effect of phase switching on _de novo_ genome assembly, Heng Li for thoughts regarding hapmer and readmer error rates, Jason Chin for thoughts regarding Peregrine-2021, and Antoine Limasset and Pierre Marjion for thoughts regarding Brutal Rewrite.

= Data Availabilty
Code documenting how analyses were conducted is available at https://github.com/jelber2/hapmers . Raw, basecalled SUP accuracy reads are available in unaligned BAM format at the following Zenodo repository https://doi.org/10.5281/zenodo.13841954 .

= Author contributions
JPE executed all analyses and wrote the manuscript. DH cultured HG002 cells. TL performed DNA extraction, library preparation, and sequencing. FL reviewed the manuscript.

// Add bibliography and create Bibiliography section
#bibliography("bibliography.bib")