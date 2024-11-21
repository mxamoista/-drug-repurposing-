#  Drug repurposing

## #1 Introduction

 Drug repurposing offers a shorter approval process than new drug development. We therefore searched large public datasets of drug-induced gene expression signatures to identify agents that might be effective against specific subtype of osteoarthritis(OA). Within this framework for OA drug  repurposing, we established  C1- and C4-specific transcriptomic signatures of OA concomitant glycosaminoglycan loss and inflammation respectively and then applied the computational drug repurposing algorithm to query the LINCS Consortium,  which includes more than a 1,000-fold scale-up of the Connectivity Map(CMap). 

## #2 Computational drug repurposing screen

Step one,we prepare signatures for drug predictions against "C1" and "C4" in the LINCS Consortium in the code of "build lincs signature.R".

Step two ,in "C1_LINCS_DESEQ2.R" and "C4_LINCS_DESEQ2.R",we first get differential genes towards C1 and C4 using DESeq2.Previous analyses in studies of  cancer drug repurposing indicate that, within the CMap or LINCS database,  compounds found to more dramatically ‘flip’ the transcriptomic sig nature of the cancer back toward a normal state were more likely to  be effective in clinical studies.  Therefore, we calculated the CMap score of the "flip" for all compounds combining the differential genes with the signatures generated in step one in the LINCS consortium against C1- and C4-specific transcriptomic signatures of OA. 

Step three,visualization of results.

![](E:\Submit\JCI\Next Generation Sequencing Workflow.png)

