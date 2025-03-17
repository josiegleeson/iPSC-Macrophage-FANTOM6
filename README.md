# iPSC-Macrophage-FANTOM6 paper

This repository contains two R scripts that were used in the iPSC-Macrophage FANTOM6 paper: (link here).

Bambu was used without a reference GTF to discover unannotated transcripts. To re-name known transcripts that were assembled using Bambu, we first ran GFFCompare and then re-named transcripts matching a known reference transcript back to their original transcript ID. Novel transcripts were also re-named so that their IDs did not contain 'Bambu', allowing us to re-run Bambu a second time alongside the reference GTF. To re-name these we used the script 'match_denovo_ids_to_reference_ids.R'.

After running Bambu the second time, the transcripts in the output GTF and counts files were re-named for consistency so that all novel transcripts began with 'BambuTx'. To re-name these transcripts we used the script 'rename_bambu_transcript_ids.R'.
