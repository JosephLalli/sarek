// Architecture Skeleton for Pangenome Integration
// This is NOT executable code. It illustrates the proposed data flow and structure.

// -----------------------------------------------------------------------------
// 1. PANGENOME ALIGNMENT SUBWORKFLOW
// -----------------------------------------------------------------------------
// Purpose: Align reads using vg giraffe.
// Key Design Principle: "Dumb" alignment. It accepts indices as input, regardless
// of whether they are global (shared) or personalized (unique).

workflow PANGENOME_ALIGN {
    take:
        ch_reads          // [ meta, [reads] ]
        ch_pangenome_idx  // [ meta, gbz, dist, min, xg ] - Can be global or per-sample
        ch_surject_ref    // [ meta, ref_fasta ] - For surjection (CRAM/BAM)

    main:
        // 1. Join reads with indices
        // Since indices are passed with 'meta', we use join on meta.id.
        // For global indices, the upstream logic must have mapped them to the sample meta.
        ch_input = ch_reads.join(ch_pangenome_idx)

        // 2. Align (produce GAM or BAM)
        VG_GIRAFFE(ch_input)
        
        // 3. Surject (if Giraffe produces GAM)
        // If Giraffe produces BAM directly (as seen in some configs), this is optional/skipped.
        // For strictness, we might prefer GAM -> Surject for better control.
        VG_SURJECT(VG_GIRAFFE.out.gam.join(ch_surject_ref))

    emit:
        bam = VG_SURJECT.out.bam // or VG_GIRAFFE.out.bam
}


// -----------------------------------------------------------------------------
// 2. PREPARE PERSONALIZED INDICES (New Subworkflow)
// -----------------------------------------------------------------------------
// Purpose: Encapsulate the complex "Personalized Genome" logic found in the legacy code.
// This extracts it from the alignment logic.

workflow PREPARE_PERSONALIZED_INDICES {
    take:
        ch_reads       // [ meta, [reads] ]
        ch_global_idx  // [ gbz, dist, min, xg ] - The base population graph

    main:
        // 1. Count Kmers
        KMC(ch_reads)

        // 2. Construct Personalized Graph
        // Inputs: Kmers + Global Graph
        MAKE_PERSONALIZED_GBZ(KMC.out.kmers, ch_global_idx)

        // 3. Indexing
        VG_INDEX(MAKE_PERSONALIZED_GBZ.out.gbz)     // -> dist
        VG_MINIMIZER(MAKE_PERSONALIZED_GBZ.out.gbz) // -> min

    emit:
        // Returns the FULL set of indices, but now specific to this sample
        indices = MAKE_PERSONALIZED_GBZ.out.gbz
                    .join(VG_INDEX.out.dist)
                    .join(VG_MINIMIZER.out.min)
                    // .join(xg) // If XG is reused or rebuilt
}


// -----------------------------------------------------------------------------
// 3. MAIN WORKFLOW INTEGRATION (Illustration)
// -----------------------------------------------------------------------------

workflow SAREK_PANGENOME {
    // ... input parsing ...

    // Decision Logic: Personalized or Standard?
    // This could be based on a params flag or meta field.
    
    if (params.use_personalized_genome) {
        
        // A. Run Preparation
        PREPARE_PERSONALIZED_INDICES(ch_reads, ch_global_indices)
        ch_alignment_indices = PREPARE_PERSONALIZED_INDICES.out.indices

    } else {
        
        // B. Use Global Indices
        // We map the global indices to every sample so they can be joined.
        ch_alignment_indices = ch_reads.map { meta, reads -> 
            [ meta, global_gbz, global_dist, global_min, global_xg ] 
        }
    }

    // C. Run Alignment (Identical call for both cases!)
    PANGENOME_ALIGN(ch_reads, ch_alignment_indices, ch_ref_fasta)
    
    // ...
}
