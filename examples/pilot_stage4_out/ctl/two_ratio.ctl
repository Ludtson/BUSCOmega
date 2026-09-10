      seqfile = __SEQFILE__
     treefile = __TREEFILE__
      outfile = __OUTFILE__

        noisy = 3      * screen output verbosity
      verbose = 0      * 0: concise output
      runmode = 0      * 0: user tree (topology fixed, codeml fits branch lengths)

      seqtype = 1      * 1: codons
    CodonFreq = 2      * 0:1/61  1:F1X4  2:F3X4  3:codon table
        clock = 0      * no clock, unrooted tree
        ndata = __NDATA__      * number of gene alignments in this run

        model = 2      * 0: one omega  1: free-ratio  2: two or more ratios (branch)
      NSsites = 0      * 0: one omega among sites
        icode = 0      * 0: universal genetic code

    fix_kappa = 0      * 0: estimate kappa
        kappa = 2.0      * initial kappa
    fix_omega = 0      * 0: estimate omega
        omega = 0.4      * initial omega

    fix_alpha = 1      * no gamma rate variation among sites (M0 / branch models)
        alpha = 0
       Malpha = 0
        ncatG = 3

    fix_blength = 0    * ignore any branch lengths on the input tree; estimate all
       method = 0      * 0: simultaneous optimisation

        getSE = 0
 RateAncestor = 0
    cleandata = 1      * 1: drop any codon column with a gap/ambiguity, once, for every model
