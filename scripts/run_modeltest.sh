#!/bin/bash

#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --export=ALL

printf ">Ced4\nREYHVDRVIKKLDEMCDLDSFFLFLHGRAGSGKSVIASQALSKSDQLIGINYDSIVWLKDSGTAPKSTFDLFTDILLMLARVVSDTDDSHSITDFINRVLSRSEDDLLNFPSVEHVTSVVLKRMICNALIDRPNTLFVFDDVVQEETIRWAQELRLRCLVTTRDVEISNAASQTCEFIEVTSLEIDECYDFLEAYGMPMPVGEKEEDVLNKTIELSSGNPATLMMFFKSCEPKTFEKMAQLNNKLESRGLVGVECITPYSYKSLAMALQRCVEVLSDEDRSALAFAVVMPPGVDIPVKLWSCVIPVD\n" > $TMPDIR/all_nbarc.fa

cat results/resistify/final/nbarc.fasta >> $TMPDIR/all_nbarc.fa

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/mafft:7.525--h031d066_1 mafft \
    --thread 8 \
    --auto \
    $TMPDIR/all_nbarc.fa \
    > $TMPDIR/all_nbarc.msa

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/modeltest-ng:0.1.7--h103dbdd_2 modeltest-ng \
  -i $TMPDIR/all_nbarc.msa \
  -d aa \
  -p 8 \
  -o $TMPDIR/modeltest-ng

singularity exec -B /mnt/:/mnt/ docker://quay.io/biocontainers/raxml-ng:1.2.2--h6d1f11b_0 raxml-ng \
  --force perf_threads \
  --threads 16 \
  --msa all_nbarc.msa \
  --model JTT+G4