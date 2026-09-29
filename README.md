# BRCA1 Variant Analysis

## Project Overview

This project analyzes BRCA1-related variants from ClinVar using Python.

The main objective is to explore how BRCA1 variants are distributed according to:

- Germline clinical classification
- Variant type
- Molecular consequence
- Pathogenic proportion by variant type

The project also compares examples of Pathogenic, Variant of Uncertain Significance (VUS), and Benign variants in order to better understand the difference between a molecular consequence and a clinical classification.

The analysis was designed as a practical bioinformatics project using real genomic variant data and common Python tools used in biological data analysis.

---

## Biological Background

BRCA1 is a tumor suppressor gene involved in DNA repair and maintenance of genomic stability.

Genetic variants affecting BRCA1 can have different biological and clinical consequences.

### Germline Clinical Classification

ClinVar classifies germline variants into several categories.

- **Pathogenic**: available evidence supports an association with disease.
- **Likely pathogenic**: strong evidence suggests pathogenicity, but the evidence is not considered definitive.
- **Uncertain significance (VUS)**: there is not enough evidence to determine whether the variant is pathogenic or benign.
- **Likely benign**: available evidence suggests that the variant is unlikely to cause disease.
- **Benign**: available evidence supports that the variant is not disease-causing.

A VUS should not automatically be considered dangerous or benign.

### Variant Type

Variant type describes what changed in the DNA sequence.

Examples include:

- Single nucleotide variant
- Deletion
- Insertion
- Duplication
- Indel

### Molecular Consequence

Molecular consequence describes the effect that a variant may have on a transcript or protein.

Examples include:

- **Missense variant**: changes an amino acid.
- **Synonymous variant**: changes the DNA sequence without changing the amino acid.
- **Nonsense variant**: introduces a premature stop codon.
- **Frameshift variant**: changes the reading frame of the coding sequence.
- **Splice donor / splice acceptor variant**: may affect RNA splicing.

It is important to distinguish between:

- Variant type
- Molecular consequence
- Clinical classification

These describe different aspects of the same variant.

A molecular consequence such as a frameshift does not automatically mean that the variant is clinically classified as Pathogenic.

### Multiple Transcripts

A gene can have several transcript isoforms because of alternative splicing.

As a result, the same genomic variant can have different molecular consequences depending on the transcript being considered.

For example:

```text
missense variant|intron variant|non-coding transcript variant