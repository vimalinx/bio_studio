# Biomni Tool Index

This document lists the available tools within the Biomni package, organized by module.


## Biochemistry (`biochemistry.py`)

| Function | Description |
| --- | --- |
| `analyze_circular_dichroism_spectra` | Analyzes circular dichroism (CD) spectroscopy data to determine secondary structure and thermal s... |
| `analyze_rna_secondary_structure_features` | Calculate numeric values for various structural features of an RNA secondary structure. |
| `analyze_protease_kinetics` | Analyze protease kinetics data from fluorogenic peptide cleavage assays. |
| `analyze_enzyme_kinetics_assay` | Performs in vitro enzyme kinetics assay and analyzes the dose-dependent effects of modulators. |
| `analyze_itc_binding_thermodynamics` | Analyzes isothermal titration calorimetry (ITC) data to determine binding affinity and thermodyna... |
| `analyze_protein_conservation` | Perform multiple sequence alignment and phylogenetic analysis to identify conserved protein regions. |

## Bioengineering (`bioengineering.py`)

| Function | Description |
| --- | --- |
| `analyze_cell_migration_metrics` | Analyze cell migration metrics from time-lapse microscopy images. |
| `perform_crispr_cas9_genome_editing` | Simulates CRISPR-Cas9 genome editing process including guide RNA design, delivery, and analysis. |
| `analyze_calcium_imaging_data` | Analyze calcium imaging data to quantify neuronal activity metrics. |
| `analyze_in_vitro_drug_release_kinetics` | Analyzes in vitro drug release kinetics from biomaterial formulations. |
| `analyze_myofiber_morphology` | Quantifies morphological properties of myofibers in microscopy images of tissue sections. |
| `decode_behavior_from_neural_trajectories` | Model neural activity trajectories and decode behavioral variables. |
| `simulate_whole_cell_ode_model` | Simulate a whole-cell model represented as a system of ordinary differential equations (ODEs). |

## Bioimaging (`bioimaging.py`)

| Function | Description |
| --- | --- |
| `split_modalities` | Convenience function for splitting modalities |
| `prepare_input_for_nnunet` | Convenience function for preparing nnUNet input |
| `segment_with_nn_unet` | Convenience function for nnUNet segmentation |
| `create_segmentation_visualization` | Convenience function for segmentation visualization |
| `preprocess_image` | Standalone image preprocessing function for Biomni integration. |
| `quick_rigid_registration` | Quick rigid registration function for Biomni integration. |
| `quick_affine_registration` | Quick affine registration function for Biomni integration. |
| `quick_deformable_registration` | Quick deformable registration function for Biomni integration. |
| `batch_register_images` | Batch registration of multiple images to a single reference. |
| `calculate_similarity_metrics` | Calculate similarity metrics between two images. |

## Biophysics (`biophysics.py`)

| Function | Description |
| --- | --- |
| `predict_protein_disorder_regions` | Predicts intrinsically disordered regions (IDRs) in a protein sequence using IUPred2A. |
| `analyze_cell_morphology_and_cytoskeleton` | Quantifies cell morphology and cytoskeletal organization from fluorescence microscopy images. |
| `analyze_tissue_deformation_flow` | Quantify tissue deformation and flow dynamics from microscopy image sequence. |

## Cancer Biology (`cancer_biology.py`)

| Function | Description |
| --- | --- |
| `analyze_ddr_network_in_cancer` | Analyze DNA Damage Response (DDR) network alterations and dependencies in cancer samples. |
| `analyze_cell_senescence_and_apoptosis` | Analyze flow cytometry data to quantify senescent and apoptotic cell populations. |
| `detect_and_annotate_somatic_mutations` | Detects and annotates somatic mutations in tumor samples compared to matched normal samples. |
| `detect_and_characterize_structural_variations` | Detects and characterizes structural variations (SVs) in genomic sequencing data. |
| `perform_gene_expression_nmf_analysis` | Performs Non-negative Matrix Factorization (NMF) on gene expression data to extract |
| `analyze_copy_number_purity_ploidy_and_focal_events` | CNVkit-based copy number workflow: CNV segmentation, purity/ploidy & HRD approximation, focal eve... |

## Cell Biology (`cell_biology.py`)

| Function | Description |
| --- | --- |
| `quantify_cell_cycle_phases_from_microscopy` | Quantify the percentage of cells in each cell cycle phase using Calcofluor white stained microsco... |
| `quantify_and_cluster_cell_motility` | Quantify cell motility features from time-lapse microscopy images and cluster cells based on moti... |
| `perform_facs_cell_sorting` | Performs Fluorescence-Activated Cell Sorting (FACS) to enrich cell populations based on fluoresce... |
| `analyze_flow_cytometry_immunophenotyping` | Analyze flow cytometry data to identify and quantify specific cell populations based on surface m... |
| `analyze_mitochondrial_morphology_and_potential` | Quantifies metrics of mitochondrial morphology and membrane potential from fluorescence microscop... |

## Database (`database.py`)

| Function | Description |
| --- | --- |
| `get_hpo_names` | Retrieve the names of given HPO terms. |
| `_query_llm_for_api` | Helper function to query LLMs for generating API calls based on natural language prompts. |
| `_query_rest_api` | General helper function to query REST APIs with consistent error handling. |
| `_query_ncbi_database` | Core function to query NCBI databases using Claude for query interpretation and NCBI eutils. |
| `_format_query_results` | A general-purpose formatter for query function results to reduce output size. |
| `query_uniprot` | Query the UniProt REST API using either natural language or a direct endpoint. |
| `query_alphafold` | Query the AlphaFold Database API for protein structure predictions. |
| `query_interpro` | Query the InterPro REST API using natural language or a direct endpoint. |
| `query_pdb` | Query the RCSB PDB database using natural language or a direct structured query. |
| `query_pdb_identifiers` | Retrieve detailed data and/or download files for PDB identifiers. |
| `query_kegg` | Take a natural language prompt and convert it to a structured KEGG API query. |
| `query_stringdb` | Query the STRING protein interaction database using natural language or direct endpoint. |
| `query_iucn` | Query the IUCN Red List API using natural language or a direct endpoint. |
| `query_paleobiology` | Query the Paleobiology Database (PBDB) API using natural language or a direct endpoint. |
| `query_jaspar` | Query the JASPAR REST API using natural language or a direct endpoint. |
| `query_worms` | Query the World Register of Marine Species (WoRMS) REST API using natural language or a direct en... |
| `query_cbioportal` | Query the cBioPortal REST API using natural language or a direct endpoint. |
| `query_clinvar` | Take a natural language prompt and convert it to a structured ClinVar query. |
| `query_geo` | Query the NCBI Gene Expression Omnibus (GEO) using natural language or a direct search term. |
| `query_dbsnp` | Query the NCBI dbSNP database using natural language or a direct search term. |
| `query_ucsc` | Query the UCSC Genome Browser API using natural language or a direct endpoint. |
| `query_ensembl` | Query the Ensembl REST API using natural language or a direct endpoint. |
| `query_opentarget` | Query the OpenTargets Platform API using natural language or a direct GraphQL query. |
| `query_monarch` | Query the Monarch Initiative API using natural language or a direct endpoint. |
| `query_openfda` | Query the OpenFDA API using natural language or direct parameters. |
| `query_gwas_catalog` | Query the GWAS Catalog API using natural language or a direct endpoint. |
| `query_gnomad` | Query gnomAD for variants in a gene using natural language or direct gene symbol. |
| `blast_sequence` | Identifies a DNA sequence using NCBI BLAST with improved error handling, timeout management, and ... |
| `query_reactome` | Query the Reactome database using natural language or a direct endpoint. |
| `query_regulomedb` | Query the RegulomeDB database using natural language or direct variant/coordinate specification. |
| `query_pride` | Query the PRIDE (PRoteomics IDEntifications) database using natural language or a direct endpoint. |
| `query_gtopdb` | Query the Guide to PHARMACOLOGY database (GtoPdb) using natural language or a direct endpoint. |
| `region_to_ccre_screen` | Given starting and ending coordinates, this function retrieves information of intersecting cCREs. |
| `get_genes_near_ccre` | Given a cCRE (Candidate cis-Regulatory Element), this function returns a string containing the |
| `query_remap` | Query the ReMap database for regulatory elements and transcription factor binding sites. |
| `query_mpd` | Query the Mouse Phenome Database (MPD) for mouse strain phenotype data. |
| `query_emdb` | Query the Electron Microscopy Data Bank (EMDB) for 3D macromolecular structures. |
| `query_synapse` | Query Synapse REST API for biomedical datasets and files. |
| `query_pubchem` | Query the PubChem PUG-REST API using natural language or a direct endpoint. |
| `query_chembl` | Query the ChEMBL REST API using natural language, direct endpoint, or specific identifiers. |
| `query_unichem` | Query the UniChem 2.0 REST API using natural language or a direct endpoint. |
| `query_clinicaltrials` | Query the ClinicalTrials.gov API v2 using natural language or a direct endpoint. |
| `query_dailymed` | Query the DailyMed RESTful API using natural language or a direct endpoint. |
| `query_quickgo` | Query the QuickGO API using natural language or a direct endpoint. |
| `query_encode` | Query the ENCODE Portal API to help users locate functional genomics data. |

## Genetics (`genetics.py`)

| Function | Description |
| --- | --- |
| `liftover_coordinates` | Perform liftover of genomic coordinates between hg19 and hg38 formats with detailed intermediate ... |
| `bayesian_finemapping_with_deep_vi` | Performs Bayesian fine-mapping from GWAS summary statistics using deep variational inference. |
| `analyze_cas9_mutation_outcomes` | Analyzes and categorizes mutations induced by Cas9 at target sites. |
| `analyze_crispr_genome_editing` | Analyzes CRISPR-Cas9 genome editing results by comparing original and edited sequences. |
| `simulate_demographic_history` | Simulate DNA sequences with specified demographic and coalescent histories using msprime. |
| `identify_transcription_factor_binding_sites` | Identifies binding sites for a specific transcription factor in a genomic sequence. |
| `fit_genomic_prediction_model` | Fit a linear mixed model for genomic prediction using genotype and phenotype data. |
| `perform_pcr_and_gel_electrophoresis` | Performs PCR amplification of a target transgene and visualizes results using agarose gel electro... |
| `analyze_protein_phylogeny` | Perform phylogenetic analysis on a set of protein sequences. |

## Genomics (`genomics.py`)

| Function | Description |
| --- | --- |
| `unsupervised_celltype_transfer_between_scRNA_datasets` | - |
| `interspecies_gene_conversion` | Convert ENSEMBL gene IDs between different species using BioMart homology mapping. |
| `_fetch_isoform_sequences` | this is not registered as a tool |
| `generate_gene_embeddings_with_ESM_models` | Generate average protein embeddings for a list of Ensembl gene IDs. |
| `annotate_celltype_scRNA` | Annotate cell types based on gene markers and transferred labels using LLM. |
| `annotate_celltype_with_panhumanpy` | Perform hierarchical cell type annotation using panhumanpy and Azimuth Neural Network. |
| `create_scvi_embeddings_scRNA` | - |
| `create_harmony_embeddings_scRNA` | - |
| `get_uce_embeddings_scRNA` | The UCE embeddings are usually our default tools to get cell embeddings, we map UCE embeddings to... |
| `map_to_ima_interpret_scRNA` | Map cell embeddings from the input dataset to the Integrated Megascale Atlas reference dataset us... |
| `get_rna_seq_archs4` | Given a gene name, this function returns the steps it performs and the max K transcripts-per-mill... |
| `get_gene_set_enrichment_analysis_supported_database_list` | - |
| `gene_set_enrichment_analysis` | Perform enrichment analysis for a list of genes, with optional background gene set and plotting f... |
| `analyze_chromatin_interactions` | Analyze chromatin interactions from Hi-C data to identify enhancer-promoter interactions and TADs. |
| `analyze_comparative_genomics_and_haplotypes` | Perform comparative genomics and haplotype analysis on multiple genome samples. |
| `perform_chipseq_peak_calling_with_macs2` | Perform ChIP-seq peak calling using MACS2 to identify genomic regions with significant binding. |
| `find_enriched_motifs_with_homer` | Find DNA sequence motifs enriched in genomic regions using the HOMER motif discovery software. |
| `analyze_genomic_region_overlap` | Analyze overlaps between two or more sets of genomic regions. |
| `generate_embeddings_with_state` | Generate State embeddings for single-cell RNA-seq data using the SE-600M model. |
| `generate_transcriptformer_embeddings` | Generate Transcriptformer embeddings for single-cell RNA-seq data. |

## Glycoengineering (`glycoengineering.py`)

| Function | Description |
| --- | --- |
| `find_n_glycosylation_motifs` | Scan a protein sequence for N-linked glycosylation sequons (N-X-[S/T]). |
| `predict_o_glycosylation_hotspots` | Heuristic O-glycosylation hotspot scoring. |
| `list_glycoengineering_resources` | Curate and summarize external glycoengineering tools and resources. |

## Immunology (`immunology.py`)

| Function | Description |
| --- | --- |
| `analyze_atac_seq_differential_accessibility` | Perform ATAC-seq peak calling and differential accessibility analysis using MACS2. |
| `analyze_bacterial_growth_curve` | Analyzes bacterial growth curve data to determine growth parameters. |
| `isolate_purify_immune_cells` | Simulates the isolation and purification of immune cells from tissue samples. |
| `estimate_cell_cycle_phase_durations` | Estimate cell cycle phase durations using dual-nucleoside pulse labeling data and mathematical mo... |
| `track_immune_cells_under_flow` | Track immune cells under flow conditions and classify their behaviors. |
| `analyze_cfse_cell_proliferation` | Analyze CFSE-labeled cell samples to quantify cell division and proliferation. |
| `analyze_cytokine_production_in_cd4_tcells` | Analyze cytokine production (IFN-γ, IL-17) in CD4+ T cells after antigen stimulation. |
| `analyze_ebv_antibody_titers` | Analyze ELISA data to quantify EBV antibody titers in plasma/serum samples. |
| `analyze_cns_lesion_histology` | Analyzes histological images of CNS lesions to quantify immune cell infiltration, |
| `analyze_immunohistochemistry_image` | Analyzes immunohistochemistry images to quantify protein expression and spatial distribution. |

## Lab Automation (`lab_automation.py`)

| Function | Description |
| --- | --- |
| `_load_pylabrobot_tutorial_content` | Load PLR tutorial/docs text from multiple sources with graceful fallback. |
| `_collect_docs_from_github_zip` | Download GitHub repo zip and extract user_guide docs for the section. |
| `_format_liquid_user_guide` | Assemble liquid-handling docs into a curated order with headings. |
| `get_pylabrobot_documentation_liquid` | Get the documentation for a specific section of the PyLabRobot tutorial. |
| `get_pylabrobot_documentation_material` | - |
| `test_pylabrobot_script` | Test a PyLabRobot script using simulation and validation. |
| `_validate_pylabrobot_imports` | Validate that all PyLabRobot imports in the script are available. |
| `_modify_script_for_testing` | Modify script to use simulation backends and enable tracking. |
| `_execute_script_safely` | Execute the modified script in a safe environment. |
| `_run_script_with_monitoring` | Run the script and monitor its execution. |
| `_create_test_result` | Create the final test result dictionary. |

## Literature (`literature.py`)

| Function | Description |
| --- | --- |
| `fetch_supplementary_info_from_doi` | Fetches supplementary information for a paper given its DOI and returns a research log. |
| `query_arxiv` | Query arXiv for papers based on the provided search query. |
| `query_scholar` | Query Google Scholar for papers based on the provided search query. |
| `query_pubmed` | Query PubMed for papers based on the provided search query. |
| `search_google` | Search using Google search. |
| `advanced_web_search_claude` | Initiate an advanced web search by launching a specialized agent to collect relevant information ... |
| `extract_url_content` | Extract the text content of a webpage using requests and BeautifulSoup. |
| `extract_pdf_content` | Extract the text content of a PDF file given its URL. |

## Microbiology (`microbiology.py`)

| Function | Description |
| --- | --- |
| `optimize_anaerobic_digestion_process` | Optimize anaerobic digestion process conditions to maximize VFA production or methane yield. |
| `analyze_arsenic_speciation_hplc_icpms` | Analyzes arsenic speciation in liquid samples using HPLC-ICP-MS technique. |
| `count_bacterial_colonies` | Count bacterial colonies from an image of agar plate using computer vision techniques. |
| `annotate_bacterial_genome` | Annotate a bacterial genome using Prokka to identify genes, proteins, and functional features. |
| `enumerate_bacterial_cfu_by_serial_dilution` | Quantify bacterial concentration (CFU/mL) using serial dilutions and spot plating. |
| `model_bacterial_growth_dynamics` | Model bacterial population dynamics over time using ordinary differential equations. |
| `quantify_biofilm_biomass_crystal_violet` | Quantifies biofilm biomass using crystal violet staining assay data. |
| `segment_and_analyze_microbial_cells` | Perform automated cell segmentation and quantify morphological metrics from fluorescence microsco... |
| `segment_cells_with_deep_learning` | Perform cell segmentation on fluorescence microscopy images using deep learning. |
| `simulate_generalized_lotka_volterra_dynamics` | Simulate microbial community dynamics using the Generalized Lotka-Volterra (gLV) model. |
| `predict_rna_secondary_structure` | Predict the secondary structure of an RNA molecule using ViennaRNA. |
| `simulate_microbial_population_dynamics` | Performs stochastic simulation of microbial population dynamics using the Gillespie algorithm. |

## Molecular Biology (`molecular_biology.py`)

| Function | Description |
| --- | --- |
| `annotate_open_reading_frames` | Find all Open Reading Frames (ORFs) in a DNA sequence using Biopython. |
| `annotate_plasmid` | Annotate a DNA sequence using pLannotate's command-line interface. |
| `get_gene_coding_sequence` | Retrieves the coding sequence(s) of a specified gene from NCBI Entrez. |
| `get_plasmid_sequence` | Unified function to retrieve plasmid sequences from either Addgene or NCBI. |
| `align_sequences` | Align short sequences (primers) to a longer sequence, allowing for one mismatch. |
| `pcr_simple` | Simulate PCR amplification with given primers and sequence. |
| `digest_sequence` | Simulates restriction enzyme digestion using Biopython's catalyze method and returns the resultin... |
| `find_restriction_sites` | Identifies restriction enzyme sites in a given DNA sequence for specified enzymes. |
| `find_restriction_enzymes` | Finds common restriction enzyme sites in a DNA sequence. |
| `find_sequence_mutations` | Compare query sequence against reference sequence to identify mutations. |
| `design_knockout_sgrna` | Design sgRNAs for CRISPR knockout by searching pre-computed sgRNA libraries. |
| `get_oligo_annealing_protocol` | Return a standard protocol for annealing oligonucleotides without phosphorylation. |
| `get_golden_gate_assembly_protocol` | Return a customized protocol for Golden Gate assembly based on the number of inserts |
| `get_bacterial_transformation_protocol` | Return a standard protocol for bacterial transformation. |
| `design_primer` | Design a single primer within the given sequence window. |
| `design_verification_primers` | Design Sanger sequencing primers to verify a specific region in a plasmid. |
| `design_golden_gate_oligos` | Design oligos for Golden Gate assembly by identifying backbone overhangs |
| `golden_gate_assembly` | Simulate Golden Gate assembly to predict final construct sequences. |

## Pathology (`pathology.py`)

| Function | Description |
| --- | --- |
| `analyze_aortic_diameter_and_geometry` | Analyze aortic diameter and geometry from cardiovascular imaging data. |
| `analyze_atp_luminescence_assay` | Analyze luminescence-based ATP assay data to determine intracellular ATP concentration. |
| `analyze_thrombus_histology` | Analyze histological images of thrombus samples stained with H&E to identify and quantify |
| `analyze_intracellular_calcium_with_rhod2` | Analyzes intracellular calcium concentration using Rhod-2 fluorescent indicator from microscopy i... |
| `quantify_corneal_nerve_fibers` | Quantify the volume/density of immunofluorescence-labeled corneal nerve fibers. |
| `segment_and_quantify_cells_in_multiplexed_images` | Segment cells and quantify protein expression levels from multichannel tissue images. |
| `analyze_bone_microct_morphometry` | Analyze bone microarchitecture parameters from 3D micro-CT images. |

## Pharmacology (`pharmacology.py`)

| Function | Description |
| --- | --- |
| `run_diffdock_with_smiles` | - |
| `docking_autodock_vina` | - |
| `run_autosite` | - |
| `retrieve_topk_repurposing_drugs_from_disease_txgnn` | This function computes TxGNN model predictions for drug repurposing. It takes in the paths to the... |
| `predict_admet_properties` | - |
| `predict_binding_affinity_protein_1d_sequence` | - |
| `analyze_accelerated_stability_of_pharmaceutical_formulations` | Analyzes the stability of pharmaceutical formulations under accelerated storage conditions. |
| `run_3d_chondrogenic_aggregate_assay` | Generates a detailed protocol for performing a 3D chondrogenic aggregate culture assay to evaluat... |
| `grade_adverse_events_using_vcog_ctcae` | Grade and monitor adverse events in animal studies using the VCOG-CTCAE standard. |
| `analyze_radiolabeled_antibody_biodistribution` | Analyze biodistribution and pharmacokinetic profile of radiolabeled antibodies. |
| `estimate_alpha_particle_radiotherapy_dosimetry` | Estimate radiation absorbed doses to tumor and normal organs for alpha-particle radiotherapeutics. |
| `perform_mwas_cyp2c19_metabolizer_status` | Perform a Methylome-wide Association Study (MWAS) to identify CpG sites significantly associated ... |
| `calculate_physicochemical_properties` | Calculate key physicochemical properties of a drug candidate molecule. |
| `analyze_xenograft_tumor_growth_inhibition` | Analyze tumor growth inhibition in xenograft models across different treatment groups. |
| `analyze_pixel_distribution` | Analyze western blot or DNA electrophoresis images and return pixel distribution statistics. |
| `find_roi_from_image` | Find the ROIs of the bands from the image which is determined by analyze_pixel_distribution funct... |
| `analyze_western_blot` | Performs densitometric analysis of Western blot images to quantify relative protein expression. |
| `_load_ddinter_data` | Load DDInter datasets from pickle files, processing if needed. |
| `_process_ddinter_data_inline` | Process DDInter CSV files into standardized pickle files. |
| `_standardize_drug_name_processing` | Standardize drug names for consistent matching during processing. |
| `_build_drug_registry_inline` | Build comprehensive drug registry from all interactions. |
| `_create_interaction_matrix_inline` | Create interaction matrix for fast lookups using standardized drug names. |
| `_create_name_mapping_inline` | Create drug name to ID mapping for fuzzy matching. |
| `_generate_ddinter_statistics_inline` | Generate statistics about the processed data. |
| `_standardize_drug_name` | Standardize drug names using fuzzy matching against DDInter database. |
| `_format_interaction_result` | Format interaction results for research log. |
| `query_drug_interactions` | Query drug-drug interactions from DDInter database. |
| `check_drug_combination_safety` | Analyze safety of a drug combination for potential interactions. |
| `analyze_interaction_mechanisms` | Analyze interaction mechanisms between two specific drugs. |
| `find_alternative_drugs_ddinter` | Find alternative drugs that don't interact with contraindicated drugs. |
| `_standardize_drug_name_fda` | Standardize drug names for FDA API queries. |
| `_apply_fda_filters` | Apply post-query filtering to FDA responses. |
| `_extract_fda_safety_signals` | Extract safety signals from adverse event data. |
| `_generate_fda_statistics` | Generate summary statistics from FDA responses. |
| `_format_adverse_event_summary` | Format adverse event data into readable summary. |
| `_format_drug_label_summary` | Format drug label information into readable summary. |
| `_format_recall_summary` | Format recall information into structured output. |
| `_format_safety_signal_summary` | Format safety signal analysis results. |
| `query_fda_adverse_events` | Query FDA adverse event reports for specific drugs. |
| `get_fda_drug_label_info` | Retrieve FDA drug label information. |
| `check_fda_drug_recalls` | Check for FDA drug recalls and enforcement actions. |
| `analyze_fda_safety_signals` | Analyze safety signals across multiple drugs. |

## Physiology (`physiology.py`)

| Function | Description |
| --- | --- |
| `reconstruct_3d_face_from_mri` | Generate a 3D model of facial anatomy from MRI scans of the head and neck. |
| `analyze_abr_waveform_p1_metrics` | Extracts P1 amplitude and latency from Auditory Brainstem Response (ABR) waveform data. |
| `analyze_ciliary_beat_frequency` | Analyze ciliary beat frequency from high-speed video microscopy data using FFT analysis. |
| `analyze_protein_colocalization` | Analyze colocalization between two fluorescently labeled proteins in microscopy images. |
| `perform_cosinor_analysis` | Performs cosinor analysis on physiological time series data to characterize circadian rhythms. |
| `calculate_brain_adc_map` | Calculate Apparent Diffusion Coefficient (ADC) map from diffusion-weighted MRI data. |
| `analyze_endolysosomal_calcium_dynamics` | Analyze calcium dynamics in endo-lysosomal compartments using ELGA/ELGA1 probe data. |
| `analyze_fatty_acid_composition_by_gc` | Analyzes fatty acid composition in tissue samples using gas chromatography data. |
| `analyze_hemodynamic_data` | Analyzes raw blood pressure data to calculate key hemodynamic parameters. |
| `simulate_thyroid_hormone_pharmacokinetics` | Simulates the transport and binding of thyroid hormones across different tissue compartments |
| `quantify_amyloid_beta_plaques` | - |

## Protocols (`protocols.py`)

| Function | Description |
| --- | --- |
| `search_protocols` | Search protocols.io for relevant protocols based on a natural language query. |
| `get_protocol_details` | Get detailed information about a specific protocol by ID. |
| `_get_protocols_directory` | Get the path to the local protocols directory. |
| `list_local_protocols` | List available protocol files in the local biomni/tool/protocols/ directory. |
| `read_local_protocol` | Read the contents of a local protocol file from biomni/tool/protocols/. |

## Support Tools (`support_tools.py`)

| Function | Description |
| --- | --- |
| `run_python_repl` | Executes the provided Python command in a persistent environment and returns the output. |
| `_capture_matplotlib_plots` | Capture any matplotlib plots that might have been generated during execution. |
| `_apply_matplotlib_patches` | Apply simple monkey patches to matplotlib functions to automatically capture plots. |
| `get_captured_plots` | Get all captured matplotlib plots. |
| `clear_captured_plots` | Clear all captured matplotlib plots. |
| `read_function_source_code` | Read the source code of a function from any module path. |
| `download_synapse_data` | Download data from Synapse using entity IDs. |

## Synthetic Biology (`synthetic_biology.py`)

| Function | Description |
| --- | --- |
| `engineer_bacterial_genome_for_therapeutic_delivery` | Engineer a bacterial genome by integrating therapeutic genetic parts for therapeutic delivery. |
| `analyze_bacterial_growth_rate` | Analyze bacterial growth data and extract growth parameters from OD600 measurements. |
| `analyze_barcode_sequencing_data` | Analyze sequencing data to extract, quantify and determine lineage relationships of barcodes. |
| `analyze_bifurcation_diagram` | Performs bifurcation analysis on a dynamical system and generates a bifurcation diagram. |
| `create_biochemical_network_sbml_model` | Generate a mathematical model of a biochemical network in SBML format. |
| `optimize_codons_for_heterologous_expression` | Analyzes and optimizes a DNA/RNA sequence for improved expression in a heterologous host organism. |
| `simulate_gene_circuit_with_growth_feedback` | Simulate gene regulatory circuit dynamics with growth feedback. |
| `identify_fas_functional_domains` | Identifies functional domains within a Fatty Acid Synthase (FAS) sequence and predicts their roles. |

## Systems Biology (`systems_biology.py`)

| Function | Description |
| --- | --- |
| `query_chatnt` | Call ChatNT to answer a question about a DNA sequence. |
| `perform_flux_balance_analysis` | Perform Flux Balance Analysis (FBA) on a genome-scale metabolic network model. |
| `model_protein_dimerization_network` | Model protein dimerization networks to find equilibrium concentrations of dimers. |
| `simulate_metabolic_network_perturbation` | Construct and simulate kinetic models of metabolic networks and analyze their responses to pertur... |
| `simulate_protein_signaling_network` | Simulate protein signaling network dynamics using ODE-based logic modeling with normalized Hill f... |
| `compare_protein_structures` | Compares two protein structures to identify structural differences and conformational changes. |
| `simulate_renin_angiotensin_system_dynamics` | Simulate the time-dependent concentrations of renin-angiotensin system (RAS) components. |

## Tool Registry (`tool_registry.py`)

| Function | Description |
| --- | --- |
| *No top-level functions found* | - |
