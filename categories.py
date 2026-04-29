#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functional Category Assigner

Description:
    This script categorizes functional annotations of BLAST hits into 
    predefined biological categories using regex-based keyword matching. 
    The categories are prioritized from most specific to most general. 
    After categorizing individual hits, it aggregates the data to determine 
    the "Majority Category" for each query (spacer). If a spacer has valid 
    categories alongside "No annotation", the latter is ignored to favor 
    biological information.

Author: Andrea Fernández
Date: April 2026
"""

import pandas as pd
from collections import Counter
import re

# --- CONFIGURATION ---
input_file = "*.tsv"
output_detailed = "*.tsv"
output_summary = "*.tsv"
# ---------------------

print(f"Loading file '{input_file}'...")
df = pd.read_csv(input_file, sep="\t", low_memory=False)

query_col = df.columns[0]
anno_col = df.columns[-1]

print(f"Using '{query_col}' as Spacer ID and '{anno_col}' as Annotation.")

# =============================================================================
# CATEGORY DICTIONARY v2 — PLASMIDS + CHROMOSOME
# Order matters: from most specific to most general.
# Top categories have priority over the ones below.
# Based on COG classification + experience with NCBI/UniProt annotations.
# =============================================================================

categories_dict = {

# -------------------------------------------------------------------------
# 1. NO ANNOTATION
# -------------------------------------------------------------------------
"No annotation": [
    r"intergenic", r"error", r"not found", r"^none$",
    r"no hits found", r"no significant", r"^na$", r"^n/a$",
    r"no match", r"no result",
],

# -------------------------------------------------------------------------
# 2. PHAGES / VIRAL ELEMENTS
# -------------------------------------------------------------------------
"Phage / Viral": [
    r"phage", r"capsid", r"tail fiber", r"tail spike", r"tail sheath",
    r"tail tube", r"tail base", r"portal", r"terminase", r"baseplate",
    r"holin", r"endolysin", r"lysin", r"prophage", r"tape measure",
    r"head.tail", r"major capsid", r"minor capsid",
    r"antirestriction", r"anti.restriction",
    r"phage integrase", r"phage repressor", r"phage antirepressor",
    r"phage shock protein",
    r"gpv\b", r"gpw\b", r"gpe\b", r"gpd\b",
    r"bacteriophage", r"lambdoid", r"lambda phage",
    r"siphoviridae", r"myoviridae", r"podoviridae",
    r"temperate phage", r"lytic", r"lysogenic",
    r"phage structural", r"phage assembly",
    r"dna circularization protein",
    r"virion", r"coat protein", r"scaffolding protein",
    r"siphovirus", 
],

# -------------------------------------------------------------------------
# 3. ANTIBIOTIC AND HEAVY METAL RESISTANCE
# -------------------------------------------------------------------------
"Resistance": [
    # Beta-lactams
    r"antibiotic resistance", r"antimicrobial resistance",
    r"lactamase", r"beta.lactamase", r"carbapenemase", r"oxacillinase",
    r"metallo.beta.lactamase", r"\bblaz\b", r"\btem\b", r"\bshv\b",
    r"\bctx.m\b", r"\bkpc\b", r"\bndm\b", r"\boxa\b",
    r"\bvim\b", r"\bimp\b",
    # Efflux pumps
    r"efflux pump", r"efflux transporter", r"\befflux\b",
    r"\bacrb\b", r"\bacra\b", r"\bacrz\b", r"\btolc\b",
    r"\bmexb\b", r"\bmexa\b", r"\bmexd\b",
    r"\bomar\b", r"\bomes\b",
    r"resistance.nodulation.division", r"\brnd\b",
    r"major facilitator superfamily.*resistance",
    # Tetracyclines
    r"tetracycline resistance", r"\btet[a-z]\b",
    # Aminoglycosides
    r"aminoglycoside", r"\baac\b", r"\baph\b", r"\bant\b",
    r"aminoglycoside.*transferase",
    # Chloramphenicol
    r"chloramphenicol acetyltransferase", r"\bcat\b.*acetyltransferase",
    r"chloramphenicol resistance",
    # Sulfonamides / Trimethoprim
    r"dihydropteroate synthase", r"\bsul[0-9]\b",
    r"dihydrofolate reductase.*resistance", r"\bdfr[a-z]\b",
    # Quinolones
    r"quinolone resistance", r"\bqnr[a-z]\b", r"gyrase.*resistance",
    # Polymyxins
    r"\bmcr\b", r"colistin resistance", r"polymyxin resistance",
    # Glycopeptides
    r"vancomycin resistance", r"\bvan[a-z]\b",
    # Fosfomycin / Rifampicin
    r"fosfomycin resistance", r"\bfos[a-z]\b",
    r"rifampin resistance", r"rifampicin resistance", r"\brifb\b",
    # Multidrug resistance
    r"multidrug resistance", r"multidrug transporter",
    r"\bmdr\b", r"\bqac\b", r"biocide resistance", r"disinfectant resistance",
    # Heavy metals
    r"heavy metal", r"metal resistance", r"metal tolerance",
    r"mercury resistance", r"\bmer[a-z]\b",
    r"arsenic resistance", r"\bars[a-z]\b",
    r"copper resistance", r"\bcop[a-z]\b", r"\bcus[a-z]\b",
    r"silver resistance", r"\bsil[a-z]\b",
    r"cadmium resistance", r"\bcad[a-z]\b",
    r"zinc resistance", r"\bznt[a-z]\b",
    r"chromate resistance", r"\bchr[a-z]\b",
    r"cation diffusion facilitator",
    r"tellurite resistance",
    r"lead resistance", r"nickel resistance",
    # Others
    r"macrolide resistance", r"lincosamide resistance",
    r"streptogramin resistance", r"\berm[a-z]\b",
    r"oxazolidinone resistance", r"\bcfr\b",
    r"tigecycline resistance",
    r"trimethoprim resistance",
    r"bleomycin resistance",
],

# -------------------------------------------------------------------------
# 4. VIRULENCE & COMPETITION
#    Includes toxin-antitoxin and restriction-modification systems
# -------------------------------------------------------------------------
"Virulence & Competition": [
    # Secretion systems
    r"type iii secretion", r"\bt3ss\b",
    r"type vi secretion", r"\bt6ss\b",
    r"\bvgrg\b", r"\bhcp\b", r"\btae\b", r"\btai\b",
    r"type ii secretion", r"\bt2ss\b",
    r"type iv pili", r"\bt4p\b",
    # Bacteriocins
    r"bacteriocin", r"colicin", r"microcin", r"klebicin",
    r"pyocin", r"bacteriocin immunity",
    # Toxins
    r"hemolysin", r"cytolysin", r"leukotoxin",
    r"cytotoxin", r"enterotoxin", r"shiga toxin", r"verotoxin",
    r"heat.labile toxin", r"heat.stable toxin",
    r"rtx toxin", r"\brtxa\b", r"\brtxb\b",
    r"pore.forming toxin",
    r"virulence factor", r"virulence.associated",
    r"effector protein", r"secreted effector",
    # Iron acquisition
    r"iron acquisition", r"siderophore",
    r"aerobactin", r"yersiniabactin", r"enterobactin", r"catecholate",
    r"\bfep[a-z]\b", r"\bfec[a-z]\b", r"\bfhu[a-z]\b",
    r"\biuc[a-z]\b", r"\biut[a-z]\b",
    r"ferric uptake", r"iron transport", r"\bfur\b.*virulence",
    # Serum resistance
    r"serum resistance", r"complement resistance",
    r"outer membrane protein.*virulence", r"invasin", r"internalin",
    # Toxin-antitoxin systems (TA systems)
    r"antitoxin", r"toxin.antitoxin",
    r"\bvapc\b", r"\bvapb\b",
    r"\brele\b", r"\bpare\b",
    r"\bmazf\b", r"\bmaze\b",
    r"\bccdb\b", r"\bccda\b",
    r"\bpemk\b", r"\bpemi\b",
    r"\bdoc\b", r"\bphd\b",
    r"\bhipa\b", r"\bhipb\b",
    r"\bdinj\b", r"\byafq\b",
    r"\bchpb\b", r"\bchps\b",
    r"\bmqsr\b", r"\bmqsa\b",
    r"\byoeb\b", r"\byefm\b",
    r"\brelb\b",
    r"\bdinq\b", r"\bsrnb\b",
    r"type i toxin", r"type ii toxin", r"type iii toxin",
    r"type i antitoxin", r"type ii antitoxin",
    r"rele.pare family", r"pare family toxin",
    r"hok.sok", r"\bhok\b", r"\bsok\b",
    r"\babce\b", r"\babce2\b",
    r"growth inhibitor", r"bacterial toxin",
    # Restriction-modification
    r"restriction.modification",
    r"type i restriction", r"type ii restriction", r"type iii restriction",
    r"restriction endonuclease", r"methyltransferase.*restriction",
    r"restriction enzyme",
    r"\bhsds\b", r"\bhsdm\b", r"\bhsdr\b",
    r"\bdcm\b methyltransferase", r"\bdam\b methyltransferase",
    r"dna methylase",
    r"dna methyltransferase",
    r"\bdam\b",          # Dam methylase (adenine)
    r"\bdcm\b",          # Dcm methylase (cytosine)
    r"n.6.*methylase",
    r"n.4.*methylase",
    r"n.6.*methyltransferase",
    r"adenine methylase",
    r"cytosine methylase",
],

# -------------------------------------------------------------------------
# 5. CRISPR & BACTERIAL DEFENSE SYSTEMS
# -------------------------------------------------------------------------
"CRISPR & Defense": [
    r"crispr", r"\bcas[0-9]+\b", r"\bcas[a-z]\b",
    r"\bcas1\b", r"\bcas2\b", r"\bcas3\b", r"\bcas4\b",
    r"\bcas5\b", r"\bcas6\b", r"\bcas7\b", r"\bcas8\b",
    r"\bcas9\b", r"\bcas10\b", r"\bcas12\b", r"\bcas13\b",
    r"crispr.associated", r"crispr array", r"crispr spacer",
    r"crispr repeat",
    r"defense system", r"anti.phage", r"antiphage",
    r"thoeris", r"zorya", r"gabija", r"druantia",
    r"\bbrex\b", r"bacteriophage exclusion",
    r"\bpycsar\b", r"\bpassa\b",
    r"abortive infection", r"\babi\b",
    r"retron", r"retron defense",
    r"dnd system", r"phosphorothioation",
    r"\bpycsar\b",
    r"defense island", r"anti.defense",
    r"\bpde\b.*defense",
    r"innate immunity",
    r"anti.crispr",
],

# -------------------------------------------------------------------------
# 6. CONJUGATION & TRANSFER
# -------------------------------------------------------------------------
"Conjugation & Transfer": [
    r"conjugal", r"conjugat",
    r"type iv secretion", r"\bt4ss\b",
    r"\bvirb[0-9]+\b", r"\bvird[0-9]+\b",
    r"\btra[a-z]\b", r"transfer protein",
    r"\bmob[a-z]\b", r"mobilization protein",
    r"relaxase", r"relaxosome",
    r"coupling protein",
    r"mating pair formation",
    r"\bpil[a-z]\b", r"\bpilus\b", r"\bpili\b",
    r"conjugative pilus",
    r"entry exclusion",
    r"\bfino\b", r"\bfinp\b", r"fertility inhibition",
    r"\borit\b", r"origin of transfer",
    r"nick site", r"nic site",
    r"\btrwb\b", r"\btrwc\b",
    r"plasmid transfer",
    r"plasmid sos inhibition",
],

# -------------------------------------------------------------------------
# 7. MOBILE ELEMENTS
# -------------------------------------------------------------------------
"Mobile element": [
    r"transpos", r"transposase",
    r"insertion sequence", r"\bis[0-9]", r"\bisq\b",
    r"\btnp[a-z]\b", r"\btnpb\b",
    r"integrase", r"site.specific integrase",
    r"recombinase", r"tyrosine recombinase", r"serine recombinase",
    r"resolvase", r"invertase",
    r"mobile element",
    r"group ii intron", r"group i intron",
    r"hnh endonuclease",
    r"retron.*transpos",
    r"iscr",
    r"class 1 integron", r"class 2 integron", r"integron integrase",
    r"genomic island", r"pathogenicity island",
    r"integrative.*conjugative", r"\bice\b.*element",
    r"tn[0-9]", r"tn[a-z][0-9]",
    r"passenger gene",
],

# -------------------------------------------------------------------------
# 8. ADHESION
# -------------------------------------------------------------------------
"Adhesion": [
    r"adhesin", r"fimbria", r"fimbriae",
    r"\bfha\b", r"filamentous hemagglutinin",
    r"\bag43\b", r"antigen 43",
    r"hemagglutinin", r"autotransporter",
    r"afimbrial adhesin", r"\bafa\b",
    r"curli", r"\bcsg[a-z]\b",
    r"type i fimbriae", r"p fimbriae", r"\bpap\b",
    r"s fimbriae", r"\bsfa\b",
    r"outer membrane autotransporter",
    r"trimeric autotransporter",
    r"biofilm", r"biofilm formation",
    r"\bbss[a-z]\b", r"\bpga[a-z]\b",
    r"\bcsg[a-z]\b", r"cellulose synthase",
],

# -------------------------------------------------------------------------
# 9. DNA REPLICATION, MAINTENANCE & REPAIR
# -------------------------------------------------------------------------
"DNA Replication, Maintenance and Repair": [
    # Plasmids
    r"\brep[a-z]?\b", r"replication initiation",
    r"replication protein", r"rep protein",
    r"initiator protein",
    r"\binc[a-z]\b", r"incompatibility",
    r"helicase", r"primase",
    r"\bssb\b", r"single.stranded dna binding",
    r"replicon", r"copy number",
    r"\bpar[a-z]\b", r"partition", r"segregation",
    r"chromosome segregation", r"plasmid maintenance",
    r"post.segregational", r"post.segregational",
    r"stability protein", r"plasmid stability",
    r"\bstb[a-z]\b", r"\bsop[a-z]\b",
    r"centromere.binding",
    # Chromosome replication
    r"dna polymerase [iii]", r"dna polymerase i\b",
    r"dna polymerase ii\b",
    r"\bpol[abi]\b", r"\bdnaa\b", r"\bdnab\b", r"\bdnac\b",
    r"\bdnae\b", r"\bdnag\b", r"\bdnan\b", r"\bdnaq\b",
    r"\bdnat\b", r"\bdnax\b", r"\bdnaz\b",
    r"replicative helicase", r"clamp loader", r"sliding clamp",
    r"okazaki fragment",
    r"topoisomerase", r"gyrase", r"topoisomerase i", r"topoisomerase iv",
    r"\bgyrb\b", r"\bgyra\b", r"\bpara\b.*topoisomerase", r"\bparc\b",
    r"dna ligase",
    r"replication-associated recombination",
    # Repair
    r"dna repair", r"mismatch repair",
    r"\bmutl\b", r"\bmuts\b", r"\bmuth\b",
    r"nucleotide excision repair",
    r"\buvra\b", r"\buvrb\b", r"\buvrc\b", r"\buvrd\b",
    r"base excision repair",
    r"\buth\b", r"\budg\b", r"\bung\b", r"\bfpg\b", r"\bogg\b",
    r"alkylation repair", r"\bada\b", r"\bogt\b", r"\balkb\b",
    r"\brecA\b", r"\breca\b",
    r"sos response protein",
    r"\bdinb\b", r"\bumuC\b", r"\bumud\b",
    r"error.prone polymerase", r"translesion synthesis",
    r"double.strand break", r"recombination repair",
    r"\brecb\b", r"\brecc\b", r"\brecd\b", r"\brecf\b",
    r"\brecg\b", r"\breco\b", r"\brecq\b", r"\brecr\b",
    r"holliday junction", r"branch migration",
    r"\brusA\b", r"\brus\b.*junction",
    r"exonuclease.*repair", r"endonuclease.*repair",
    r"\bsbia\b",
    r"xthA", r"exodeoxyribonuclease",
    r"very short patch repair",
    r"\bdcd\b",
    r"\bdgtpase\b",
    r"\bdgt\b",
    r"deoxyguanosinetriphosphate",
    r"uracil-dna glycosylase", r"uracil dna glycosylase",
    # Chromosome organization
    r"chromosome compaction", r"chromosome organization",
    r"nucleoid", r"nucleoid.associated",
    r"dna condensation", r"chromosome segregation",
    r"smc complex", r"\bsmcc\b", r"\bscpab\b",
    r"structural maintenance",
    r"histone.like", r"\bhu\b.*dna", r"\bifhb\b", r"\bfis\b.*nucleoid",
    r"genome integrity",
    r"ter.*replication", r"terminus.*replication",
    r"\bdif\b.*site", r"\bxer[cd]\b",
],

# -------------------------------------------------------------------------
# 11. CELL DIVISION & MORPHOLOGY
# -------------------------------------------------------------------------
"Cell division & morphology": [
    r"cell division",
    r"\bftsz\b", r"\bftsa\b", r"\bftsi\b", r"\bftsq\b",
    r"\bftsw\b", r"\bftsl\b", r"\bftsb\b", r"\bftse\b",
    r"\bftsn\b", r"\bftsk\b", r"\bftsx\b", r"\bftsyy\b",
    r"divisome", r"z.ring", r"septal ring",
    r"septum formation", r"cell constriction",
    r"\bminc\b", r"\bmind\b", r"\bmine\b", r"min system",
    r"nucleoid occlusion", r"\bnoc\b", r"\bslma\b",
    r"\bsulA\b", r"filamentation",
    r"\bmreb\b", r"\bmrec\b", r"\bmred\b", r"\brodz\b",
    r"actin.like", r"tubulin.like",
    r"cell elongation", r"rod shape",
    r"penicillin.binding protein", r"\bpbp[0-9]\b",
    r"\bponb\b", r"\bpona\b",
    r"cell shape", r"cell morphology",
    r"polar landmark",
    r"septation",
],

# -------------------------------------------------------------------------
# 12. CELL WALL & ENVELOPE
# -------------------------------------------------------------------------
"Cell wall & Envelope": [
    # Peptidoglycan
    r"peptidoglycan", r"murein",
    r"\bmura\b", r"\bmurb\b", r"\bmurc\b", r"\bmurd\b",
    r"\bmure\b", r"\bmurf\b", r"\bmurg\b", r"\bmurh\b",
    r"\bmuri\b", r"\bmurj\b",
    r"muropeptide", r"n.acetylmuramic", r"n.acetylglucosamine",
    r"d.alanine", r"d.alanyl", r"d.ala.d.ala",
    r"transpeptidase", r"transglycosylase",
    r"lytic transglycosylase", r"muramidase",
    r"amidase.*cell wall", r"cell wall amidase",
    r"\bltg[a-z]\b",
    # LPS
    r"lipopolysaccharide", r"lps biosynthesis", r"lps assembly",
    r"o.antigen", r"lipid a", r"kdo",
    r"\bwzy\b", r"\bwzx\b", r"\bwcaa\b", r"\brfb\b",
    r"lipooligosaccharide",
    # Outer membrane
    r"outer membrane protein", r"\bomp[a-z]\b",
    r"lipoprotein.*outer", r"\blpp\b", r"\bpal\b",
    r"\bbama\b", r"\bbamb\b", r"\bbamc\b", r"\bbamd\b", r"\bbame\b",
    r"bam complex",
    r"\btolb\b", r"\btolc\b", r"\btold\b", r"\btolq\b", r"\btolr\b",
    r"tol.pal system",
    # Capsule and exopolysaccharides
    r"capsule", r"capsular polysaccharide",
    r"colanic acid", r"\bwca\b",
    r"exopolysaccharide",
    r"cellulose synthase",
    # Teichoic acids (Gram+)
    r"teichoic acid", r"wall teichoic", r"lipoteichoic",
    r"\btar[a-z]\b", r"\bltaf\b",
    # General
    r"cell wall synthesis", r"cell wall biosynthesis",
    r"cell envelope",
    r"outer membrane biogenesis",
    r"membrane integrity",
],

# -------------------------------------------------------------------------
# 13. STRESS & SOS RESPONSE
# -------------------------------------------------------------------------
"Stress response": [
    # Oxidative stress
    r"oxidative stress", r"oxidative damage",
    r"superoxide dismutase", r"\bsod[ab]\b",
    r"catalase", r"\bkata\b", r"\bkatg\b", r"\bkatb\b",
    r"peroxidase", r"glutathione peroxidase", r"alkyl hydroperoxide",
    r"\bahpc\b", r"\bahpf\b",
    r"thioredoxin", r"thioredoxin reductase", r"\btrxa\b", r"\btrxb\b",
    r"glutaredoxin", r"\bgrxa\b",
    r"oxyR", r"\boxyr\b", r"\bsoxr\b", r"\bsoxs\b",
    r"envelope stress",
    # Heat stress / Stress chaperones
    r"heat shock", r"heat stress",
    r"\bgroel\b", r"\bgroes\b", r"\bdnak\b", r"\bdnaj\b",
    r"\bhspa\b", r"\bhsp[0-9]+\b",
    r"\bhtpg\b", r"\bclpb\b",
    r"small heat shock", r"\bshsp\b",
    r"\bsigma.*heat", r"\brpoh\b", r"\bsig32\b",
    # Cold stress
    r"cold shock", r"cold-shock", r"\bcsp[a-z]\b",
    r"cold.inducible", r"low.temperature",
    # Osmotic stress
    r"osmotic stress", r"osmotic shock",
    r"\bosmc\b", r"\bosma\b", r"\bkdp[a-z]\b",
    r"compatible solute", r"osmoprotectant",
    r"betaine", r"ectoine", r"trehalose",
    r"\btrea\b", r"\botsa\b",
    # pH stress
    r"acid resistance", r"acid stress",
    r"acid tolerance", r"ph stress",
    r"\bgad[ab]\b", r"\baro\b.*acid",
    # Universal stress proteins
    r"universal stress protein", r"\busp[a-z]\b",
    r"stationary phase", r"\brpof\b", r"\bsigS\b",
    r"\bdps\b", r"starvation",
    r"\bsigB\b", r"general stress",
    # SOS response
    r"\bsosa\b.*response", r"sos response",
    r"\bdinf\b", r"\bdini\b", r"\bsulb\b",
    # Stringent response
    r"stringent response", r"\bspota\b", r"\brela\b",
    r"(p)ppgpp", r"guanosine tetraphosphate",
],

# -------------------------------------------------------------------------
# 15. MOTILITY & CHEMOTAXIS
# -------------------------------------------------------------------------
"Motility & Chemotaxis": [
    r"flagell", r"motility",
    r"flagellar motor", r"flagellar hook",
    r"\bflga\b", r"\bflgb\b", r"\bflgc\b", r"\bflgd\b",
    r"\bflge\b", r"\bflgf\b", r"\bflgg\b", r"\bflgh\b",
    r"\bflgi\b", r"\bflgj\b", r"\bflgk\b", r"\bflgl\b",
    r"\bflic\b", r"\bflid\b", r"\bflie\b",
    r"\bflha\b", r"\bflhb\b", r"\bflhc\b", r"\bflhd\b",
    r"\bmota\b", r"\bmotb\b",
    r"chemotaxis", r"chemoreceptor",
    r"\bchea\b", r"\bcheb\b", r"\bchec\b",
    r"\bchew\b", r"\bchey\b", r"\bchez\b",
    r"methyl.accepting chemotaxis",
    r"mcpl", r"tar receptor", r"tsr receptor",
    r"gliding motility", r"twitching motility",
    r"type ivb pili.*motility",
],

# -------------------------------------------------------------------------
# 16. PROTEIN SECRETION & TRANSLOCATION
# -------------------------------------------------------------------------
"Protein secretion & translocation": [
    r"sec pathway", r"sec translocon",
    r"\bseca\b", r"\bsecb\b", r"\bsecd\b", r"\bsece\b",
    r"\bsecf\b", r"\bsecg\b", r"\bsecy\b",
    r"signal peptide", r"signal peptidase",
    r"\bspase\b", r"\bsigl\b",
    r"tat pathway", r"twin.arginine",
    r"\btata\b", r"\btatb\b", r"\btatc\b",
    r"type i secretion", r"\bt1ss\b",
    r"type ii secretion", r"\bt2ss\b",
    r"\bgsp[a-z]\b", r"\bxcp[a-z]\b",
    r"type v secretion", r"autotransporter",
    r"signal recognition particle", r"\bsrp\b",
    r"\bffh\b", r"\bftsY\b",
    r"protein export", r"protein translocation",
    r"\byaeg\b", r"\bptsa\b",
],

# -------------------------------------------------------------------------
# 17. TRANSCRIPTION & REGULATION
# -------------------------------------------------------------------------
"Transcription & Regulation": [
    r"transcription factor", r"transcriptional regulator", r"transcription regulator",
    r"transcriptional activator", r"transcriptional repressor",
    r"rna polymerase",
    r"sigma factor", r"sigma.70", r"sigma.32", r"sigma.54",
    r"sigma.28", r"sigma.38", r"sigma.24",
    r"anti.sigma factor",
    r"antitermination", r"antiterminator",
    # Two-component systems
    r"two.component", r"response regulator", r"sensor kinase",
    r"sensor histidine kinase", r"phosphorelay",
    r"\bphop\b", r"\bphoq\b", r"\bphob\b", r"\bphor\b",
    r"\benvz\b", r"\bompR\b",
    r"\bnark\b.*regulator", r"\bnarr\b", r"\bnarl\b",
    # Global transcription factors / nucleoid
    r"helix.turn.helix", r"winged.helix",
    r"lysr.type", r"arac.type", r"gntr.type", r"marr.type",
    r"tetr.type", r"laci.type",
    r"laci family", r"lysr family", r"arac family",
    r"gntr family", r"iclr family", r"xre family",
    r"\bfis\b", r"\bihf\b", r"\bh.ns\b", r"\bhns\b", r"\bhu\b protein",
    r"\bcrp\b", r"\bcap\b.*activator",
    r"\bfnr\b", r"\barc[ab]\b",
    r"transcription termination",
    r"\brho\b.*transcription",
    r"transcription.*\brho\b",
    r"\bnusa\b", r"\bnusb\b", r"\bnusg\b",  # other termination/antitermination factors
    r"termination factor",
    r"rho.dependent",
    r"rho.independent",
    # Quorum sensing
    r"quorum sensing", r"autoinducer",
    r"\bluxr\b", r"\bluxn\b", r"\bluxs\b", r"\bluxI\b",
    r"acyl.homoserine lactone", r"\bahl\b",
    r"diffusible signal factor",
    # General regulators
    r"global regulator",
    r"sos response.*regul", r"\blexa\b",
    r"\bcfp\b", r"\bdgkr\b",
    r"c.di.gmp.*regul", r"cyclic di.gmp.*regul",
],

# -------------------------------------------------------------------------
# 18. TRANSLATION & RIBOSOME
# -------------------------------------------------------------------------
"Translation & Ribosome": [
    r"ribosomal protein", r"ribosome",
    r"\btrna\b", r"\brrna\b",
    r"aminoacyl.trna", r"trna synthetase",
    r"elongation factor", r"initiation factor", r"release factor",
    r"peptidyl transferase", r"\b50s\b", r"\b30s\b",
    r"\b16s\b", r"\b23s\b", r"\b5s\b",
    r"ribosome assembly", r"ribosome maturation",
    r"\bnusa\b", r"\bnusb\b", r"\bnusg\b",
    r"translation factor",
    r"peptide chain release", r"ribosome recycling",
    r"\bera\b.*gtpase", r"\brbfa\b", r"\brbgd\b",
    r"rrna methyltransferase", r"rrna modification",
    r"trna modification", r"trna processing",
    r"\btrnae\b", r"\btrnaf\b",
    r"pseudouridine synthase.*trna",
    r"tmrna", r"\bssra\b.*tmrna",
    r"\bsmpb\b",
    r"translational throttle",
    r"\betta\b",
    r"abc.f.*ribosom",
    r"ribosome.*abc.f",
],

# -------------------------------------------------------------------------
# 19. METABOLISM & TRANSPORT
# -------------------------------------------------------------------------
"Metabolism & Transport": [
    # Amino acids
    r"amino acid.*biosynthesis", r"amino acid.*degradation",
    r"amino acid.*metabolism",
    r"aspartate.*biosynthesis", r"glutamate.*biosynthesis",
    r"arginine biosynthesis", r"\barg[a-h]\b",
    r"lysine biosynthesis", r"lysine degradation",
    r"threonine biosynthesis", r"isoleucine biosynthesis",
    r"valine biosynthesis", r"leucine biosynthesis",
    r"\bilv[a-z]\b", r"\bleu[a-z]\b",
    r"tryptophan biosynthesis", r"\btrp[a-z]\b",
    r"phenylalanine biosynthesis", r"tyrosine biosynthesis",
    r"\bphe[a-z]\b", r"\btyr[a-z]\b",
    r"histidine biosynthesis", r"\bhis[a-z]\b",
    r"serine biosynthesis", r"glycine biosynthesis",
    r"\bser[abc]\b", r"\bgly[a-z]\b",
    r"proline biosynthesis", r"\bpro[abc]\b",
    r"methionine biosynthesis", r"\bmet[a-z]\b",
    r"cysteine biosynthesis", r"\bcys[a-z]\b",
    r"glutamine synthetase", r"\bglna\b",
    r"glutamate dehydrogenase", r"\bgdh\b",
    r"aspartate aminotransferase",
    r"ornithine", r"citrulline", r"urea cycle",
    r"branched.chain amino acid",
    r"aromatic amino acid",
    r"\baro[a-z]\b.*biosynthesis",
    r"aminotransferase", r"transaminase",
    r"symporter",
    
    # Glycolysis / gluconeogenesis
    r"glycolysis", r"gluconeogenesis",
    r"\bpfk\b", r"\bpyk\b", r"\beno\b.*enolase",
    r"phosphoglycerate", r"\bgapdh\b",
    r"phosphofructokinase", r"pyruvate kinase",
    r"enolase", r"phosphoglycerate mutase",
    r"triosephosphate isomerase",
    r"glycogen phosphorilase",
    # TCA Cycle
    r"tca cycle", r"citric acid cycle", r"krebs cycle",
    r"citrate synthase", r"\bacs\b.*acetyl",
    r"isocitrate dehydrogenase", r"\bicd\b",
    r"alpha.ketoglutarate", r"oxoglutarate",
    r"succinate dehydrogenase", r"\bsdh[a-z]\b",
    r"fumarate", r"malate dehydrogenase",
    r"oxaloacetate",
    r"hydrogenase",
    r"aconitate hydratase",
    # Pentose phosphate pathway
    r"pentose phosphate", r"transketolase", r"transaldolase",
    r"glucose.6.phosphate dehydrogenase",
    r"\bzwf\b", r"\bgnd\b",
    # Acetate / propionate metabolism
    r"acetyl.coa", r"acetate kinase",
    r"\bpta\b.*phosphotransacetylase", r"\bacs\b",
    r"propionate", r"methylcitrate",
    # Anaerobic respiration
    r"anaerobic respiration", r"nitrate respiration",
    r"fumarate respiration", r"formate dehydrogenase",
    r"\bfrd[a-z]\b", r"\bnar[a-z]\b",
    r"\bnir[a-z]\b", r"nitrite reductase",
    r"\bnap[a-z]\b",
    # Respiratory chain
    r"electron transport chain", r"respiratory chain",
    r"cytochrome", r"cytochrome oxidase",
    r"cytochrome bc1", r"ubiquinol",
    r"\bndh\b", r"\bnadh dehydrogenase\b",
    r"\bsdh\b", r"succinate ubiquinone",
    r"\bcyo[abcd]\b", r"\bcyd[ab]\b",
    r"atp synthase", r"\batpa\b", r"\batpb\b",
    r"\batpf\b", r"\batph\b",
    # Sugars
    r"sugar.*metabolism", r"glucose.*transport",
    r"phosphotransferase system", r"\bpts\b",
    r"\bpep\b.*phosphotransferase",
    r"lactose.*metabolism", r"\blac[abcz]\b",
    r"maltose.*metabolism", r"\bmal[a-z]\b",
    r"galactose.*metabolism", r"\bgal[a-z]\b",
    r"galactosidase",
    r"arabinose.*metabolism", r"\bara[a-z]\b",
    r"ribose.*metabolism", r"\brib[a-z]\b",
    r"xylose.*metabolism", r"\bxyl[ab]\b",
    r"xylosidase",
    r"gluconate.*metabolism", r"\bgnt[a-z]\b",
    # Fermentation
    r"fermentation", r"alcohol dehydrogenase",
    r"acetaldehyde", r"lactate dehydrogenase",
    r"pyruvate.*fermentation", r"mixed acid fermentation",
    r"\badh[a-z]\b", r"\bldh\b",
    r"acetoin", r"butanediol", r"formate",
    r"phosphoenolpyruvate",
    # Nucleotides
    r"purine biosynthesis", r"\bpur[a-z]\b",
    r"pyrimidine biosynthesis", r"\bpyr[a-z]\b",
    r"inosine", r"adenosine", r"guanosine",
    r"thymidylate synthase", r"\bthya\b",
    r"dihydrofolate reductase", r"\bdhfr\b", r"\bfola\b",
    r"ribonucleotide reductase", r"\bnrda\b", r"\bnrdb\b",
    r"nucleotide kinase", r"adenylate kinase",
    r"\badk\b",
    r"nucleotide.*metabolis",
    r"deoxyribose",
    r"phosphoribosyl",
    r"nucleoside.*hydrolase",
    r"\bdrm\b", r"\bdeo[abcd]\b",
    r"nucleoside phosphorylase",
    r"purine phosphorylase",
    r"purine nucleoside phosphorylase",
    r"pyrimidine phosphorylase",
    r"pyrimidine nucleoside phosphorylase",
    r"thymidine phosphorylase",
    r"uridine phosphorylase",
    r"uridine/thymidine phosphorylase",
    r"adenosine phosphorylase",
    r"guanosine phosphorylase",
    r"inosine phosphorylase",
    r"xanthosine phosphorylase",
    r"deoxyadenosine phosphorylase",
    r"deoxyguanosine phosphorylase",
    r"deoxyinosine phosphorylase",
    r"deoxyuridine phosphorylase",
    r"methylthioadenosine phosphorylase",
    r"polynucleotide phosphorylase",
    r"nucleosidase",
    # Cofactors
    r"coenzyme", r"cobalamin", r"\bcob[a-z]\b",
    r"biotin", r"\bbio[abcdf]\b",
    r"folate", r"tetrahydrofolate", r"\bfol[abcde]\b",
    r"thiamine", r"\bthia\b", r"\bthib\b", r"\bthic\b",
    r"riboflavin", r"\brib[a-z]\b.*flavin",
    r"heme biosynthesis", r"\bheme\b.*synth",
    r"porphyrin", r"\bhem[abcdefgln]\b",
    r"nicotinamide", r"\bnad\b.*biosynthesis", r"\bpnca\b",
    r"\bnad[abce]\b", r"\bnics\b",
    r"pantothenate", r"\bpan[bcd]\b",
    r"lipoic acid", r"\blip[ab]\b",
    r"menaquinone", r"ubiquinone biosynthesis",
    r"\bmen[abcde]\b", r"\bubi[abcefgh]\b",
    r"pyridoxal phosphate", r"pyridoxine", r"\bpdb\b",
    r"molybdenum cofactor", r"\bmog[abc]\b",
    r"iron.sulfur cluster", r"\bisc[absu]\b", r"\bsuf[abcde]\b",
    r"glutathione biosynthesis", r"\bgsh\b.*biosynthesis",
    r"\bgsh[ab]\b",
    # Lipids
    r"fatty acid biosynthesis", r"fatty acid synthesis",
    r"fatty acid.*metabolism",
    r"\bfab[abdefghiz]\b",
    r"acyl carrier protein", r"\bacp\b",
    r"malonyl", r"acetyl.coa carboxylase",
    r"fatty acid oxidation", r"beta.oxidation",
    r"\bfad[abijl]\b",
    r"phospholipid biosynthesis", r"phospholipid metabolism",
    r"\bpls[abc]\b", r"\bpsd\b.*phospholipid", r"\bpss[ab]\b",
    r"phosphatidylethanolamine", r"phosphatidylglycerol",
    r"cardiolipin",
    r"\bcls[ab]\b",
    r"lipid a biosynthesis",
    r"isoprenoid", r"terpenoid",
    r"mevalonate", r"mep.dxp pathway",
    r"\bdxs\b", r"\bisp[a-z]\b",
    r"hopanoid", r"hopene",

    # Other transporters
    r"transporter", r"permease",
    r"abc transporter", r"atp.binding cassette",
    r"\babc\b.*transport",
    r"mfs transporter", r"major facilitator",
    r"pts.*transporter",
    r"channel protein",
    r"porin", r"outer membrane channel",
    r"solute.binding protein", r"substrate.binding protein",
    r"atpase.*transport",
    r"secondary transporter",
    r"ion transporter", r"cation transporter", r"anion transporter",
    r"potassium transporter", r"\bkup\b", r"\bkdp[abcde]\b",
    r"sodium transporter",
    r"phosphate transporter", r"\bpst[abcs]\b",
    r"sulfate transporter", r"\bcys[aptu]\b",
    r"iron transporter", r"iron uptake", r"\bfeo[ab]\b", r"\bfep[abcg]\b",
    r"zinc transporter", r"manganese transporter",
    r"magnesium transporter", r"\bcora\b", r"\bmgt[ab]\b",
    r"calcium transporter",
    r"drug transporter",
    r"sugar transporter", r"amino acid transporter",
    r"dipeptide transporter", r"oligopeptide transporter",
    r"\bapp[abcdf]\b", r"\bopt[abc]\b",
    r"nucleoside transporter",
    r"vitamin transporter",
    r"copper transporter",
    
    # Inorganic and redox
    r"sulfur assimilation", r"sulfur metabolism",
    r"\bcys[jnqr]\b", r"sulfite reductase", r"sulfate adenylyl",
    r"nitrogen assimilation", r"nitrogen metabolism",
    r"nitrogen fixation", r"\bnif[a-z]\b",
    r"nitrogenase", r"dinitrogenase",
    r"ammonia assimilation",
    r"phosphate assimilation", r"\bphob\b.*phosphate",
    r"\bphna\b", r"alkaline phosphatase",
    r"ferredoxin", r"flavodoxin",
    r"thioredoxin.*redox",
    r"periplasmic binding protein.*inorganic",
    
    # Rest of metabolism
    r"kinase", r"dehydrogenase", r"synthase", r"synthetase",
    r"reductase", r"oxygenase", r"oxidase", r"oxidoreductase",
    r"hydrolase", r"transferase", r"isomerase", r"lyase",
    r"phosphatase", r"epimerase", r"mutase", r"racemase",
    r"acyltransferase", r"methyltransferase",
    r"metabolism", r"biosynthesis", r"catabolism",
    r"abc transporter", r"mfs transporter",
    r"atpase", r"atp.binding",
    r"solute.binding", r"substrate.binding",
    r"outer membrane channel",
    r"methylase",
    r"phosphodiester glycosidase",
    r"putrescine",
],

# -------------------------------------------------------------------------
# 26. CHAPERONES & PROTEOLYSIS
# -------------------------------------------------------------------------
"Chaperones & Proteolysis": [
    r"chaperone", r"chaperonin",
    r"\bgroel\b", r"\bgroes\b",
    r"\bdnak\b", r"\bdnaj\b", r"\bgrpe\b",
    r"\bhspa\b", r"\bhsp[0-9]+\b",
    r"\bhtpg\b",
    r"trigger factor", r"\btig\b.*chaperone",
    r"\bsur[a-z]\b.*chaperone", r"\bskp\b.*chaperone",
    r"prefoldin",
    r"protein folding",
    # Proteolysis / proteases
    r"protease", r"peptidase", r"serine protease",
    r"metalloprotease", r"cysteine protease",
    r"\bclp[abcpqrsx]\b",
    r"\blon\b.*protease",
    r"\bhslu\b", r"\bhslv\b",
    r"\bfts[h]\b.*protease",
    r"aaa..*protease",
    r"signal peptidase", r"leader peptidase",
    r"propeptide", r"protein maturation protease",
    r"\bdegp\b", r"\bdegq\b", r"\bdeg[a-z]\b",
    r"\bprc\b.*protease",
    r"periplasmic protease",
    r"intracellular protease",
    r"ubiquitin.like",
    r"pupylation", r"\bpafA\b",
],

# -------------------------------------------------------------------------
# 27. POST-TRANSLATIONAL MODIFICATION
# -------------------------------------------------------------------------
"Post-translational modification": [
    r"protein modification",
    r"phosphorylation", r"acetylation", r"methylation.*protein",
    r"protein kinase", r"serine.threonine kinase",
    r"adenylation", r"glutamylation",
    r"\bpka\b", r"\bhanks\b kinase",
    r"protein phosphatase",
    r"protein acetylase", r"protein methylase",
    r"formyl", r"n.terminal",
    r"signal peptide cleavage",
    r"lipoprotein signal peptidase",
    r"disulfide bond", r"\bdsbA\b", r"\bdsbB\b", r"\bdsbC\b", r"\bdsbD\b",
    r"sulfhydryl oxidase",
    r"protein sumoylation",
    r"selenocysteine",
    r"pyroglutamate",
],

# -------------------------------------------------------------------------
# 28. SIGNALING & 2ND MESSENGERS
# -------------------------------------------------------------------------
"Signalling & 2nd messengers": [
    r"signal transduction",
    r"c.di.gmp", r"cyclic di.gmp",
    r"diguanylate cyclase", r"\bdgc\b",
    r"phosphodiesterase.*c.di.gmp", r"\bpde\b.*signalling",
    r"\bggdef\b", r"\beal\b.*domain",
    r"camp.*signalling", r"cyclic amp",
    r"\badenylate cyclase\b",
    r"second messenger",
    r"serine.threonine.*kinase.*signall",
    r"(p)ppGpp signall",
    r"sensory histidine kinase",
    r"phosphotransfer",
    r"\bche[yz]\b.*signal",
    r"methyl.accepting chemotaxis.*signal",
],

# -------------------------------------------------------------------------
# 29. SPORULATION & DIFFERENTIATION
# -------------------------------------------------------------------------
"Sporulation & Differentiation": [
    r"\bspore", r"\bsporulation", r"endospore",
    r"forespore", r"cortex", r"cortex.*spore", r"spore coat",
    r"\bspo0[a-z]\b", r"\bsig[efgk]\b", r"\bsig[efgk]\b.*sporulation",
    r"\bsign\b.*spore", r"\bkine\b.*spore",
    r"germination",
    r"dipicolinic acid",
    r"aerial mycelium", r"substrate mycelium",
    r"colony morphology",
    r"developmental regulation",
    r"stationary phase.*morphol",
],

# -------------------------------------------------------------------------
# 32. HYPOTHETICAL PROTEINS
# -------------------------------------------------------------------------
"Hypothetical protein": [
    r"hypothetical", r"uncharacterized", r"unknown function",
    r"\bduf[0-9]*\b", r"domain of unknown function",
    r"putative protein", r"putative uncharacterized",
    r"predicted protein",
    r"small protein",
    r"orf[0-9]",
    r"\by[a-z]{3,4}\b",  # Genes like yaaA, ybhQ, etc.
],
}


# =============================================================================
# CLASSIFICATION LOGIC
# =============================================================================

def assign_category(product_name):
    """Assigns a category based on the product name matching predefined regex patterns."""
    if pd.isna(product_name):
        return "No annotation"
    
    prod_lower = str(product_name).lower()
    
    for category, keywords in categories_dict.items():
        for kw in keywords:
            if re.search(kw, prod_lower):
                return category
    
    return "Others"


print("Categorizing each hit...")
df['Functional_Category'] = df[anno_col].apply(assign_category)

df.to_csv(output_detailed, sep="\t", index=False)
print(f"-> Detailed file saved: {output_detailed}")


# =============================================================================
# CALCULATE MAJORITY CATEGORY
# =============================================================================

print("\nCalculating majority category per spacer...")

def get_majority_category(x):
    """Calculates the most common category for a given group, excluding 'No annotation' if possible."""
    counts = Counter(x)
    
    if len(counts) == 1 and "No annotation" in counts:
        return "No annotation"
    
    if "No annotation" in counts and len(counts) > 1:
        del counts["No annotation"]
    
    return counts.most_common(1)[0][0]

summary_df = df.groupby(query_col)['Functional_Category'].agg(get_majority_category).reset_index()
summary_df.columns = ['Query_ID', 'Majority_Category']

summary_df.to_csv(output_summary, sep="\t", index=False)
print(f"-> Summary saved: {output_summary}")
print("\nProcess successfully completed!")