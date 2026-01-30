import re
from typing import Dict, List, Tuple

def classify_gene(gene_name: str) -> Dict[str, str]:
    name = gene_name.strip()
    
    # Ensembl gene ID
    if re.match(r'^ENSG\d{11}$', name):
        return { 
            
            'type': 'Ensembl Gene ID',
            'description': 'Ensembl database gene identifier',
            'example': 'ENSG00000141510 (TP53)',
            'pattern': 'ENSG followed by 11 digits'
        }
    
    # Ensembl transcript ID
    if re.match(r'^ENST\d{11}$', name):
        return { 
            
            'type': 'Ensembl Transcript ID',
            'description': 'Ensembl database transcript identifier',
            'example': 'ENST00000269305',
            'pattern': 'ENST followed by 11 digits'
        }
    
    # LINC (Long Intergenic Non-Coding RNA)
    if re.match(r'^LINC\d{5,}$', name, re.IGNORECASE):
        return { 
            
            'type': 'LINC - Long Intergenic Non-Coding RNA',
            'description': 'Systematic naming for long non-coding RNAs between genes',
            'example': 'LINC00634',
            'pattern': 'LINC followed by 5+ digits'
        }
    
    # MIR (microRNA)
    if re.match(r'^MIR\d+[A-Z]*\d*$', name, re.IGNORECASE):
        return { 
            
            'type': 'MicroRNA (miRNA)',
            'description': 'MicroRNA gene symbol',
            'example': 'MIR21, MIR125A',
            'pattern': 'MIR followed by numbers, optionally letters and more numbers'
        }
    
    # Clone-based identifiers (RP, AC, CT, etc.)
    if re.match(r'^(RP|AC|AL|AP|CR|CT|CU|CH|KB)\d+[-.][\dA-Z]+\.\d+$', name, re.IGNORECASE):
        return { 
            
            'type': 'Clone-based Identifier',
            'description': 'Gene identified by BAC/PAC clone location (Human Genome Project era)',
            'example': 'RP11-458J1.1, AC092580.4',
            'pattern': 'Library prefix + clone ID + gene number (e.g., RP11-458J1.1)'
        }
    
    # LOC genes (uncharacterized locus)
    if re.match(r'^LOC\d+$', name, re.IGNORECASE):
        return { 
            
            'type': 'LOC - Uncharacterized Locus',
            'description': 'Placeholder name for genes without official symbols',
            'example': 'LOC100287792',
            'pattern': 'LOC followed by digits'
        }
    
    # Pseudogenes
    if re.search(r'P\d+$', name) or re.search(r'PS\d+$', name):
        return { 
            
            'type': 'Pseudogene',
            'description': 'Non-functional gene copy',
            'example': 'PTENP1, GAPDHP1',
            'pattern': 'Gene symbol ending in P# or PS#'
        }
    
    # Small nucleolar RNA (snoRNA)
    if re.match(r'^SNOR[AD]\d+[A-Z]*$', name, re.IGNORECASE):
        return { 
            
            'type': 'Small Nucleolar RNA (snoRNA)',
            'description': 'Small RNA involved in rRNA modification',
            'example': 'SNORD116, SNORA74A',
            'pattern': 'SNORD or SNORA followed by numbers'
        }
    
    # Small nuclear RNA (snRNA)
    if re.match(r'^RNU\d+[A-Z]*\d*$', name, re.IGNORECASE):
        return { 
            
            'type': 'Small Nuclear RNA (snRNA)',
            'description': 'Small RNA in splicing and other nuclear processes',
            'example': 'RNU1-1, RNU6-1',
            'pattern': 'RNU followed by numbers'
        }
    
    # Long non-coding RNA (general)
    if re.search(r'(RNA|lnc|NCRNA)', name, re.IGNORECASE) and not re.match(r'^[A-Z0-9]{2,6}$', name):
        return { 
            
            'type': 'Long Non-Coding RNA (lncRNA)',
            'description': 'Non-protein coding RNA transcript',
            'example': 'MALAT1, NEAT1',
            'pattern': 'Various, often contains RNA/lnc in name'
        }
    
    # RefSeq identifiers
    if re.match(r'^(NM|NR|XM|XR)_\d+$', name):
        return { 
            
            'type': 'RefSeq Identifier',
            'description': 'NCBI Reference Sequence database identifier',
            'example': 'NM_000546 (TP53)',
            'pattern': 'NM/NR/XM/XR followed by underscore and digits'
        }
    
    # Standard HGNC protein-coding gene
    if re.match(r'^[A-Z][A-Z0-9]{1,10}$', name):
        return { 
            
            'type': 'HGNC Protein-Coding Gene Symbol',
            'description': 'Standard HUGO Gene Nomenclature Committee approved symbol',
            'example': 'TP53, BRCA1, C1QC',
            'pattern': 'Uppercase letters and numbers, typically 2-6 characters'
        }
    
    # Mitochondrial genes
    if re.match(r'^MT-', name, re.IGNORECASE):
        return { 
            
            'type': 'Mitochondrial Gene',
            'description': 'Gene encoded in mitochondrial genome',
            'example': 'MT-CO1, MT-ND1',
            'pattern': 'MT- prefix followed by gene symbol'
        }
    
    # Unknown/Other
    return { 
        
        'type': 'Unknown/Other',
        'description': 'Does not match standard naming conventions',
        'example': 'N/A',
        'pattern': 'Custom or non-standard format'
    }